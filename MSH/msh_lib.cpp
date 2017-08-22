/**************************************************************************
MSHLib - Object:
Task:
Programing:
08/2005 OK Encapsulated from mshlib
**************************************************************************/

#include "math.h"
// C++
#include <string>
#include <vector>

// FileIO
#include "GMSHInterface.h"

// GEOLib
#include "geo_lib.h"
#include "files0.h"
// MSHLib
#include "msh_lib.h"
// PCSLib
#include "mathlib.h"
#include "rf_mmp_new.h"                           //OK411

#ifdef RFW_FRACTURE
#include"rf_mmp_new.h"
#endif
#ifdef USE_TOKENBUF
#include "tokenbuf.h"
#endif

// WW extern void RFConfigRenumber(void);
#ifndef NEW_EQS                                   //WW. 07.11.2008
extern void ConfigRenumberProperties(void);
#endif
extern int ReadRFIFile(std::string g_strFileNameBase);
#include "rf_pcs.h"
#include "gs_project.h"

std::vector<Mesh_Group::CFEMesh*>fem_msh_vector;

#define FEM_FILE_EXTENSION ".msh"

double msh_x_min,msh_x_max;                       //OK
double msh_y_min,msh_y_max;                       //OK
double msh_z_min,msh_z_max;                       //OK
double msh_x_mid,msh_y_mid,msh_z_mid;             //OK

#define MSH_SIZE 1e5

using namespace Math_Group;

/**************************************************************************
FEMLib-Method:
Task:
Programing:
04/2005 OK Implementation
**************************************************************************/
void MSHDelete(std::string m_msh_name)
{
   CFEMesh* m_fem_msh = NULL;
   size_t fem_msh_vector_size = fem_msh_vector.size();
   for(size_t i=0;i<fem_msh_vector_size;i++)
   {
      m_fem_msh = fem_msh_vector[i];
      if(m_fem_msh->pcs_name.compare(m_msh_name)==0)
      {
         delete m_fem_msh;
         fem_msh_vector.erase((fem_msh_vector.begin()+i));
      }
   }
}


/**************************************************************************
FEMLib-Method:
03/2005 OK Implementation
05/2005 TK modified
05/2006 TK modified
**************************************************************************/
void FEMDeleteAll()
{
   for(int i=0;i<(int)fem_msh_vector.size();i++)
   {
      delete fem_msh_vector[i];
      fem_msh_vector[i]=NULL;
   }
   fem_msh_vector.clear();
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2005 OK Implementation
08/2005 WW Topology construction and rfi compatible
10/2005 OK BINARY
08/2010 KR deleted binary mesh read
**************************************************************************/

CFEMesh* FEMRead(const std::string &file_base_name, GEOLIB::GEOObjects* geo_obj, std::string* unique_name)
{
   //----------------------------------------------------------------------
   CFEMesh *fem_msh (NULL);
   char line[MAX_ZEILE];
   std::string sub_line;
   std::string line_string;
   std::ios::pos_type position;
   //========================================================================
   // File handling
   std::string msh_file_name_ascii = file_base_name + FEM_FILE_EXTENSION;

   // test if this is a GMSH mesh
   if (FileIO::GMSHInterface::isGMSHMeshFile (msh_file_name_ascii))
   {
      // create mesh object
      fem_msh = new CFEMesh();
      // read mesh
      GMSH2MSH(msh_file_name_ascii.c_str(), fem_msh);
      return fem_msh;
   }

   std::ifstream msh_file_ascii;

#ifdef USE_TOKENBUF
   TokenBuf *tokenbuf;
#endif
   std::cout << "MSHRead: ";
   std::cout << "ASCII file" << std::endl;
   msh_file_ascii.open(msh_file_name_ascii.data(),std::ios::in);
   if (!msh_file_ascii.good()) return NULL;

   //----------------------------------------------------------------------
   // RFI - WW
   bool rfiMesh = true;
   getline(msh_file_ascii, line_string);          // The first line
   if(line_string.find("#FEM_MSH")!=std::string::npos)
      rfiMesh = false;
   //	else	// KR: included error message
   //	{
   //		std::cout << "Error in CFEMesh::FEMRead() - The file \"" << file_base_name << "\" is not an OpenGeoSys mesh file." << std::endl;
   //		return NULL;
   //	}

   msh_file_ascii.seekg(0L,std::ios::beg);

#ifdef USE_TOKENBUF
   tokenbuf = new TokenBuf(msh_file_ascii, 10485760);
#endif

   if (rfiMesh)
   {
      fem_msh = new CFEMesh(geo_obj, unique_name);
      Read_RFI(msh_file_ascii, fem_msh);
      //KR fem_msh_vector.push_back(fem_msh);
      msh_file_ascii.close();
      return fem_msh;
   }

#ifdef USE_TOKENBUF
   while(!tokenbuf->done())
   {
      tokenbuf->get_non_empty_line(line, MAX_ZEILE);
      line_string = line;
      if(line_string.find("#STOP") != std::string::npos)
         return true;

                                                  // mesh
      if(line_string.find("#FEM_MSH") != string::npos)
      {
         m_fem_msh = new CFEMesh();
         m_fem_msh->Read(tokenbuf);
         fem_msh_vector.push_back(m_fem_msh);
      }
   }
   delete tokenbuf;
#else

   while(!msh_file_ascii.eof())
   {
      msh_file_ascii.getline(line,MAX_ZEILE);
      line_string = line;
      if(line_string.find("#STOP")!=std::string::npos)
         return fem_msh;
      //..................................................................
                                                  // keyword found
      if(line_string.find("#FEM_MSH")!=std::string::npos)
      {
         fem_msh = new CFEMesh(geo_obj, unique_name);
         position = fem_msh->Read(&msh_file_ascii);
         //fem_msh_vector.push_back(m_fem_msh);
         msh_file_ascii.seekg(position,std::ios::beg);
      }                                           // keyword found
   }                                              // eof
#endif

   msh_file_ascii.close();

   //========================================================================
   return fem_msh;
}


/**************************************************************************
MSHLib-Method: Read rfi file ()
Task:
Programing:
08/2005 WW Re-implememtation
**************************************************************************/
void Read_RFI(std::istream& msh_file,CFEMesh* m_msh)
{
   long id;
   long i=0;
   int NumNodes=0;
   int NumElements=0;
   int End = 1;
   double x,y,z;
   std::string strbuffer;

   CNode* node = NULL;
   CElem* elem = NULL;
   //----------------------------------------------------------------------
   while (End)
   {
      getline(msh_file, strbuffer);               // The first line
      msh_file>>i>>NumNodes>>NumElements>>std::ws;
      //....................................................................
      // Node data
      for(i=0;i<NumNodes;i++)
      {
         msh_file>>id>>x>>y>>z>>std::ws;
         node = new CNode(id,x,y,z);
         m_msh->nod_vector.push_back(node);
      }
      for(i=0;i<NumElements; i++)
      {
         elem = new CElem(i);
         elem->Read(msh_file, 1);
         m_msh->ele_vector.push_back(elem);
      }
      End =0;
   }
}


/**************************************************************************
MSHLib-Method:
02/2006 WW Implementation
**************************************************************************/
void CompleteMesh()
{
   for(int i=0;i<(int)fem_msh_vector.size(); i++)
   {
      fem_msh_vector[i]->ConstructGrid();
      fem_msh_vector[i]->FillTransformMatrix();
   }
}


/**************************************************************************
FEMLib-Method:
Task: Master write functionn
Programing:
03/2005 OK Implementation
10/2005 OK BINARY
last modification:
08/2010	KR binary case deleted
**************************************************************************/
void MSHWrite(std::string file_base_name)
{
   CFEMesh* m_fem_msh = NULL;
   std::string sub_line;
   std::string line_string;
   //----------------------------------------------------------------------
   // File handling
   std::string fem_msh_file_name = file_base_name + FEM_FILE_EXTENSION;
   std::fstream fem_msh_file;
   std::string msh_file_test_name = file_base_name + "_test" + FEM_FILE_EXTENSION;
   std::fstream msh_file_test;
   for(size_t i=0; i<fem_msh_vector.size(); i++)
   {
      m_fem_msh = fem_msh_vector[i];
   }
   fem_msh_file.open(fem_msh_file_name.c_str(),std::ios::trunc|std::ios::out);
   if(!fem_msh_file.good()) return;
   fem_msh_file.setf(std::ios::scientific,std::ios::floatfield);
   fem_msh_file.precision(12);

   //----------------------------------------------------------------------
   for(size_t i=0; i<fem_msh_vector.size(); i++)
   {
      m_fem_msh = fem_msh_vector[i];
      m_fem_msh->Write(&fem_msh_file);
   }
   //----------------------------------------------------------------------
   fem_msh_file << "#STOP";
   fem_msh_file.close();
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2005 OK Implementation
last modification:
**************************************************************************/
CFEMesh* FEMGet(const std::string &msh_name)
{
   size_t no_msh = fem_msh_vector.size();
   // If there is only one msh file available, use it for all process. WW
   if (no_msh==1) return fem_msh_vector[0];       //WW
   for(size_t i=0;i<no_msh;i++)
   {
      if ( fem_msh_vector[i]->pcs_name.compare(msh_name) == 0 )
         return fem_msh_vector[i];
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
09/2004 OK Implementation
**************************************************************************/
void MSHCalcMinMaxMidCoordinates()
{
   double m_dXMin1 = 1.e+19;
   double m_dXMax1 = -1.e+19;
   double m_dYMin1 = 1.e+19;
   double m_dYMax1 = -1.e+19;
   double m_dZMin1 = 1.e+19;
   double m_dZMax1 = -1.e+19;
   double value;
   CFEMesh* m_msh = NULL;
   //----------------------------------------------------------------------
   for(int j=0;j<(int)fem_msh_vector.size();j++)
   {
      m_msh = fem_msh_vector[j];
      for(long i=0;i<(long)m_msh->nod_vector.size();i++)
      {
         value = m_msh->nod_vector[i]->X();
         if(value<m_dXMin1) m_dXMin1 = value;
         if(value>m_dXMax1) m_dXMax1 = value;
         value = m_msh->nod_vector[i]->Y();
         if(value<m_dYMin1) m_dYMin1 = value;
         if(value>m_dYMax1) m_dYMax1 = value;
         value = m_msh->nod_vector[i]->Z();
         if(value<m_dZMin1) m_dZMin1 = value;
         if(value>m_dZMax1) m_dZMax1 = value;
         //..................................................................
         // Shrink a bit
         msh_x_min = m_dXMin1 - 0.05*(m_dXMax1-m_dXMin1);
         msh_x_max = m_dXMax1 + 0.05*(m_dXMax1-m_dXMin1);
         msh_y_min = m_dYMin1 - 0.05*(m_dYMax1-m_dYMin1);
         msh_y_max = m_dYMax1 + 0.05*(m_dYMax1-m_dYMin1);
         msh_z_min = m_dZMin1 - 0.05*(m_dZMax1-m_dZMin1);
         msh_z_max = m_dZMax1 + 0.05*(m_dZMax1-m_dZMin1);
      }
   }
   //----------------------------------------------------------------------
   msh_x_mid = 0.5*(msh_x_min+msh_x_max);
   msh_y_mid = 0.5*(msh_y_min+msh_y_max);
   msh_z_mid = 0.5*(msh_z_min+msh_z_max);
   //----------------------------------------------------------------------
}


/**************************************************************************
GeoLib-Method:
Task:
Programing:
04/2004 OK Implementation
01/2005 OK File handling
09/2005 OK MSH ToDo
last modification:
**************************************************************************/
void MSHWriteVOL2TEC(std::string m_msh_name)
{
   long i,j;
   CGLVolume *m_vol = NULL;
   std::vector<CGLVolume*>::const_iterator p_vol;
   std::string name("VOLUMES");
   std::vector<Surface*>::const_iterator p_sfc;
   double x,y,z;
   //  CGLPoint m_point; // TF
   std::ios::pos_type position;
   int vol_number = -1;
   Surface* m_sfc = NULL;
   //--------------------------------------------------------------------
   CFEMesh* m_msh (FEMGet(m_msh_name));
   if(!m_msh)
      return;
   long no_nodes = (long)m_msh->nod_vector.size();
   long ep_layer = (long)m_msh->ele_vector.size() / m_msh->getNumberOfMeshLayers();
   //--------------------------------------------------------------------
   // File handling
   std::string tec_path;
   CGSProject* m_gsp = GSPGetMember("gli");
   if(m_gsp)
      tec_path = m_gsp->path;
   //======================================================================
   p_vol = volume_vector.begin();
   while(p_vol!=volume_vector.end())
   {
      m_vol = *p_vol;
      if(m_vol->layer==0)                         //OK
      {
         p_vol++;
         continue;
      }
      p_sfc = m_vol->surface_vector.begin();
      m_sfc = *p_sfc;
      if(!m_sfc)
         return;
      //--------------------------------------------------------------------
      long jb = (m_vol->layer-1)*ep_layer;
      long je = jb + ep_layer;
      vol_number++;
      //--------------------------------------------------------------------
      // file handling
      std::string vol_file_name = tec_path + "VOL_" + m_vol->name + TEC_FILE_EXTENSION;
      std::fstream vol_file (vol_file_name.data(),std::ios::trunc|std::ios::out);
      vol_file.setf(std::ios::scientific,std::ios::floatfield);
      vol_file.precision(12);
      if (!vol_file.good()) return;
      vol_file.seekg(0L,std::ios::beg);
      //--------------------------------------------------------------------
      vol_file << "VARIABLES = X,Y,Z,VOL" << std::endl;
      //--------------------------------------------------------------------
      long no_mat_elements = 0;
      CElem* m_ele = NULL;
      vec<long>node_indeces(6);
      for(i=jb;i<je;i++)
      {
         m_ele = m_msh->ele_vector[i];
         if (m_ele->GetElementType()==MshElemType::PRISM)
         {
            m_ele->GetNodeIndeces(node_indeces);
            //nodes = m_msh->ele_vector[i]->nodes;
            x=0.0; y=0.0; z=0.0;
            for(j=0;j<6;j++)
            {
               x += m_msh->nod_vector[node_indeces[j]]->X();
               y += m_msh->nod_vector[node_indeces[j]]->Y();
               z += m_msh->nod_vector[node_indeces[j]]->Z();
            }
            x /= double(6);
            y /= double(6);
            z /= double(6);
            // TF m_sfc->PointInSurface(&m_point) returns always false
            //        m_point.x = x;
            //        m_point.y = y;
            //        m_point.z = z;
            //        if(m_sfc->PointInSurface(&m_point)){
            //          no_mat_elements++;
            //        }
         }
      }
      //--------------------------------------------------------------------
      position = vol_file.tellg();
      vol_file << "ZONE T = " << m_vol->name << ", " \
         << "N = " << no_nodes << ", " \
         << "E = " << no_mat_elements << ", " \
         << "F = FEPOINT" << ", " << "ET = BRICK" << std::endl;
      for(i=0;i<no_nodes;i++)
      {
         vol_file \
            << m_msh->nod_vector[i]->X() << " " << m_msh->nod_vector[i]->Y() << " " << m_msh->nod_vector[i]->Z() << " " << vol_number << std::endl;
      }
      for(long i=jb;i<je;i++)
      {
         m_ele = m_msh->ele_vector[i];
         if(m_ele->GetElementType()==MshElemType::PRISM)
         {
            m_ele->GetNodeIndeces(node_indeces);
            x=0.0; y=0.0; z=0.0;
            for(j=0;j<6;j++)
            {
               x += m_msh->nod_vector[node_indeces[j]]->X();
               y += m_msh->nod_vector[node_indeces[j]]->Y();
               z += m_msh->nod_vector[node_indeces[j]]->Z();
            }
            x /= double(6);
            y /= double(6);
            z /= double(6);
            // TF m_sfc->PointInSurface(&m_point) returns always false
            //        m_point.x = x;
            //        m_point.y = y;
            //        m_point.z = z;
            //        if(m_sfc->PointInSurface(&m_point)){
            //          vol_file
            //            << node_indeces[0]+1 << " " << node_indeces[0]+1 << " " << node_indeces[1]+1 << " " << node_indeces[2]+1 << " "
            //            << node_indeces[3]+1 << " " << node_indeces[3]+1 << " " << node_indeces[4]+1 << " " << node_indeces[5]+1 << std::endl;
            //        }
         }
      }
      ++p_vol;
      //======================================================================
   }
}


/**************************************************************************
GeoLib-Method:
Task:
Programing:
04/2005 OK Implementation
11/2005 OK OO-ELE
**************************************************************************/
void MSHWriteTecplot()
{
   MshElemType::type ele_type = MshElemType::INVALID;
   long no_nodes;
   long no_elements;
   std::string delimiter(", ");
   long i;
   CElem* m_ele = NULL;
   vec<long> node_indeces(8);
   //----------------------------------------------------------------------
   // File handling
   std::string file_path = "MSH";
   CGSProject* m_gsp = NULL;
   m_gsp = GSPGetMember("msh");
   if (m_gsp)
      file_path = m_gsp->path + "MSH";
   //----------------------------------------------------------------------
   CFEMesh* m_msh = NULL;
   for (int j = 0; j < (int) fem_msh_vector.size(); j++)
   {
      m_msh = fem_msh_vector[j];
      no_nodes = (long) m_msh->nod_vector.size();
      no_elements = (long) m_msh->ele_vector.size();
      // Test ele_type
      if (no_elements > 0)
      {
         m_ele = m_msh->ele_vector[0];
         ele_type = m_ele->GetElementType();
      }
      // File handling
      std::string msh_file_name = file_path + "_" + m_msh->pcs_name
         + TEC_FILE_EXTENSION;
      std::fstream msh_file(msh_file_name.data(), std::ios::trunc | std::ios::out);
      msh_file.setf(std::ios::scientific, std::ios::floatfield);
      msh_file.precision(12);
      if (!msh_file.good())
         return;
      msh_file.seekg(0L, std::ios::beg);
      msh_file << "VARIABLES = X,Y,Z" << std::endl;
      msh_file << "ZONE T = " << m_msh->pcs_name << delimiter << "N = "
         << no_nodes << delimiter << "E = " << no_elements << delimiter;
      msh_file << "F = FEPOINT" << delimiter;
      switch (ele_type)
      {
         //..................................................................
         case MshElemType::LINE:
            msh_file << "ET = QUADRILATERAL" << std::endl;
            for (i = 0; i < no_nodes; i++)
            {
               msh_file << m_msh->nod_vector[i]->X() << " "
                  << m_msh->nod_vector[i]->Y() << " "
                  << m_msh->nod_vector[i]->Z() << std::endl;
            }
            for (i = 0; i < no_elements; i++)
            {
               m_ele = m_msh->ele_vector[i];
               m_ele->GetNodeIndeces(node_indeces);
               msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1
                  << " " << node_indeces[1] + 1 << " " << node_indeces[0]
                  + 1 << std::endl;
            }
            break;
            //..................................................................
         case MshElemType::QUAD:
            msh_file << "ET = QUADRILATERAL" << std::endl;
            for (i = 0; i < no_nodes; i++)
            {
               msh_file << m_msh->nod_vector[i]->X() << " "
                  << m_msh->nod_vector[i]->Y() << " "
                  << m_msh->nod_vector[i]->Z() << std::endl;
            }
            for (i = 0; i < no_elements; i++)
            {
               m_ele = m_msh->ele_vector[i];
               m_ele->GetNodeIndeces(node_indeces);
               msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1
                  << " " << node_indeces[2] + 1 << " " << node_indeces[3]
                  + 1 << std::endl;
            }
            break;
            //..................................................................
         case MshElemType::HEXAHEDRON:
            msh_file << "ET = BRICK" << std::endl;
            for (i = 0; i < no_nodes; i++)
            {
               msh_file << m_msh->nod_vector[i]->X() << " "
                  << m_msh->nod_vector[i]->Y() << " "
                  << m_msh->nod_vector[i]->Z() << std::endl;
            }
            for (i = 0; i < no_elements; i++)
            {
               m_ele = m_msh->ele_vector[i];
               m_ele->GetNodeIndeces(node_indeces);
               msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1
                  << " " << node_indeces[2] + 1 << " " << node_indeces[3]
                  + 1 << " " << node_indeces[4] + 1 << " "
                  << node_indeces[5] + 1 << " " << node_indeces[6] + 1
                  << " " << node_indeces[7] + 1 << std::endl;
            }
            break;
            //..................................................................
         case MshElemType::TRIANGLE:
            msh_file << "ET = TRIANGLE" << std::endl;
            for (i = 0; i < no_nodes; i++)
            {
               msh_file << m_msh->nod_vector[i]->X() << " "
                  << m_msh->nod_vector[i]->Y() << " "
                  << m_msh->nod_vector[i]->Z() << std::endl;
            }
            for (i = 0; i < no_elements; i++)
            {
               m_ele = m_msh->ele_vector[i];
               m_ele->GetNodeIndeces(node_indeces);
               msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1
                  << " " << node_indeces[2] + 1 << std::endl;
            }
            break;
            //..................................................................
         case MshElemType::TETRAHEDRON:
            msh_file << "ET = TETRAHEDRON" << std::endl;
            for (i = 0; i < no_nodes; i++)
            {
               msh_file << m_msh->nod_vector[i]->X() << " "
                  << m_msh->nod_vector[i]->Y() << " "
                  << m_msh->nod_vector[i]->Z() << std::endl;
            }
            for (i = 0; i < no_elements; i++)
            {
               m_ele = m_msh->ele_vector[i];
               m_ele->GetNodeIndeces(node_indeces);
               msh_file << node_indeces[0] + 1 << " " << node_indeces[1] + 1
                  << " " << node_indeces[2] + 1 << " " << node_indeces[3]
                  + 1 << std::endl;
            }
            break;
            //..................................................................
         case MshElemType::PRISM:
            msh_file << "ET = BRICK" << std::endl;
            for (i = 0; i < no_nodes; i++)
            {
               msh_file << m_msh->nod_vector[i]->X() << " "
                  << m_msh->nod_vector[i]->Y() << " "
                  << m_msh->nod_vector[i]->Z() << std::endl;
            }
            for (i = 0; i < no_elements; i++)
            {
               m_ele = m_msh->ele_vector[i];
               m_ele->GetNodeIndeces(node_indeces);
               if (m_ele->GetElementType() == MshElemType::PRISM)
               {
                  msh_file << node_indeces[0] + 1 << " " << node_indeces[1]
                     + 1 << " " << node_indeces[2] + 1 << " "
                     << node_indeces[2] + 1 << " " << node_indeces[3]
                     + 1 << " " << node_indeces[4] + 1 << " "
                     << node_indeces[5] + 1 << " " << node_indeces[5]
                     + 1 << std::endl;
               }
               if (m_ele->GetElementType() == MshElemType::HEXAHEDRON)
               {
                  msh_file << node_indeces[0] + 1 << " " << node_indeces[1]
                     + 1 << " " << node_indeces[2] + 1 << " "
                     << node_indeces[3] + 1 << " " << node_indeces[4]
                     + 1 << " " << node_indeces[5] + 1 << " "
                     << node_indeces[6] + 1 << " " << node_indeces[7]
                     + 1 << std::endl;
               }
            }
            break;
         default:
            std::cerr << "MSHWriteTecplot MshElemType not handled" << std::endl;
      }
   }
}


/**************************************************************************
GeoLib-Method:
Task:
Programing:
04/2005 OK Implementation
11/2005 OK OO-ELE
**************************************************************************/
void MSHLayerWriteTecplot()
{
   MshElemType::type ele_type = MshElemType::INVALID;
   long no_nodes;
   long no_elements;
   std::string delimiter(", ");
   CElem* m_ele = NULL;
   vec<long> node_indeces(8);
   std::string no_layer_str;
   char no_layer_char[3];
   //----------------------------------------------------------------------
   // File handling
   std::string file_path = "MSH";
   CGSProject* m_gsp = NULL;
   m_gsp = GSPGetMember("msh");
   if (m_gsp)
      file_path = m_gsp->path;
   //----------------------------------------------------------------------
   CFEMesh* m_msh = NULL;
   for (int j = 0; j < (int) fem_msh_vector.size(); j++)
   {
      m_msh = fem_msh_vector[j];
      for (size_t k = 0; k < m_msh->getNumberOfMeshLayers(); k++)
      {
         sprintf(no_layer_char, "%lu", k + 1);
         no_layer_str = no_layer_char;
         no_nodes = (long) m_msh->nod_vector.size() / (m_msh->getNumberOfMeshLayers()
            + 1);
         no_elements = (long) m_msh->ele_vector.size() / m_msh->getNumberOfMeshLayers();
         // Test ele_type
         if (no_elements > 0)
         {
            m_ele = m_msh->ele_vector[0];
            ele_type = m_ele->GetElementType();
         }
         // File handling
         std::string msh_file_name = file_path + "MSH_LAYER" + no_layer_str + "_"
            + m_msh->pcs_name + TEC_FILE_EXTENSION;
         std::fstream msh_file(msh_file_name.data(), std::ios::trunc | std::ios::out);
         msh_file.setf(std::ios::scientific, std::ios::floatfield);
         msh_file.precision(12);
         if (!msh_file.good())
            return;
         msh_file.seekg(0L, std::ios::beg);
         msh_file << "VARIABLES = X,Y,Z" << std::endl;
         msh_file << "ZONE T = " << m_msh->pcs_name << delimiter << "N = "
            << (long) m_msh->nod_vector.size() << delimiter << "E = "
            << no_elements << delimiter;
         msh_file << "F = FEPOINT" << delimiter;
         switch (ele_type)
         {
            //..................................................................
            case MshElemType::LINE:
               msh_file << "ET = QUADRILATERAL" << std::endl;
               for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
               {
                  msh_file << m_msh->nod_vector[i]->X() << " "
                     << m_msh->nod_vector[i]->Y() << " "
                     << m_msh->nod_vector[i]->Z() << std::endl;
               }
               for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
               {
                  m_ele = m_msh->ele_vector[i];
                  m_ele->GetNodeIndeces(node_indeces);
                  msh_file << node_indeces[0] + 1 << " " << node_indeces[1]
                     + 1 << " " << node_indeces[1] + 1 << " "
                     << node_indeces[0] + 1 << std::endl;
               }
               break;
               //..................................................................
            case MshElemType::QUAD:
               msh_file << "ET = QUADRILATERAL" << std::endl;
               for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
               {
                  msh_file << m_msh->nod_vector[i]->X() << " "
                     << m_msh->nod_vector[i]->Y() << " "
                     << m_msh->nod_vector[i]->Z() << std::endl;
               }
               for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
               {
                  m_ele = m_msh->ele_vector[i];
                  m_ele->GetNodeIndeces(node_indeces);
                  msh_file << node_indeces[0] + 1 << " " << node_indeces[1]
                     + 1 << " " << node_indeces[2] + 1 << " "
                     << node_indeces[3] + 1 << std::endl;
               }
               break;
               //..................................................................
            case MshElemType::HEXAHEDRON:
               msh_file << "ET = BRICK" << std::endl;
               for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
               {
                  msh_file << m_msh->nod_vector[i]->X() << " "
                     << m_msh->nod_vector[i]->Y() << " "
                     << m_msh->nod_vector[i]->Z() << std::endl;
               }
               for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
               {
                  m_ele = m_msh->ele_vector[i];
                  m_ele->GetNodeIndeces(node_indeces);
                  msh_file << node_indeces[0] + 1 << " " << node_indeces[1]
                     + 1 << " " << node_indeces[2] + 1 << " "
                     << node_indeces[3] + 1 << " " << node_indeces[4]
                     + 1 << " " << node_indeces[5] + 1 << " "
                     << node_indeces[6] + 1 << " " << node_indeces[7]
                     + 1 << std::endl;
               }
               break;
               //..................................................................
            case MshElemType::TRIANGLE:
               msh_file << "ET = TRIANGLE" << std::endl;
               for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
               {
                  msh_file << m_msh->nod_vector[i]->X() << " "
                     << m_msh->nod_vector[i]->Y() << " "
                     << m_msh->nod_vector[i]->Z() << std::endl;
               }
               for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
               {
                  m_ele = m_msh->ele_vector[i];
                  m_ele->GetNodeIndeces(node_indeces);
                  msh_file << node_indeces[0] + 1 << " " << node_indeces[1]
                     + 1 << " " << node_indeces[2] + 1 << std::endl;
               }
               break;
               //..................................................................
            case MshElemType::TETRAHEDRON:
               msh_file << "ET = TETRAHEDRON" << std::endl;
               for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
               {
                  msh_file << m_msh->nod_vector[i]->X() << " "
                     << m_msh->nod_vector[i]->Y() << " "
                     << m_msh->nod_vector[i]->Z() << std::endl;
               }
               for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
               {
                  m_ele = m_msh->ele_vector[i];
                  m_ele->GetNodeIndeces(node_indeces);
                  msh_file << node_indeces[0] + 1 << " " << node_indeces[1]
                     + 1 << " " << node_indeces[2] + 1 << " "
                     << node_indeces[3] + 1 << std::endl;
               }
               break;
               //..................................................................
            case MshElemType::PRISM:
               msh_file << "ET = BRICK" << std::endl;
               for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
               {
                  msh_file << m_msh->nod_vector[i]->X() << " "
                     << m_msh->nod_vector[i]->Y() << " "
                     << m_msh->nod_vector[i]->Z() << std::endl;
               }
               for (size_t i = k * no_elements; i < (k + 1) * no_elements; i++)
               {
                  m_ele = m_msh->ele_vector[i];
                  m_ele->GetNodeIndeces(node_indeces);
                  msh_file << node_indeces[0] + 1 << " " << node_indeces[1]
                     + 1 << " " << node_indeces[2] + 1 << " "
                     << node_indeces[2] + 1 << " " << node_indeces[3]
                     + 1 << " " << node_indeces[4] + 1 << " "
                     << node_indeces[5] + 1 << " " << node_indeces[5]
                     + 1 << std::endl;
               }
               break;
            default:
               std::cerr << "MSHLayerWriteTecplot MshElemType not handled" << std::endl;
         }
      }                                           // layer
   }
}


/**************************************************************************
MSHLib-Method:
12/2005 OK Implementation
07/2007 OK PCS
**************************************************************************/
CFEMesh* MSHGet(const std::string &geo_name)
{
   CFEMesh* m_msh = NULL;
   for(int i=0;i<(int)fem_msh_vector.size();i++)
   {
      m_msh = fem_msh_vector[i];
      if(m_msh->geo_name.compare(geo_name)==0)
      {
         return m_msh;
      }
      if(m_msh->pcs_name.compare(geo_name)==0)
      {
         return m_msh;
      }
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
12/2005 OK Implementation
**************************************************************************/
CFEMesh* MSHGet(const std::string &pcs_type_name,const std::string &geo_name)
{
   CFEMesh* m_msh = NULL;
   for(int i=0;i<(int)fem_msh_vector.size();i++)
   {
      m_msh = fem_msh_vector[i];
      if((m_msh->pcs_name.compare(pcs_type_name)==0)&&\
         (m_msh->geo_name.compare(geo_name)==0)) \
      { \
         return m_msh;
      }
   }
   return NULL;
}


/**************************************************************************
PCSLib-Method:
12/2005 OK Implementation
**************************************************************************/
CFEMesh* MSHGetGEO(std::string geo_name)
{
   int no_msh = (int)fem_msh_vector.size();
   // If there is only one msh file available, use it for all process. WW
   if(no_msh==1)
      return fem_msh_vector[0];                   //WW
   //----------------------------------------------------------------------
   CFEMesh* m_msh = NULL;
   for(int i=0;i<no_msh;i++)
   {
      m_msh = fem_msh_vector[i];
      if(m_msh->geo_name.compare(geo_name)==0)
         return m_msh;
   }
   //----------------------------------------------------------------------
   return NULL;
}


#ifdef RFW_FRACTURE
/**************************************************************************************
 ROCKFLOW - Funktion: ELEGetCommonNodes
 Aufgabe:
   Get the nodes shared by 2 elements
 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: const CElem* elem1 and elem2  :  element pointers
   R: vec<CNode*> nodes  :  a vector containing pointers to the matching (i.e. common) nodes
 Ergebnis:
   Returns the value true if elements are actually neighbors(i.e. have nodes in common),
   otherwise returns false
 Programmaenderungen:
05/2005 RFW Implementierung
05/2006 RFW �nderung
***************************************************************************************/
bool MSHGetCommonNodes(CElem* elem1, CElem* elem2, vector<CNode*>& nodes)
{
   nodes.clear();
   bool neighbor = false;
   vec<CNode*> nodes1, nodes2;
   int numnodes1, numnodes2;
   if(elem1!=NULL && elem2!=NULL)                 //RFW 19/05/2005
   {
      elem1->GetNodes(nodes1);
      elem2->GetNodes(nodes2);
      numnodes1 = elem1->GetVertexNumber();
      numnodes2 = elem2->GetVertexNumber();

      for(int i=0; i!=numnodes1; ++i)
      {
         for(int j=0; j!=numnodes2; ++j)
         {
            if(nodes1[i] == nodes2[j])
               nodes.push_back(nodes1[i]);
         }
      }
      if(nodes.size() == 0)
         return neighbor;
      else
      {
         neighbor = true;
         return neighbor;
      }
   }
   else
      return neighbor;
}


/*************************************************************************
 ROCKFLOW - Funktion: MSHSetFractureElements
 Aufgabe:
   This function gives each fracture element a weight and an aperture direction.  It
   does this by using polylines that define the top and bottom of the 2D fracture.  For
   the function to work properly, the points defining the polyline must be evenly spaced.
   Function should only be called once, after configuration.
 Formalparameter:
 Ergebnis:
 Programmaenderungen:
   11/2005 RFW Implementierung
***************************************************************************/
void MSHSetFractureElements(void)
{
   //** TEST *************
   std::cout << "Location 0 MSHSetFractureElements\n";
   //** TEST *************
   int group, frac_exists=0;
   CElem *elem = NULL;
   CGLPolyline *frac_top = NULL, *frac_bot = NULL;
   CMediumProperties *m_mmp = NULL;
   CFEMesh* m_msh = NULL;
   group = -1;                                    //WW
   vector<long> bound_elements;
   vector<double> node_x, node_y;
   double x, y, tri[6], dx_avg, dy_avg, d_max, *gravity_center;
   //The following vector is a strange beast.  The first index indicates the fracture being considered, the second
   //index indicates the segment of the fracture being considered, and the final index contains the element numbers
   //of those elements making up the given fracture segment.
   vector<vector<vector<int> > > segment_elements;
   vector<vector<bool> > segment_on_boundary;
   vector<vector<double> > segment_dx;
   vector<vector<double> > segment_dy;
   vector<int> fractures;

   //finding the material group for the fractures, there should only be 1
   for(int i=0; i<(int)mmp_vector.size(); ++i)
   {
      m_mmp = mmp_vector[i];
      //for(int j=0; j<(int)m_mmp->relative_permeability_function.size(); ++j)
      //{
      if (m_mmp->frac_num != 0)
      {
         group = i;
         frac_exists = 1;
      }
      //}
   }

   if(frac_exists)
   {
      m_mmp = mmp_vector[group];
      //frac_segments.resize(m_mmp->frac_num);
      segment_elements.resize(m_mmp->frac_num);
      segment_on_boundary.resize(m_mmp->frac_num);
      segment_dx.resize(m_mmp->frac_num);
      segment_dy.resize(m_mmp->frac_num);

      //grabbing the nodes of the boundary polyline, will be used later---------------------START
      CGLPolyline *bound_ply = NULL;
      vector<long> poly_nodes;

      bound_ply = GEOGetPLYByName("BOUNDARY");
      m_msh = fem_msh_vector[0];                  //this isn't great

      if(bound_ply)
      {
         m_msh->GetNODOnPLY(bound_ply, poly_nodes);
      }
      else
      {
         std::cout <<"Error 1 in ELESetBoundaryElements, no BOUNDARY polyline.\n";
         abort();
      }
      //grabbing the nodes of the boundary polyline-------------------------------------------------END

      for(long j=0; j<m_mmp->frac_num; ++j)       //loop1, over all fractures
      {
         std::string polyname_top = m_mmp->frac_names[j] + "_top";
         std::string polyname_bot = m_mmp->frac_names[j] + "_bot";
         frac_top = GEOGetPLYByName(polyname_top);
         frac_bot = GEOGetPLYByName(polyname_bot);
         //frac_segments[j].resize((long)frac_top->point_vector.size()-1);
         segment_elements[j].resize((long)frac_top->point_vector.size()-1);
         segment_on_boundary[j].resize((long)frac_top->point_vector.size()-1);
         segment_dx[j].resize((long)frac_top->point_vector.size()-1);
         segment_dy[j].resize((long)frac_top->point_vector.size()-1);

         if(!frac_top || !frac_bot)
         {
            std::cout << "Error 1 in MSHSetFractureElements: fracture polyline "<<m_mmp->frac_names[j]<<" not found.\n";
            abort();
         }
                                                  //loop2, over fracture segments
         for(long k=0; k<((long)frac_top->point_vector.size()-1); ++k)
         {
            node_x.resize(4);                     //of course, these are not really nodes, but I'm using code from ELEWhatElemIsPointIn
            node_y.resize(4);
            //4 points defining the fracture segment
            node_x[0]=frac_bot->point_vector[k]->x;        node_y[0]=frac_bot->point_vector[k]->y;
            node_x[1]=frac_bot->point_vector[k+1]->x;    node_y[1]=frac_bot->point_vector[k+1]->y;
            node_x[2]=frac_top->point_vector[k+1]->x;    node_y[2]=frac_top->point_vector[k+1]->y;
            node_x[3]=frac_top->point_vector[k]->x;        node_y[3]=frac_top->point_vector[k]->y;

            //2 trangular subareas defining the fracture segment
            double Area1, Area2;
            Area1 = ( (node_x[1]*node_y[2] - node_x[2]*node_y[1])
               + (node_x[2]*node_y[0] - node_x[0]*node_y[2])
               + (node_x[0]*node_y[1] - node_x[1]*node_y[0]) )/2;
            Area2 = ( (node_x[3]*node_y[2] - node_x[2]*node_y[3])
               + (node_x[2]*node_y[0] - node_x[0]*node_y[2])
               + (node_x[0]*node_y[3] - node_x[3]*node_y[0]) )/2;

                                                  //loop 3, all elements
            for (int l = 0; l < (int)m_msh->ele_vector.size(); l++)
            {
               elem = m_msh->ele_vector[l];
               if(elem->GetPatchIndex()==group)   //this is simply to make things go faster
               {

                  gravity_center = elem->GetGravityCenter();
                  x = gravity_center[0];          //x = ELEGetEleMidPoint(l,0);
                  y = gravity_center[1];          //y = ELEGetEleMidPoint(l,1);

                  //calculate triangular coordinates for both triangles making up fracture segment
                  //first triangle
                  tri[0] = ( (node_x[1]*node_y[2]-node_x[2]*node_y[1]) + (node_y[1]-node_y[2])*x + (node_x[2]-node_x[1])*y )/(2*Area1);
                  tri[1] = ( (node_x[2]*node_y[0]-node_x[0]*node_y[2])+ (node_y[2]-node_y[0])*x + (node_x[0]-node_x[2])*y )/(2*Area1);
                  tri[2] = ( (node_x[0]*node_y[1]-node_x[1]*node_y[0]) + (node_y[0]-node_y[1])*x + (node_x[1]-node_x[0])*y )/(2*Area1);
                  //second triangle
                  tri[3] = ( (node_x[3]*node_y[2]-node_x[2]*node_y[3]) + (node_y[3]-node_y[2])*x + (node_x[2]-node_x[3])*y )/(2*Area2);
                  tri[4] = ( (node_x[2]*node_y[0]-node_x[0]*node_y[2])+ (node_y[2]-node_y[0])*x + (node_x[0]-node_x[2])*y )/(2*Area2);
                  tri[5] = ( (node_x[0]*node_y[3]-node_x[3]*node_y[0]) + (node_y[0]-node_y[3])*x + (node_x[3]-node_x[0])*y )/(2*Area2);

                  if((tri[0]>=-0.00000001 && tri[1]>=-0.00000001 && tri[2]>=-0.00000001)
                     || (tri[3]>=-0.00000001 && tri[4]>=-0.00000001 && tri[5]>=-0.00000001))
                  {
                     elem->SetFracNum(j);

                     //what is the appropriate search dirction for aperture searches?
                     dx_avg = ( (node_x[1]-node_x[0]) + (node_x[2]-node_x[3]) )/2;
                     dy_avg = ( (node_y[1]-node_y[0]) + (node_y[2]-node_y[3]) )/2;
                     d_max = max(abs(dx_avg), abs(dy_avg));
                     elem->SetFracDx(-1*dy_avg/d_max);
                     elem->SetFracDy(dx_avg/d_max);

                     //increment the number of elements in the fracture segment, for weight calculations
                     //frac_segments[j][k] += 1;
                     segment_elements[j][k].push_back( l );

                     //check if element is on boundary, if so, the segment_on_boundary gets a value of true, otherwise false
                     vector<CNode*> match_nodes;
                     vec<CNode*> elem_nodes;
                     elem->GetNodes(elem_nodes);
                     int elem_num_nodes = elem->GetVertexNumber();

                     for(long m=0; m!=(long)poly_nodes.size(); ++m)
                     {
                        for(int n=0; n!=elem_num_nodes; ++n)
                        {
                                                  //NOT SURE THIS WILL WORK!, probably will
                           if(elem_nodes[n] == m_msh->nod_vector[poly_nodes[m]])
                           {
                              match_nodes.push_back(elem_nodes[n]);
                           }
                        }
                     }
                     if(match_nodes.size()==2)
                     {
                        segment_on_boundary[j][k] = true;
                        dx_avg = match_nodes[0]->X() - match_nodes[1]->X();
                        dy_avg = match_nodes[0]->Y() - match_nodes[1]->Y();
                        d_max = max(abs(dx_avg), abs(dy_avg));
                                                  //no longer the inverese of the average slope as above, rather the orientation of the boundary
                        segment_dx[j][k] = dx_avg/d_max;
                        segment_dy[j][k] = dy_avg/d_max;
                     }

                  }
               }                                  //if GroupNumber
            }                                     //loop 3, all elements
         }                                        //loop2, over fracture segments
      }                                           //loop1, over all fractures

      for(long j=0; j<m_mmp->frac_num; ++j)       //loop4, over all fractures
      {
         double Weight=0, seg_length=0;
         vector<double> point_x, point_y;
         point_x.resize(2); point_y.resize(2);
         std::string polyname_top = m_mmp->frac_names[j] + "_top";
         std::string polyname_bot = m_mmp->frac_names[j] + "_bot";
         frac_top = GEOGetPLYByName(polyname_top);
         frac_bot = GEOGetPLYByName(polyname_bot);

         //calculating polyline length
         double top_poly_length = frac_top->CalcPolylineLength();

                                                  //loop5, over fracture segments
         for(long k=0; k<((long)frac_top->point_vector.size()-1); ++k)
         {
            //calculating segment length
            point_x[0]=frac_top->point_vector[k+1]->x;    point_y[0]=frac_top->point_vector[k+1]->y;
            point_x[1]=frac_top->point_vector[k]->x;        point_y[1]=frac_top->point_vector[k]->y;
            seg_length = sqrt(   pow( (point_x[1]-point_x[0]), 2 ) + pow( (point_y[1]-point_y[0]), 2 )   );

                                                  //loop6, over elements in segment
            for(long l=0; l<(long)segment_elements[j][k].size(); ++l)
            {
               //the weight that each element gets is the total weight of 1, divided by the number of segments
               //in which the fracture is divided, and further divided by the number of elements in the given segment
               //WW                long test1 = ((long)frac_top->point_vector.size()-1);
               //WW                long test2 = (long)segment_elements[j][k].size();
               Weight = seg_length / top_poly_length / (double)segment_elements[j][k].size();

               elem = m_msh->ele_vector[segment_elements[j][k][l]];
               //elem->in_frac = true; // RFW 18/11/2005
               elem->SetFrac(Weight);
               //if element is part of boundary segment, reassign it's search directions appropriately
               if(segment_on_boundary[j][k] == true)
               {
                  elem->SetFracDx(segment_dx[j][k]);
                  elem->SetFracDy(segment_dy[j][k]);
               }

            }                                     //loop6, over elements in segment
         }                                        //loop5, over fracture segments
      }                                           //loop4, over all fractures
   }                                              //if frac_exists

}


/*************************************************************************
 ROCKFLOW - Funktion: MSHResetFractureElements
 Aufgabe:
  Resets the aperture calculation for the next time, this means that CalculateFracAperture
  function will not be called multiple times each timestep
 Formalparameter:
 Ergebnis:
 Programmaenderungen:
   11/2005 RFW Implementierung
***************************************************************************/
void MSHResetFractureElements(void)
{
   CElem *elem = NULL;
   CFEMesh* m_msh = NULL;
   m_msh = fem_msh_vector[0];                     //this isn't great

                                                  //loop 3, all elements
   for (int i = 0; i < (int)m_msh->ele_vector.size(); i++)
   {
      elem = m_msh->ele_vector[i];
      elem->ApertureIsNotSet();
      elem->PermeabilityIsNotSet();
   }

   for(int i = 0; i<(int)mmp_vector.size(); i++)
   {
      mmp_vector[i]->fracs_set=0;
   }

}


/**************************************************************************
GeoSys-FEM-Method:
Task: Find what element a set of given coords is in.  Works for triangles and quadrilatera in a 2d system.
Programming:
      RFW 04.2005 index is an initial guess of an element close to
      the one containing the given point.
**************************************************************************/
long MSHWhatElemIsPointIn(double x, double y, long index)
{
   CElem* elem = NULL;                            //PointElementNow
   vec<CElem*> neighbors(10);
   CFEMesh* m_msh = NULL;
   m_msh = fem_msh_vector[0];                     // Something must be done later on here.
   //need a pointer for node
   bool inside = false, in_domain = false;
   int num_face, num_nodes;
   double Area1, Area2, tri[6];
   vector<double> node_x, node_y;
   long return_value=0, count=0;                  //perhaps return_value shouldn't be initialized
   vector<long> index_guess;
   index_guess.push_back(index);
   long number_of_elem = (long)m_msh->ele_vector.size();

   //----------------------------------------------point in model domain?-------------------------------------------------------START
   //boundary polyline defining the model boundaries, used to check if point is inside or outside model domain
   //this is not a very elegant solution, needs to be improved.  Current assumption, quadratisch model domain.
   CGLPolyline *bound_ply = NULL;
   vector<double> bnode_x, bnode_y;
   bnode_x.resize(4);  bnode_y.resize(4);         //of course, these are not really nodes, rather points
   bound_ply = GEOGetPLYByName("BOUNDARY");
   bnode_x[0]=bound_ply->point_vector[0]->x;       bnode_y[0]=bound_ply->point_vector[0]->y;
   bnode_x[1]=bound_ply->point_vector[1]->x;        bnode_y[1]=bound_ply->point_vector[1]->y;
   bnode_x[2]=bound_ply->point_vector[2]->x;        bnode_y[2]=bound_ply->point_vector[2]->y;
   bnode_x[3]=bound_ply->point_vector[3]->x;        bnode_y[3]=bound_ply->point_vector[3]->y;

   Area1 = ( (bnode_x[1]*bnode_y[2] - bnode_x[2]*bnode_y[1])
      + (bnode_x[2]*bnode_y[0] - bnode_x[0]*bnode_y[2])
      + (bnode_x[0]*bnode_y[1] - bnode_x[1]*bnode_y[0]) )/2;
   Area2 = ( (bnode_x[3]*bnode_y[2] - bnode_x[2]*bnode_y[3])
      + (bnode_x[2]*bnode_y[0] - bnode_x[0]*bnode_y[2])
      + (bnode_x[0]*bnode_y[3] - bnode_x[3]*bnode_y[0]) )/2;
   //calculate triangular coordinates for both triangles making up the element
   //first triangle
   tri[0] = ( (bnode_x[1]*bnode_y[2]-bnode_x[2]*bnode_y[1]) + (bnode_y[1]-bnode_y[2])*x + (bnode_x[2]-bnode_x[1])*y )/(2*Area1);
   tri[1] = ( (bnode_x[2]*bnode_y[0]-bnode_x[0]*bnode_y[2])+ (bnode_y[2]-bnode_y[0])*x + (bnode_x[0]-bnode_x[2])*y )/(2*Area1);
   tri[2] = ( (bnode_x[0]*bnode_y[1]-bnode_x[1]*bnode_y[0]) + (bnode_y[0]-bnode_y[1])*x + (bnode_x[1]-bnode_x[0])*y )/(2*Area1);
   //second triangle
   tri[3] = ( (bnode_x[3]*bnode_y[2]-bnode_x[2]*bnode_y[3]) + (bnode_y[3]-bnode_y[2])*x + (bnode_x[2]-bnode_x[3])*y )/(2*Area2);
   tri[4] = ( (bnode_x[2]*bnode_y[0]-bnode_x[0]*bnode_y[2])+ (bnode_y[2]-bnode_y[0])*x + (bnode_x[0]-bnode_x[2])*y )/(2*Area2);
   tri[5] = ( (bnode_x[0]*bnode_y[3]-bnode_x[3]*bnode_y[0]) + (bnode_y[0]-bnode_y[3])*x + (bnode_x[3]-bnode_x[0])*y )/(2*Area2);

   if((tri[0]>=-0.0000001 && tri[1]>=-0.0000001 && tri[2]>=-0.0000001)
                                                  //perhaps this should be +ve??
      || (tri[3]>=-0.0000001 && tri[4]>=-0.0000001 && tri[5]>=-0.0000001))
   {
      in_domain = true;
   }
   else                                           //point is not in model domain
   {
      return_value = -100;
   }

   //----------------------------------------------point in model domain?---------------------------------------------------------END

   //----------------------------------------------which element is point in?----------------------------------------------------START
   while (!inside && count < (number_of_elem+2) && in_domain)
   {
      count++;
      CNode* node;
                                                  //check to avoid infinite loops
      if ((long)index_guess.size() > number_of_elem)
      {
         std::cout<<"Error in ELEWhatElemIsPointIn, point not in any element";
         return_value = -100;
         break;
      }
      for (long i=0; i!=(long)index_guess.size(); ++i)
      {
         elem = m_msh->ele_vector[index_guess[i]];
         num_nodes = elem->GetVertexNumber();
         node_x.resize(num_nodes);
         node_y.resize(num_nodes);
         for(int j=0; j!=num_nodes; ++j)
         {
            node = elem->GetNode(j);
            node_x[j] = node->X_displaced();
            node_y[j] = node->Y_displaced();
         }
         //----------------------------------------------triangles------------------------------------------------------------
         if(num_nodes==3)
         {
            //calculate area of triangle
            Area1 = ( (node_x[1]*node_y[2] - node_x[2]*node_y[1])
               + (node_x[2]*node_y[0] - node_x[0]*node_y[2])
               + (node_x[0]*node_y[1] - node_x[1]*node_y[0]) )/2;
            //calculate triangular coordinates
            tri[0] = ( (node_x[1]*node_y[2]-node_x[2]*node_y[1]) + (node_y[1]-node_y[2])*x + (node_x[2]-node_x[1])*y )/(2*Area1);
            tri[1] = ( (node_x[2]*node_y[0]-node_x[0]*node_y[2])+ (node_y[2]-node_y[0])*x + (node_x[0]-node_x[2])*y )/(2*Area1);
            tri[2] = ( (node_x[0]*node_y[1]-node_x[1]*node_y[0]) + (node_y[0]-node_y[1])*x + (node_x[1]-node_x[0])*y )/(2*Area1);

                                                  //perhaps this should be +ve??
            if(tri[0]>=-0.0000001 && tri[1]>=-0.0000001 && tri[2]>=-0.0000001 )
            {
               inside = true;
               return_value = index_guess[i];
               break;
            }
         }
         //----------------------------------------------quadrilaterals------------------------------------------------------------
         else if(num_nodes == 4)
         {
            Area1 = ( (node_x[1]*node_y[2] - node_x[2]*node_y[1])
               + (node_x[2]*node_y[0] - node_x[0]*node_y[2])
               + (node_x[0]*node_y[1] - node_x[1]*node_y[0]) )/2;
            Area2 = ( (node_x[3]*node_y[2] - node_x[2]*node_y[3])
               + (node_x[2]*node_y[0] - node_x[0]*node_y[2])
               + (node_x[0]*node_y[3] - node_x[3]*node_y[0]) )/2;
            //calculate triangular coordinates for both triangles making up the element
            //first triangle
            tri[0] = ( (node_x[1]*node_y[2]-node_x[2]*node_y[1]) + (node_y[1]-node_y[2])*x + (node_x[2]-node_x[1])*y )/(2*Area1);
            tri[1] = ( (node_x[2]*node_y[0]-node_x[0]*node_y[2])+ (node_y[2]-node_y[0])*x + (node_x[0]-node_x[2])*y )/(2*Area1);
            tri[2] = ( (node_x[0]*node_y[1]-node_x[1]*node_y[0]) + (node_y[0]-node_y[1])*x + (node_x[1]-node_x[0])*y )/(2*Area1);
            //second triangle
            tri[3] = ( (node_x[3]*node_y[2]-node_x[2]*node_y[3]) + (node_y[3]-node_y[2])*x + (node_x[2]-node_x[3])*y )/(2*Area2);
            tri[4] = ( (node_x[2]*node_y[0]-node_x[0]*node_y[2])+ (node_y[2]-node_y[0])*x + (node_x[0]-node_x[2])*y )/(2*Area2);
            tri[5] = ( (node_x[0]*node_y[3]-node_x[3]*node_y[0]) + (node_y[0]-node_y[3])*x + (node_x[3]-node_x[0])*y )/(2*Area2);

            if((tri[0]>=-0.0000001 && tri[1]>=-0.0000001 && tri[2]>=-0.0000001)
                                                  //perhaps this should be +ve??
               || (tri[3]>=-0.0000001 && tri[4]>=-0.0000001 && tri[5]>=-0.0000001))
            {
               inside = true;
               return_value = index_guess[i];
               break;
            }
         }
      }                                           //end of for loop over i

      if(!inside)                                 //Add the neighboring elements to the search
      {
         bool already_there;
         vector<long> index_guess_old;
         index_guess_old.resize(index_guess.size());
         //for (long j=0; j!=(long)index_guess.size(); ++j)
         //    index_guess_old[j] = index_guess[j];
         index_guess_old = index_guess;           //not quite sure this works
         //CLEAR OLD VALUES FROM INDEX_GUESS AT THIS POINT?
         for (long j=0; j!=(long)index_guess_old.size(); ++j)
         {

            elem = m_msh->ele_vector[index_guess_old[j]];
            elem->GetNeighbors(neighbors);
            num_face = elem->GetFacesNumber();

            for(long k =0; k!=num_face; ++k)
            {
               if(neighbors[k]->GetIndex() >= 0)
               {
                  already_there = false;
                  for(long l=0; l!=(long)index_guess.size(); ++l)
                  {
                     if(index_guess[l] == neighbors[k]->GetIndex())
                        already_there = true;
                  }
                  if(!already_there)
                     index_guess.push_back(neighbors[k]->GetIndex());
               }
            }

         }
      }                                           // end of if !inside

   }                                              //end of while
   //----------------------------------------------which element is point in?------------------------------------------------------END

   if(count<=number_of_elem)
      return return_value;
   else
   {
      std::cout<<"Error 1 in ELEWhatElemIsPointIn. Starting element: "<<index<<"\n";
      abort();
   }
}
#endif

/**************************************************************************
MSHLib-Method:
07/2007 OK Implementation
**************************************************************************/
bool CompleteMesh(std::string pcs_name)
{
   bool succeed = false;
   for(int i=0;i<(int)fem_msh_vector.size(); i++)
   {
      if(fem_msh_vector[i]->pcs_name.compare(pcs_name)==0)
      {
         fem_msh_vector[i]->ConstructGrid();
         fem_msh_vector[i]->FillTransformMatrix();
         succeed = true;
      }
   }
   return succeed;
}


/**************************************************************************/
/* ROCKFLOW - Function: MSHGetNextNode
 */
/* Task:
   Find the next node to the starting node in user defined direction
                                                                          */
/* Parameter: (I: Input; R: Return; X: Both)
   I: startnode, direction
                                                                          */
/* Return:
   nextnode
                                                                          */
/* Programming:
   09/2002     MB         First Version
   08/2005     MB
                                                                          */
/**************************************************************************/
long MSHGetNextNode (long startnode, CFEMesh* m_msh)
{
   size_t NumberOfNodes (m_msh->nod_vector.size());
   long NumberOfNodesPerLayer = NumberOfNodes / (m_msh->getNumberOfMeshLayers () + 1);
   return startnode + NumberOfNodesPerLayer;
}


/**************************************************************************/
/* ROCKFLOW - Function: MSHSelectFreeSurfaceNodes
 */
/* Task:
   Selection of free surface nodes, i.e. Setting free surface node flag = 1
   for the uppermost row and free surface node flag = 2 for the lowermost
   row. (Lowermost row is not moving)
                                                                          */
/* Parameter: (I: Input; R: Return; X: Both)   - void -
 */
/* Return:
   - void -
                                                                          */
/* Programming:
   03/2003     MB   First Version
   09/2004     MB   PCS
   08/2005     MB msh
                                                                          */
/**************************************************************************/
void MSHSelectFreeSurfaceNodes (CFEMesh* m_msh)
{
   // Number of nodes per node layer
	size_t NumberOfNodesPerLayer = m_msh->nod_vector.size()
			/ (m_msh->getNumberOfMeshLayers() + 1);
	size_t no_unconfined_layer = 0;
	// create array with nodes in vertical column
	size_t *strang(new size_t[m_msh->getNumberOfMeshLayers()]);
	for (size_t i = 0; i < NumberOfNodesPerLayer; i++) {
		if (m_msh->nod_vector[i]->free_surface == 4) {
			size_t nextnode = i;
			no_unconfined_layer = 0;
			for (size_t j = 0; j < m_msh->getNumberOfMeshLayers(); j++) {
				//				strang = (long*) Realloc(strang,(j+1)*sizeof(long));
				strang[j] = nextnode;
				size_t startnode = nextnode;
				nextnode = MSHGetNextNode(startnode, m_msh);
				if (m_msh->nod_vector[nextnode]->free_surface == 4) {
					strang[j + 1] = nextnode;
					no_unconfined_layer++;
				} else {
					continue;
				}
			}
		} //endif free_surface==4

		// mark start of vertical column with 1 - end of column with 2
		// this is than used in MSHMoveNODUcFlow
		m_msh->nod_vector[strang[0]]->free_surface = 1;
		m_msh->nod_vector[strang[no_unconfined_layer]]->free_surface = 2;
	} /*endfor*/
	delete[] strang;
}


/**************************************************************************
FEMLib-Method:
Task: Searches mobile nodes and sets node->free_surface = 4
Programing:
09/2004 OK / MB Implementation
05/2005 OK Bugfix
07/2005 MB MMP keyword
08/2005 MB m_pcs
**************************************************************************/
void MSHDefineMobile(CRFProcess*m_pcs)
{
   long* mobile_nodes = NULL;
   long no_mobile_nodes = -1;
   long i;
   CMediumProperties *m_mat_mp = NULL;

   //----------------------------------------------------------------------
   // Define mobile MSH nodes
   //----------------------------------------------------------------------
   // MMP Groups
   if(mmp_vector.size()==0) return;
   ////Schleife �ber alle Gruppen
   for(i=0;i<(long) mmp_vector.size();i++)
   {
      m_mat_mp = mmp_vector[i];

      //WW    int test = m_pcs->m_msh->GetMaxElementDim();
      //m_pcs->m_msh->cross_section

      //if (m_mat_mp->unconfined_flow_group ==1 && m_pcs->m_msh->GetMaxElementDim() == 3){
      if ((m_mat_mp->unconfined_flow_group ==1 && m_pcs->m_msh->GetMaxElementDim() == 3) || m_pcs->m_msh->hasCrossSection())
      {

         //if (m_mat_mp->unconfined_flow_group ==1){
         //if (m_pcs->m_msh->cross_section){
         //....................................................................
         //DOMAIN
         if(m_mat_mp->geo_type_name.find("DOMAIN")!=std::string::npos)
         {
            //CGLDomain *m_domain = NULL;
            //m_domain = m_domain->Get(m_mat_mp->geo_name);
            //mobile_nodes = m_domain->GetPointsIn(&no_mobile_nodes);
            //ToDo einlesen von domains ????
            for(i=0; i < (long)m_pcs->m_msh->nod_vector.size(); i++)
            {
               mobile_nodes = (long *) Realloc(mobile_nodes,sizeof(long)*(i+1));
               mobile_nodes[i] = i;
            }
            no_mobile_nodes = (long)m_pcs->m_msh->nod_vector.size();
         }
         //....................................................................
         //LAYER
         if(m_mat_mp->geo_type_name.find("LAYER")!=std::string::npos)
         {
            //TODO next version, change to msh file !!!
         }
         //....................................................................
         //SURFACE
         if(m_mat_mp->geo_type_name.find("SURFACE")!=std::string::npos)
         {
            Surface *m_surface = NULL;
                                                  //CC
            m_surface = GEOGetSFCByName(m_mat_mp->geo_name);
                                                  //CC
            mobile_nodes = GetPointsIn(m_surface,&no_mobile_nodes);
         }
         //....................................................................
         //VOLUME
         if(m_mat_mp->geo_type_name.find("VOLUME")!=std::string::npos)
         {
            CGLVolume *m_volume = NULL;
                                                  //CC 10/05
            m_volume = GEOGetVOL(m_mat_mp->geo_name);
            //ToDo TK
            //OK411 mobile_nodes =GetPointsInVolume(m_volume,&no_mobile_nodes);//CC 10/05
         }

      }                                           //end if unconfined flow group

   }                                              //end for mmp vector

   //----------------------------------------------------------------------
   // Set mobile MSH nodes flag
   for(i=0;i<no_mobile_nodes;i++)
   {
      m_pcs->m_msh->nod_vector[i]->free_surface = 4;
   }
   //----------------------------------------------------------------------
   if (no_mobile_nodes > 0)
   {
      m_pcs->mobile_nodes_flag = 1;
      MSHSelectFreeSurfaceNodes(m_pcs->m_msh);
   }
}


/**************************************************************************
   ROCKFLOW - Function: MSHGetNodesInColumn

   Task:
   Gets nodes of a column searching downward from startnode.

   Parameter: (I: Input; R: Return; X: Both)
           I: long node, int anz_zeilen

   Return:
           *long strang

Programming:
09/2002   MB   First Version
08/2005   MB   m_msh
**************************************************************************/
long* MSHGetNodesInColumn(long nextnode, int anz_zeilen, CFEMesh* m_msh)
{
   int i;
   long startnode;
   long *strang=NULL;

   for (i = 0; i< anz_zeilen +1; i++)
   {
      strang = (long*) Realloc(strang,(i+1)*sizeof(long));
      strang[i] = nextnode;
      startnode = nextnode;
      //nextnode = MSHGetNextNode (startnode, direction);
      nextnode = MSHGetNextNode (startnode, m_msh);
   }
   return strang;
}


/**************************************************************************/
/* ROCKFLOW - Function: MSHMoveNODUcFlow
 */
/* Task:
   Moves free surface nodes according to the pressure distribution
                                                                          */
/* Parameter: (I: Input; R: Return; X: Both)   - void -
 */
/* Return:
   - void -
                                                                          */
/* Programming:
   09/2002     MB       First Version
   05/2003     MB       verallgemeinert f�r Prismen und Vierecke
   09/2004     MB       Methode vereinfacht
   09/2004     MB       PCS
  08/2005      MB       m_msh                                                                   */
/**************************************************************************/
void MSHMoveNODUcFlow (CRFProcess*m_pcs)
{
   long nextnode = -1;
   long startnode;
   int index, index2;
   long node;
   int anz_zeilen = 0;
   int i;
   double spanne_ges;
   double spanne_rel;
   long* strang=NULL;
   double head =0.0;
   int xxflag;
   int nidy;
   // Number of nodes per node layer
   const long NumberOfNodesPerLayer (m_pcs->m_msh->nod_vector.size() / (m_pcs->m_msh->getNumberOfMeshLayers() + 1));
   double MinThickness = 1e-1;                    //OKMB
   double z_bottom;                               //OKMB

   for (node = 0; node < NumberOfNodesPerLayer; node++)
   {

      index = m_pcs->m_msh->nod_vector[node]->free_surface;
      if (index == 1)
      {
         /* Z�hlen der Zeilen (-> anz_zeilen) */
         anz_zeilen = 0;
         xxflag = 0;
         nextnode = node;
         do
         {
            startnode = nextnode;
            nextnode = MSHGetNextNode (startnode, m_pcs->m_msh);

            /* Test2: Geh�rt der n�chste Knoten zu unterer Reihe ==> Abbruch */
            index2 = m_pcs->m_msh->nod_vector[nextnode]->free_surface;

            if (index2 == 2)
            {
               xxflag = 1;
            }
            anz_zeilen++;                         /* Anzahl der beweglichen Zeilen (ohne die feste untere Zeile) */
         } while (xxflag != 1);
         /** Ende Z�hlen der Zeilen zwischen den oberen free surface node etc... und den Unteren **/

         /* Die Knoten unterhalb eines Free Surface Knotens bilden einen Strang */
         /* Die Knoten eines Stranges werden zwischengespeichert */
         strang = MSHGetNodesInColumn(node, anz_zeilen, m_pcs->m_msh);

         /* Die Knoten eines Stranges werden entsprechend der neuen Druckverteilung  verformt */
         /* Standrohrspiegelh�he bestimmen */
         nidy = m_pcs->GetNodeValueIndex("HEAD")+1;
         if (GetRFProcessDensityFlow())           /* mit Dichteunterschiede */
         {
            //OK_MOD     head = MODCalcHeadInColumn_MB(strang, anz_zeilen);
         }
         else                                     /* ohne Dichteunterschiede */
         {
            head = m_pcs->GetNodeValue(strang[0],nidy);
         }

         /* nicht �ber surface elevation */
         CRFProcess* m_pcs_OLF = NULL;
         m_pcs_OLF = PCSGet("OVERLAND_FLOW");
         double SurfaceZ;

         if(m_pcs_OLF!=NULL)
         {
            SurfaceZ = m_pcs_OLF->m_msh->nod_vector[strang[0]]->Z();
            if (head > SurfaceZ)
            {
               head = SurfaceZ;
            }
         }

         /* Set minimum thickness */
         z_bottom = m_pcs->m_msh->nod_vector[strang[anz_zeilen]]->Z();
         if(head - z_bottom < MinThickness)
            head = z_bottom + MinThickness;

         /* Berechnung der Differenz */
         spanne_ges = head - z_bottom;
         spanne_rel = spanne_ges / anz_zeilen;
         m_pcs->m_msh->nod_vector[strang[0]]->SetZ(head);

         if(spanne_ges != 0)
         {
            /* Setzen der neuen Z-Werte entlang eines Stranges */
            for (i = 1; i< anz_zeilen; i++)       /* Schleife �ber Anzahl der Zeilen */
            {
               m_pcs->m_msh->nod_vector[strang[i]]->SetZ(head -i * spanne_rel);
            }
         }

         strang = (long*) Free(strang);
      }                                           /*endif index ==1 */
   }                                              /* end for Schleife �ber alle Knoten */
}


/**************************************************************************
FEMLib-Method:
Task: Searches mobile nodes and sets node->free_surface = 4
Programing:
09/2004 OK / MB Implementation
05/2005 OK Bugfix
07/2005 MB MMP keyword
08/2005 MB m_pcs
01/2006 OK LAYER
**************************************************************************/
void CFEMesh::DefineMobileNodes(CRFProcess*m_pcs)
{
   long* mobile_nodes = NULL;
   long no_mobile_nodes = -1;
   long i,j;
   //----------------------------------------------------------------------
   // Define mobile MSH nodes
   //----------------------------------------------------------------------
   //......................................................................
   //DOMAIN
   if(m_pcs->geo_type.find("DOMAIN")!=std::string::npos)
   {
      for(i=0;i<(long)nod_vector.size();i++)
      {
         mobile_nodes = (long *) Realloc(mobile_nodes,sizeof(long)*(i+1));
         mobile_nodes[i] = i;
      }
      no_mobile_nodes = (long)m_pcs->m_msh->nod_vector.size();
   }
   //......................................................................
   //LAYER
   if(m_pcs->geo_type.find("LAYER")!=std::string::npos)
   {
      std::string m_string;
      long no_nodes_per_layer = (long)nod_vector.size() / (getNumberOfMeshLayers()+1);
      int pos = 0;
      int layer_start=0,layer_end=0;
      if(m_pcs->geo_type_name.find("-")!=std::string::npos)
      {
         pos = m_pcs->geo_type_name.find("-")!=std::string::npos;
         m_string = m_pcs->geo_type_name.substr(0,pos);
         layer_start = strtol(m_string.c_str(),NULL,0);
         m_string = m_pcs->geo_type_name.substr(pos+1,std::string::npos);
         layer_end = strtol(m_string.c_str(),NULL,0);
      }
      else
      {
         layer_start = strtol(m_pcs->geo_type_name.c_str(),NULL,0);
         layer_end = layer_start;
      }
      int no_layers = layer_end-layer_start+1;
      no_mobile_nodes = (no_layers+1)*no_nodes_per_layer;
      mobile_nodes = new long[no_mobile_nodes];
      for(i=0;i<no_layers+1;i++)
      {
         for(j=0;j<no_nodes_per_layer;j++)
         {
            mobile_nodes[i*no_nodes_per_layer+j] = j + (layer_start-1+i)*no_nodes_per_layer;
         }
      }
   }
   //......................................................................
   //SURFACE
   if(m_pcs->geo_type.find("SURFACE")!=std::string::npos)
   {
      Surface *m_sfc = NULL;
                                                  //CC
      m_sfc = GEOGetSFCByName(m_pcs->geo_type_name);
      if(m_sfc)
                                                  //CC
         mobile_nodes = GetPointsIn(m_sfc,&no_mobile_nodes);
      else
         std::cout << "Warning in CFEMesh::DefineMobileNodes - no GEO data" << std::endl;
   }
   //......................................................................
   //VOLUME
   /*OK411
     if(m_pcs->geo_type.find("VOLUME")!=std::string::npos)
     {
       CGLVolume *m_vol = NULL;
       m_vol = GEOGetVOL(m_pcs->geo_type_name);//CC 10/05
       if(m_vol)
         mobile_nodes = GetPointsInVolume(m_vol,&no_mobile_nodes);//CC 10/05
       else
         std::cout << "Warning in CFEMesh::DefineMobileNodes - no GEO data" << endl;
     }
   */
   //----------------------------------------------------------------------
   // Set mobile MSH nodes flag
   //----------------------------------------------------------------------
   for(i=0;i<(long)nod_vector.size();i++)
   {
      nod_vector[i]->free_surface = -1;
   }
   for(i=0;i<no_mobile_nodes;i++)
   {
      nod_vector[i]->free_surface = 4;
      //nod_vector[mobile_nodes[i]]->free_surface = 4;
   }
   //----------------------------------------------------------------------
   if (no_mobile_nodes > 0)
   {
      m_pcs->mobile_nodes_flag = 1;
      MSHSelectFreeSurfaceNodes(this);
   }
   //----------------------------------------------------------------------
   delete [] mobile_nodes;
   mobile_nodes = NULL;
}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
09/2004 OK Implementation
ToDo evtl. vector<CGLPoint>
08/2005 CC Modification: CGLPoint* e_pnt - Move from GeoLib to MshLib
**************************************************************************/
void MSHGetNodesClose(std::vector<long>&msh_point_vector,CGLPoint* e_pnt)
{
   e_pnt = e_pnt;
   msh_point_vector.size();
   /*OK411
     long i;
     CGLPoint m_pnt;
     // Node loop
     for (i=0;i<NodeListSize();i++) {
       if (GetNode(i)==NULL) continue;
       m_pnt.x = GetNodeX(i);
       m_pnt.y = GetNodeY(i);
       m_pnt.z = GetNodeZ(i);
       if(e_pnt->PointDis(&m_pnt)<=(e_pnt->epsilon+MKleinsteZahl))
         msh_point_vector.push_back(i);
   }
   */
}


/*************************************************************************
  ROCKFLOW - Function: MSHGetNodesClose
  Task: Searching grid points which are close to a polyline
  Programming:
   10/2002 OK Encapsulated from ExecuteSourceSinkMethod11 (CT)
   01/2003 OK Test
  last modified: 20.01.2003 OK
   08/2005 CC Modification Move from GeoLib to MSHLib
**************************************************************************/
long* MSHGetNodesClose(long *number_of_nodes, CGLPolyline* m_ply)
{
   long *nodes_all = NULL;
   m_ply = m_ply;
   number_of_nodes = number_of_nodes;
   /*OK411
     long j,k,l;
     double pt1[3],line1[3],line2[3],pt0[3];
     double mult_eps = 1.0;
     double dist1p,dist2p,*length,laenge;
     long anz_relevant = 0;
     typedef struct {
        long knoten;
        long abschnitt;
        double laenge;
     } INFO;
   INFO *relevant=NULL;
   int weiter;
   double w1,w2;
   long knoten_help;
   double laenge_help;
   double gesamte_laenge = 0.;
   long polyline_point_vector_size;

   m_ply->sbuffer.clear();
   m_ply->ibuffer.clear();

   if (m_ply) {

   length = (double*) Malloc(sizeof(double) *(long)m_ply->point_vector.size());

   pt0[0] = m_ply->point_vector[0]->x;
   pt0[1] = m_ply->point_vector[0]->y;
   pt0[2] = m_ply->point_vector[0]->z;

   polyline_point_vector_size =(long)m_ply->point_vector.size();
   for (k=0;k<polyline_point_vector_size-1;k++) {
   line1[0] = m_ply->point_vector[k]->x;
   line1[1] = m_ply->point_vector[k]->y;
   line1[2] = m_ply->point_vector[k]->z;
   line2[0] = m_ply->point_vector[k+1]->x;
   line2[1] = m_ply->point_vector[k+1]->y;
   line2[2] = m_ply->point_vector[k+1]->z;
   length[k] = MCalcDistancePointToPoint(line2, line1);
   gesamte_laenge += length[k];
   }

   // Wiederholen bis zumindest ein Knoten gefunden wurde
   while(anz_relevant==0) {

   for (j=0;j<NodeListSize();j++) {
   if (GetNode(j)==NULL) continue;

   polyline_point_vector_size =(long)m_ply->point_vector.size();
   for (k=0;k<polyline_point_vector_size-1;k++) {

   pt1[0] = GetNodeX(j);
   pt1[1] = GetNodeY(j);
   pt1[2] = GetNodeZ(j);

   line1[0] = m_ply->point_vector[k]->x;
   line1[1] = m_ply->point_vector[k]->y;
   line1[2] = m_ply->point_vector[k]->z;
   line2[0] = m_ply->point_vector[k+1]->x;
   line2[1] = m_ply->point_vector[k+1]->y;
   line2[2] = m_ply->point_vector[k+1]->z;

   if ( MCalcDistancePointToLine(pt1,line1,line2) <= mult_eps*m_ply->epsilon ) {
   MCalcProjectionOfPointOnLine(pt1,line1,line2,pt1);
   dist1p = MCalcDistancePointToPoint(line1, pt1);
   dist2p = MCalcDistancePointToPoint(line2, pt1);
   if ((dist1p+dist2p-length[k]) <=  mult_eps*m_ply->epsilon ) {

   // For boundary conditions. WW
   m_ply->sbuffer.push_back(dist1p);
   m_ply->ibuffer.push_back(k);
   // ---------------------------

   anz_relevant++;
   nodes_all = (long *) Realloc(nodes_all,sizeof(long)*anz_relevant);
   relevant = (INFO *) Realloc(relevant, sizeof(INFO) * anz_relevant);
   nodes_all[anz_relevant-1] = j;
   laenge = 0.;
   for (l=0; l < k; l++)
   laenge += length[l];
   relevant[anz_relevant-1].knoten = j;
   relevant[anz_relevant-1].laenge = laenge + dist1p;
   k =(long)m_ply->point_vector.size();
   }
   }
   }
   }
   if(anz_relevant==0) mult_eps *=2.;
   }

   if (mult_eps > 1.)
   std::cout << "!!! Epsilon increased in sources!" << endl;

   do {
   weiter = 0;
   for (k=0;k<anz_relevant-1;k++) {
   w1=relevant[k].laenge;
   w2=relevant[k+1].laenge;
   if (w1>w2) { // Die Eintraege vertauschen
   knoten_help = relevant[k].knoten;
   laenge_help = relevant[k].laenge;
   relevant[k].knoten = relevant[k+1].knoten;
   relevant[k].laenge = relevant[k+1].laenge;
   relevant[k+1].knoten = knoten_help;
   relevant[k+1].laenge = laenge_help;
   weiter=1;
   }
   }
   } while (weiter);

   relevant = (INFO*) Free(relevant);
   *number_of_nodes = anz_relevant;
   }
   */
   return nodes_all;
}


/**************************************************************************
GeoLib-Method: GetPointsIn
Task:
Programing:
01/2004 OK Implementation
08/2005 CC Modification Move from Geolib to Mshlib
**************************************************************************/
long* GetPointsIn(Surface* m_sfc,long* number_of_nodes)
{
   long *nodes = NULL;
   number_of_nodes = number_of_nodes;
   m_sfc = m_sfc;
   /*OK411
     long i;
     double *xp=NULL,*yp=NULL,*zp=NULL;
     long anz_relevant = 0;
     CGLPoint m_pnt;
     // Inside polygon
     if(!m_sfc->polygon_point_vector.empty()) {
       xp = (double*) Malloc(((long)m_sfc->polygon_point_vector.size())*sizeof(double));
       yp = (double*) Malloc(((long)m_sfc->polygon_point_vector.size())*sizeof(double));
       zp = (double*) Malloc(((long)m_sfc->polygon_point_vector.size())*sizeof(double));
       long polygon_point_vector_length = (long)m_sfc->polygon_point_vector.size();
   for(i=0;i<polygon_point_vector_length;i++) {
   xp[i] = m_sfc->polygon_point_vector[i]->x;
   yp[i] = m_sfc->polygon_point_vector[i]->y;
   zp[i] = m_sfc->polygon_point_vector[i]->z;
   }

   //-----------------------------------------------------------------
   for(i=0;i<NodeListSize();i++) {
   if (GetNode(i)==NULL) continue;
   m_pnt.x = GetNodeX(i);
   m_pnt.y = GetNodeY(i);
   m_pnt.z = GetNodeZ(i);
   if(m_pnt.IsInsidePolygonPlain(
   xp,yp,zp,\
   (long)m_sfc->polygon_point_vector.size())) {
   anz_relevant++;
   nodes = (long *) Realloc(nodes,sizeof(long)*anz_relevant);
   nodes[anz_relevant-1] = i;
   }
   }
   }
   // Destructions
   // nodes extern
   xp = (double*) Free(xp);
   yp = (double*) Free(yp);
   zp = (double*) Free(zp);
   //
   *number_of_nodes = anz_relevant;
   */
   return nodes;
}


/**************************************************************************
GeoLib-Method: GEOGetVolume
Task: Get element nodes in a material domain
Programing:
10/2004 WW Implementation
**************************************************************************/
void GEOGetNodesInMaterialDomain(const int MatIndex, std::vector<long>& Nodes)
{
   (void)MatIndex;
   (void)Nodes;
   /*OK411
   MatIndex;
   Nodes.size();
    long index, *element_nodes;
    int i, j, Size, nn, order = 2;
    const int L_Nodes = GetLowOrderNodeNumber();
    bool exist;
    if(L_Nodes==NodeListSize()) order = 1;
    if(L_Nodes==0) order = 1;

    Nodes.resize(0);
   nn = 0;
   for (index=0;index<ElListSize();index++)
   {
   if (ElGetElement(index)!=NULL)
   {  // Eelement exist
   if (ElGetElementActiveState(index))
   {  // Element active
   if(order==1) nn = NumbersOfElementNode(index);
   if(order==2) nn = NumbersOfElementNodeHQ(index);
   if(ElGetElementGroupNumber(index)==MatIndex)
   {
   Size = (int)Nodes.size();
   element_nodes = ElGetElementNodes(index);
   for(i=0; i<nn; i++)
   {
   exist = false;
   for(j=0; j<Size; j++)
   {
   if(element_nodes[i]==Nodes[j])
   {
   exist = true;
   break;
   }
   }
   if(!exist) Nodes.push_back(element_nodes[i]);
   }
   }
   }
   }
   }
   */
}


/**************************************************************************
GeoLib-Method: GEOGetVolume
Task: Get element nodes in a material domain
Programing:
10/2004 WW Implementation
**************************************************************************/
//using Mesh_Group::CElem;
void GEOGetNodesInMaterialDomain(CFEMesh* m_msh, const int MatIndex, std::vector<long>& Nodes, bool Order)
{
   int i, nn;
   long e, j, Size;
   CElem* elem = NULL;
   bool exist;
   nn =0;
   Size = 0;
   Nodes.resize(0);
   for (e = 0; e < (long)m_msh->ele_vector.size(); e++)
   {
      elem = m_msh->ele_vector[e];
      if (elem->GetMark())                        // Marked for use
      {
         nn = elem->GetNodesNumber(Order);
         if(elem->GetPatchIndex()==MatIndex)
         {
            Size = (int)Nodes.size();
            for(i=0; i<nn; i++)
            {
               exist = false;
               for(j=0; j<Size; j++)
               {
                  if(elem->GetNodeIndex(i)==Nodes[j])
                  {
                     exist = true;
                     break;
                  }
               }
               if(!exist) Nodes.push_back(elem->GetNodeIndex(i));
            }
         }
      }                                           // if
   }                                              //For
}


/**************************************************************************
GeoLib-Method:
Task:
Programing:
01/2004 OK Implementation based on algorithm originally by CT
09/2005 CC Move from geolib to mshlib
**************************************************************************/
void SetRFIPointsClose(CGLLine* m_lin)
{
   m_lin = m_lin;
   /*OK411
     long j;
     double pt1[3],pt2[3],line1[3],line2[3];
     double mult_eps = 1.0;
     double dist1p,dist2p,length;
     long anz_relevant;
     double *dist;
     typedef struct {
        long knoten;
        long abschnitt;
        double laenge;
   } INFO;
   INFO *relevant=NULL;
   long knoten_help;
   double laenge_help;
   double w1,w2;
   int weiter;
   //----------------------------------------------------------------------
   // Tests
   if(!ELEListExists()) {
   return;
   }
   if(m_lin->point1<0)
   return;
   if(m_lin->point2<0)
   return;
   //----------------------------------------------------------------------
   // Initializations
   anz_relevant = 0;
   m_lin->msh_nodes = NULL;
   CGLPoint *m_point=NULL;
   m_point = GEOGetPointById(m_lin->point1);//CC
   line1[0] = m_point->x;
   line1[1] = m_point->y;
   line1[2] = m_point->z;
   m_point = GEOGetPointById(m_lin->point2);//CC
   line2[0] = m_point->x;
   line2[1] = m_point->y;
   line2[2] = m_point->z;
   length = MCalcDistancePointToPoint(line2,line1);
   //----------------------------------------------------------------------
   // Repeat untill at least one node is found
   while(anz_relevant==0) {
   // NOD list loop
   for (j=0;j<NodeListSize();j++) {
   if (GetNode(j)==NULL) continue;
   pt1[0] = GetNodeX(j);
   pt1[1] = GetNodeY(j);
   pt1[2] = GetNodeZ(j);
   // Is MSH point near to line
   if ( MCalcDistancePointToLine(pt1,line1,line2) <= mult_eps*m_lin->epsilon ) {
   // Calc projection of pt1 to line and use this in the following
   MCalcProjectionOfPointOnLine(pt1,line1,line2,pt1);
   // Abstand des Punktes zum ersten Punkt des Polygonabschnitts
   dist1p = MCalcDistancePointToPoint(line1, pt1);
   // Abstand des Punktes zum zweiten Punkt des Polygonabschnitts
   dist2p = MCalcDistancePointToPoint(line2, pt1);
   // Ist der Knoten innerhalb des Intervalls?
   if ((dist1p+dist2p-length) <=  mult_eps*m_lin->epsilon ) {
   anz_relevant++;
   // Feld anpassen
   m_lin->msh_nodes = (long *) Realloc(m_lin->msh_nodes,sizeof(long)*anz_relevant);
   relevant = (INFO *) Realloc(relevant, sizeof(INFO)*anz_relevant);
   // Ablegen von Knotennummer und Position
   m_lin->msh_nodes[anz_relevant-1] = j;
   }
   } // endif
   } // Ende Schleife ueber Knoten
   if(anz_relevant==0) mult_eps *=2.;
   } // Ende Schleife Wiederholungen
   if (mult_eps > 1.)
   printf("Warning: Epsilon increased in CGLLine::SetRFIPointsClose.");
   m_lin->no_msh_nodes = anz_relevant;
   //----------------------------------------------------------------------
   // Sort MSH nodes, beginning from first line point
   //......................................................................
   // Calc distances from first line point
   m_point = GEOGetPointById(m_lin->point1);//CC
   pt1[0] = m_point->x;
   pt1[1] = m_point->y;
   pt1[2] = m_point->z;
   dist = (double*) Malloc(sizeof(double)*m_lin->no_msh_nodes);
   for(j=0;j<m_lin->no_msh_nodes;j++) {
   pt2[0] = GetNodeX(m_lin->msh_nodes[j]);
   pt2[1] = GetNodeY(m_lin->msh_nodes[j]);
   pt2[2] = GetNodeZ(m_lin->msh_nodes[j]);
   dist[j] = MCalcDistancePointToPoint(pt1,pt2);
   }
   //......................................................................
   // Sorting by pair switching
   do {
   weiter = 0;
   for (j=0;j<m_lin->no_msh_nodes-1;j++) {
   w1=dist[j];
   w2=dist[j+1];
   if (w1>w2) { // Die Eintraege vertauschen
   knoten_help = m_lin->msh_nodes[j];
   laenge_help = dist[j];
   m_lin->msh_nodes[j] = m_lin->msh_nodes[j+1];
   dist[j] = dist[j+1];
   m_lin->msh_nodes[j+1] = knoten_help;
   dist[j+1] = laenge_help;
   weiter=1;
   }
   }
   } while (weiter);
   //----------------------------------------------------------------------
   // Destructions
   dist = (double*) Free(dist);
   relevant = (INFO*) Free(relevant);
   */
}


/**************************************************************************
GeoSys-GUI Function
Programing:
01/2004 OK Implementation
09/2005 CC Modification No MFC function
**************************************************************************/
//bool IsPointInSurface(Surface* m_fsc, CGLPoint *m_point)
//{
//  bool ok = false;
//  m_point = m_point;
//  if(!m_fsc)
//    return false; //OK
//  return ok;
//}
