/**************************************************************************
FEMLib - Object: Initial Conditions IC
Task:
Programing:
08/2004 OK Implementation
12/2005 OK Restart
last modified
**************************************************************************/

#include "makros.h"
#include "gs_project.h"
// C++ STL
#include <list>
#include <iostream>
using namespace std;
// FEM-Makros
//#include "files0.h"
#include "mathlib.h"
#include "files0.h"
// Base
//#include "StringTools.h"
// GEOLIB
#include "GEOObjects.h"
// MSHLib
#include "msh_lib.h"
// FEMLib
#include "rf_ic_new.h"
#include "rf_pcs.h"
#include "rf_node.h"
#include "rfmat_cp.h"

//==========================================================================
vector<CInitialConditionGroup*>ic_group_vector;
vector<CInitialCondition*>ic_vector;

/**************************************************************************
FEMLib-Method:
Task:
Programing:
08/2004 OK Implementation
**************************************************************************/
CInitialCondition::CInitialCondition()
{
   //  geo_type_name = "DOMAIN";
   this->setProcessDistributionType(FiniteElement::CONSTANT);
   // HS: not needed, removed.
   // m_node = new CNodeValue();
   // m_node->node_value = 0.0;
   SubNumber=0;
   this->setProcess(NULL);                        //OK
   a0=b0=c0=d0=NULL;                              //WW
   m_msh = NULL;                                  //OK
   m_node = NULL;                                 //HS
}


/**************************************************************************
FEMLib-Method:
Task: BC deconstructor
Programing:
04/2004 OK Implementation
**************************************************************************/
CInitialCondition::~CInitialCondition(void)
{
   if ( !node_value_vector.empty() )
   {
      for (int i=0; i < (int)node_value_vector.size(); i++)
         delete node_value_vector[i];
   }
   node_value_vector.clear();
   delete m_node;                                 //HS
   //WW
   if(a0) delete [] a0;
   if(b0) delete [] b0;
   if(c0) delete [] c0;
   if(d0) delete [] d0;
   a0=b0=c0=d0=NULL;
}


/**************************************************************************
FEMLib-Method:
Task: IC read function
Programing:
08/2004 OK Implementation
01/2005 OK Boolean type
01/2005 OK Destruct before read
05/2010 TF reformated, restructured, signature changed, use new GEOLIB data structures
**************************************************************************/
bool ICRead(const std::string& file_base_name,
const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
   // File handling
   std::string ic_file_name = file_base_name + IC_FILE_EXTENSION;
   std::ifstream ic_file(ic_file_name.data(), std::ios::in);
   if (!ic_file.good())
   {
      std::cout << "WARNING: ICRead: No initial conditions !" << std::endl;
      return false;
   }

   char line[MAX_ZEILE];
   std::string line_string;
   std::ios::pos_type position;

   // Keyword loop
   std::cout << "ICRead" << std::endl;
   while (!ic_file.eof())
   {
      ic_file.getline(line, MAX_ZEILE);
      line_string = line;
      if (line_string.find("#STOP") != string::npos)
         return true;

                                                  // keyword found
      if (line_string.find("#INITIAL_CONDITION") != std::string::npos)
      {
         CInitialCondition *ic = new CInitialCondition();
         std::ios::pos_type pos (ic_file.tellg());
         position = ic->Read(&ic_file, geo_obj, unique_name);
         if (pos != position)
         {
            ic_vector.push_back(ic);
         }
         else
         {
            std::cerr << "WARNING: in ICRead: could not read initial condition" << std::endl;
            delete ic;
         }
         ic = NULL;
         ic_file.seekg(position, ios::beg);
      }                                           // keyword found
   }                                              // eof
   return true;
}


/**************************************************************************
FEMLib-Method: ICWrite
Task: master write function
Programing:
04/2004 OK Implementation
last modification:
**************************************************************************/
void ICWrite(string base_file_name)
{
   CInitialCondition *m_ic = NULL;
   string sub_line;
   string line_string;
   long i;
   //========================================================================
   // File handling
   string ic_file_name = base_file_name + IC_FILE_EXTENSION;
   fstream ic_file (ic_file_name.data(),ios::trunc|ios::out);
   ic_file.setf(ios::scientific,ios::floatfield);
   ic_file.precision(12);
   if (!ic_file.good()) return;
   ic_file.seekg(0L,ios::beg);
   //========================================================================
   ic_file << "GeoSys-IC: Initial Conditions ------------------------------------------------\n";
   //========================================================================
   // IC list
   long no_ic =(long)ic_vector.size();
   for(i=0;i<no_ic;i++)
   {
      m_ic = ic_vector[i];
      m_ic->Write(&ic_file);
   }
   ic_file << "#STOP";
   ic_file.close();
}


/**************************************************************************
FEMLib-Method:
Task: ST read function
Programing:
08/2004 OK Implementation (DOMAIN)
10/2004 OK POLYLINE
11/2004 CMCD String streaming, Gradient
01/2005 WW Initialize by sub-domain
05/2005 OK PRIMARY_VARIABLE
12/2005 OK RESTART
07/2006 WW Read data by expression
06/2010 TF changed signature to use the new GEOLIB data structures
**************************************************************************/
ios::pos_type CInitialCondition::Read(std::ifstream *ic_file,
const GEOLIB::GEOObjects& geo_obj, const std::string& unique_geo_name)
{
   string line_string;
   std::stringstream in;
   ios::pos_type position;

   bool new_keyword = false;
   int ibuf (0);
   double d_buf (0.0);

   // read loop
   while (!new_keyword)
   {
      position = ic_file->tellg();
      line_string = GetLineFromFile1(ic_file);
      if (line_string.size() < 1)
         break;
      if (line_string.find("#") != string::npos)
      {
         new_keyword = true;
         break;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$PCS_TYPE") != string::npos)
      {
         in.str(GetLineFromFile1(ic_file));
         std::string tmp;
         in >> tmp;                               // pcs_type_name;
         this->setProcessType(convertProcessType(tmp));
         in.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$PRIMARY_VARIABLE") != string::npos)
      {
         in.str(GetLineFromFile1(ic_file));
         std::string tmp;
         in >> tmp;                               // pcs_pv_name;
         if ( this->getProcessType() == MASS_TRANSPORT )
         {
             // HS set the pointer to MCP based on component name.
             // first do a check whether this name is existing and unique. 
             if ( cp_name_2_idx.count( tmp ) == 1 )
             {
                 setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess() );
                 setProcessPrimaryVariable( CONCENTRATION );
             }
             else
             {
                 DisplayErrorMsg("Error: In reading IC file, the input component names are not found in MCP file!!!");
                 exit(1);             
             }
         }
         else
         {
             setProcess( PCSGet( this->getProcessType() ) );
             setProcessPrimaryVariable (convertPrimaryVariable (tmp));
         }
         in.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$COMP_NAME") != string::npos)
      {
         in.str(GetLineFromFile1(ic_file));
         std::string tmp;
         in >> tmp;
         // HS set the pointer to MCP based on component name.
         // first do a check whether this name is existing and unique. 
         if ( cp_name_2_idx.count( tmp ) == 1 )
         {
             setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess() );
             setProcessPrimaryVariable( CONCENTRATION );
         }
         else
         {
             DisplayErrorMsg("Error: In reading BC file, the input component names are not found in MCP file!!!");
             exit(1);             
         }
         in.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$DIS_TYPE") != string::npos)
      {
         in.str(GetLineFromFile1(ic_file));
         std::string tmp;
         in >> tmp;                               // dis_type_name;
         this->setProcessDistributionType(FiniteElement::convertDisType(tmp));

         //if (dis_type_name.find("CONSTANT") != string::npos) {
         if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
         {
            m_node = new CNodeValue();
            in >> m_node->node_value;
            node_value_vector.push_back(m_node);
            m_node = NULL;
         }
         //if (dis_type_name.find("GRADIENT") != string::npos) {
         if (this->getProcessDistributionType() == FiniteElement::GRADIENT)
         {
            in >> gradient_ref_depth;             //CMCD
            in >> gradient_ref_depth_value;       //CMCD
            in >> gradient_ref_depth_gradient;    //CMCD
         }
         //if (dis_type_name.find("RESTART") != string::npos) { //OK
         if (this->getProcessDistributionType() == FiniteElement::RESTART)
         {
            in >> rfr_file_name;
         }
         in.clear();
         continue;
      }
      //....................................................................
                                                  //subkeyword found
      if (line_string.find("$GEO_TYPE") != string::npos)
      {
         in.str(GetLineFromFile1(ic_file));
         std::string geo_type_name;
         in >> geo_type_name;
         if (geo_type_name.find("POINT") != string::npos)
         {
            in >> geo_name;
            const GEOLIB::Point *pnt ((geo_obj.getPointVecObj(unique_geo_name))->getElementByName (geo_name));
            if (pnt == NULL)
            {
               std::cerr << "error in CInitialCondition::Read: point name \"" << geo_name << "\" not found!" << std::endl;
               exit (1);
            }
            setGeoType (GEOLIB::POINT);
            setGeoObj (pnt);

            in.clear();
            geo_name = "";                        // REMOVE CANDIDATE
         }
         if (geo_type_name.find("POLYLINE") != string::npos)
         {
            in >> geo_name;
            const GEOLIB::Polyline *ply ((geo_obj.getPolylineVecObj(unique_geo_name))->getElementByName (geo_name));
            if (ply == NULL)
            {
               std::cerr << "error in CInitialCondition::Read: polyline name \"" << geo_name << "\" not found!" << std::endl;
               exit (1);
            }
            setGeoType (GEOLIB::POLYLINE);
            setGeoObj (ply);

         }
         if (geo_type_name.find("SURFACE") != string::npos)
         {
            setGeoType (GEOLIB::SURFACE);
            //				// TF 07/2010 - get the surface vector and get the surface ID
            //				if (!((geo_obj.getSurfaceVecObj(unique_name))->getElementIDByName(
            //						geo_name, _geo_obj_idx))) {
            //					std::cerr
            //							<< "ERROR: CInitialCondition::Read: surface name not found!"
            //							<< std::endl;
            //					exit(1);
            //				}
            in >> geo_name;
         }
         if (geo_type_name.find("VOLUME") != string::npos)
         {
            setGeoType (GEOLIB::VOLUME);
         }
         if (geo_type_name.find("DOMAIN") != string::npos)
         {
            setGeoType (GEOLIB::GEODOMAIN);
            //  Give initial condition by patches of domain. WW
            if (geo_type_name.find("SUB") != string::npos)
            {
               *ic_file >> SubNumber;
               if (convertPrimaryVariableToString(this->getProcessPrimaryVariable()).find("STRESS") != string::npos)
                  //     ||dis_type_name.find("FUNCTION")!=string::npos) //01.07.2008 WW
               {
                  string str_buff;
                  vector<string> tokens;
                  stringstream buff;
                  char *pch;
                  char seps[] = "+\n";
                  char seps1[] = "*";
                  double f_buff;
                  a0 = new double[SubNumber];
                  b0 = new double[SubNumber];
                  c0 = new double[SubNumber];
                  d0 = new double[SubNumber];
                  for (size_t i = 0; i < SubNumber; i++)
                  {
                     a0[i] = b0[i] = c0[i] = d0[i] = 0.0;
                     *ic_file >> ibuf >> str_buff >> ws;
                     subdom_index.push_back(ibuf);
                     pch = strtok(const_cast<char*> (str_buff.c_str()),
                        seps);
                     buff << pch;
                     buff >> a0[i];
                     buff.clear();
                     while (pch != NULL)
                     {
                        pch = strtok(NULL, seps);
                        if (pch == NULL)
                           break;
                        string token = pch;
                        tokens.push_back(token);
                     }
                     for (size_t k=0; k < tokens.size(); k++)
                     {
                        pch = strtok(
                           const_cast<char*> (tokens[k].c_str()),
                           seps1);
                        buff << pch;
                        buff >> f_buff;
                        buff.clear();
                        pch = strtok(NULL, seps1);
                        switch (pch[0])
                        {
                           case 'x':
                              b0[i] = f_buff;
                              break;
                           case 'y':
                              c0[i] = f_buff;
                              break;
                           case 'z':
                              d0[i] = f_buff;
                              break;
                        }
                     }
                     tokens.clear();
                  }
               }
               else
               {
                  for (size_t i = 0; i < SubNumber; i++)
                  {
                     *ic_file >> ibuf >> d_buf;
                     subdom_index.push_back(ibuf);
                     subdom_ic.push_back(d_buf);
                  }
               }
            }
         }
         in.clear();
         continue;
      }
   }                                              // Schleife ueber alle Phasen bzw. Komponenten
   return position;
}


/**************************************************************************
FEMLib-Method:
Task: write function
Programing:
02/2004 OK Implementation
01/2005 OK Extensions
05/2005 OK PRIMARY_VARIABLE
12/2008 NW GRADIENT keyword
last modification:
**************************************************************************/
void CInitialCondition::Write(fstream* ic_file) const
{
   //KEYWORD
   *ic_file  << "#INITIAL_CONDITION" << endl;
   //--------------------------------------------------------------------
   //NAME+NUMBER
   *ic_file << " $PCS_TYPE" << endl;
   *ic_file << "  ";
   *ic_file << convertProcessTypeToString(this->getProcessType()) << endl;
   *ic_file << " $PRIMARY_VARIABLE" << endl;
   *ic_file << "  ";
   *ic_file << convertPrimaryVariableToString(this->getProcessPrimaryVariable()) << endl;
   //--------------------------------------------------------------------
   //GEO_TYPE
   *ic_file << " $GEO_TYPE" << endl;
   *ic_file << "  " << getGeoTypeAsString() << " CInitialCondition::Write ToDo write name of GeoObject " << std::endl;

   //--------------------------------------------------------------------
   //DIS_TYPE
   *ic_file << " $DIS_TYPE" << endl;
   *ic_file << "  ";
   *ic_file << convertDisTypeToString(this->getProcessDistributionType()) << endl;
   *ic_file << " ";
   if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
   {
      CNodeValue* m_node = NULL;
      int node_value_vector_size = (int)node_value_vector.size();
      if(node_value_vector_size>0)
      {
         m_node = node_value_vector[0];
         *ic_file << m_node->node_value;
      }
   }
   else if (this->getProcessDistributionType() == FiniteElement::GRADIENT)
   {
      *ic_file << this->gradient_ref_depth << " ";
      *ic_file << this->gradient_ref_depth_value << " ";
      *ic_file << this->gradient_ref_depth_gradient;
   }
   //--------------------------------------------------------------------
   *ic_file << endl;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
08/2004 OK Implementation
**************************************************************************/
void CInitialCondition::Set(int nidx)
{
   switch(getGeoType())
   {
      case GEOLIB::POINT:
         SetPoint(nidx);
         std::cout << "WARNING: CInitialCondition::Set - ToDo" << endl;
         break;
      case GEOLIB::POLYLINE:
         SetPolyline(nidx);
         break;
      case GEOLIB::SURFACE:
         SetSurface(nidx);
         break;
      case GEOLIB::VOLUME:
         std::cout << "WARNING: CInitialCondition::Set - ToDo" << endl;
         break;
      case GEOLIB::GEODOMAIN:
         SetDomain(nidx);
         break;
      case GEOLIB::INVALID:
         std::cout << "WARNING: CInitialCondition::Set - invalid geo type" << endl;
         break;
   }

   /* KR not used in benchmarks
      // Direct assign by node indeces 17.11.2009 PCH
      if(dis_type_name.find("DIRECT")!=string::npos)
     {
         SetByNodeIndex(nidx);
      }
   */
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
11/2009 PCH Implementation
09/2010 KR cleaned up code
**************************************************************************/
void CInitialCondition::SetByNodeIndex(int nidx)
{
   string line_string;
   std::stringstream in;
   long node_index;
   double node_val;

   // File handling
   ifstream d_file (fname.c_str(),ios::in);
   if (!d_file.is_open())
   {
      cout << "! Error in direct node source terms: Could not find file " << fname << endl;
      abort();
   }

   while (!d_file.eof())
   {
      line_string = GetLineFromFile1(&d_file);
      if(line_string.find("#STOP")!=string::npos)
         break;

      in.str(line_string);
      in >> node_index >> node_val;
      in.clear();

      this->getProcess()->SetNodeValue(node_index,nidx,node_val);
   }
}


/**************************************************************************
FEMLib-Method:
Task: set polyline values
Programing:
06/2006 MX Implementation
last modification:
**************************************************************************/
void CInitialCondition::SetPoint(int nidx)
{
   if (m_msh && this->getProcessDistributionType() == FiniteElement::CONSTANT)
   {
      this->getProcess()->SetNodeValue(this->getProcess()->m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(getGeoObj())), nidx,
         node_value_vector[0]->node_value);
   }
   else
   {
      std::cerr << "Error in CInitialCondition::SetPoint - point: "
         << *(static_cast<const GEOLIB::Point*>(getGeoObj()))
         << " not found" << endl;
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
08/2004 OK Implementation
**************************************************************************/
void CInitialCondition::SetEle(int nidx)
{
   switch (getGeoType())
   {
      case GEOLIB::POINT:
         std::cout << "Warning CInitialCondition::Set - ToDo" << std::endl;
         break;
      case GEOLIB::POLYLINE:
         // SetPolyline(nidx);
         std::cout << "Warning CInitialCondition::Set - ToDo" << std::endl;
         break;
      case GEOLIB::SURFACE:
         std::cout << "Warning CInitialCondition::Set - ToDo" << std::endl;
         break;
      case GEOLIB::VOLUME:
         std::cout << "Warning CInitialCondition::Set - ToDo" << std::endl;
         break;
      case GEOLIB::GEODOMAIN:
         SetDomainEle(nidx);
         break;
      case GEOLIB::INVALID:
         std::cout << "WARNING: CInitialCondition::SetEle - invalid geo type" << endl;
         break;
   }
}


/**************************************************************************
FEMLib-Method:
Task: set polyline values
Programing:
10/2004 OK Implementation
07/2006 OK/MX MSH
**************************************************************************/
void CInitialCondition::SetPolyline(int nidx)
{
   if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
   {
      //		CGLPolyline* m_polyline = polyline_vector[getGeoObjIdx()];
      //		if (m_polyline) {
      //			std::vector<long> nodes_vector;
      //			m_msh->GetNODOnPLY(m_polyline, nodes_vector);
      //			for (size_t i = 0; i < nodes_vector.size(); i++) {
      //				m_pcs->SetNodeValue(nodes_vector[i], nidx, node_value_vector[0]->node_value);
      //			}
      //		} else {
      //			std::cout << "Error in CInitialCondition::SetPolyline - polyline: "
      //					<< geo_name << " not found" << std::endl;
      //		}
      if (getGeoObj())
      {
         std::vector<long> nodes_vector;
         m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(getGeoObj()), nodes_vector);
         for (size_t i = 0; i < nodes_vector.size(); i++)
         {
            this->getProcess()->SetNodeValue(nodes_vector[i], nidx, node_value_vector[0]->node_value);
         }
      }
      else
      {
         std::cout << "Error in CInitialCondition::SetPolyline - polyline: "
            << geo_name << " not found" << std::endl;
      }
   }
   //	if (dis_type_name.find("LINEAR") != string::npos) {
   //	}
}


/**************************************************************************
FEMLib-Method:
Task: set surface values
Programing:
02/2008 JOD Implementation
**************************************************************************/
void CInitialCondition::SetSurface(int nidx)
{
   double value, node_depth;
   vector<long>sfc_nod_vector;
   Surface* m_sfc = GEOGetSFCByName(geo_name);

   if(m_sfc && m_msh)
   {

      m_msh->GetNODOnSFC(m_sfc, sfc_nod_vector);

      if(this->getProcessDistributionType() == FiniteElement::CONSTANT)
      {
         for(size_t i = 0; i < sfc_nod_vector.size(); i++)
         {
            value = node_value_vector[0]->node_value;
            this->getProcess()->SetNodeValue(sfc_nod_vector[i],nidx, value);
         }                                        // end surface nodes

      }                                           // end constant
      else if(this->getProcessDistributionType() == FiniteElement::GRADIENT)
      {
         int onZ = m_msh->GetCoordinateFlag()%10;
         long msh_node;

         for (size_t i = 0; i < sfc_nod_vector.size(); i++)
         {
            msh_node = sfc_nod_vector[i];
            if(onZ == 1)                          //2D
               node_depth = m_msh->nod_vector[msh_node]->Y();
            else if(onZ == 2)                     //3D
               node_depth = m_msh->nod_vector[msh_node]->Z();
            else
            {
               cout << "Error in CInitialCondition::SetSurface - dis_type: " << convertDisTypeToString(this->getProcessDistributionType()) << "don't know If 2D or 3D"<<  endl;
               node_depth = 0;
            }
            value = ((gradient_ref_depth_gradient) * (gradient_ref_depth-node_depth)) + gradient_ref_depth_value;
            this->getProcess()->SetNodeValue(msh_node,nidx, value);
         }                                        // end surface nodes

      }                                           // end gradient
      else
         cout << "Error in CInitialCondition::SetSurface - dis_type: " << convertDisTypeToString(this->getProcessDistributionType()) << " not found" << endl;

   }                                              // end m_sfc
   else
      cout << "Error in CInitialCondition::SetSurface - surface: " << geo_name << " not found" << endl;

}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
04/2004 OK Implementation
09/2010 KR cleaned up code
**************************************************************************/
CInitialConditionGroup* GetInitialConditionGroup(string this_pcs_name)
{
   CInitialConditionGroup *m_ic_group = NULL;
   size_t no_ic_groups = ic_group_vector.size();
   for(size_t i=0;i<no_ic_groups;i++)
   {
      m_ic_group = ic_group_vector[i];
      if (m_ic_group->pcs_type_name.find(this_pcs_name) != string::npos)
         return m_ic_group;
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method:
Task: set domain values
Programing:
02/2004 OK Implementation
10/2004 CMCD Gradient
02/2005 MB Water Depth // ToDo
03/2005 WW By patches of domain //ToDo
06/2005 OK OVERLAND_FLOW
12/2005 OK RESTART
07/2008 WW Patch-wise polynomal distribution
**************************************************************************/
void CInitialCondition::SetDomain(int nidx)
{
   int k, onZ;
   long i;
   double node_val, node_depth;
   vector<long> nodes_vector;
   CFEMesh* m_msh = this->getProcess()->m_msh;
   onZ = 0;
   if (m_msh)                                     //OK
      onZ = m_msh->GetCoordinateFlag() % 10;
   node_depth = 0.0;
   k = 0;
   if (SubNumber == 0)
   {
      if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
      {
         // if (this->getProcess()->pcs_type_name.compare("OVERLAND_FLOW") == 0)
         if (this->getProcess()->getProcessType() == OVERLAND_FLOW)
         {
                                                  //OK MSH
            for (i = 0; i < this->getProcess()->m_msh->GetNodesNumber(false); i++)
            {
               node_val = node_value_vector[0]->node_value
                  + this->getProcess()->m_msh->nod_vector[i]->Z();
               this->getProcess()->SetNodeValue(i, nidx, node_val);
            }
         }
         else
         {
            //................................................................
            node_val = node_value_vector[0]->node_value;
                                                  //OK MSH
            for (i = 0; i < m_msh->GetNodesNumber(true); i++)
               this->getProcess()->SetNodeValue(i, nidx, node_val);
            //................................................................
         }
      }
      //--------------------------------------------------------------------
                                                  // Remove unused stuff by WW
      if (this->getProcessDistributionType() == FiniteElement::GRADIENT)
      {
                                                  //WW
         for (i = 0; i < m_msh->GetNodesNumber(true); i++)
         {
            if (onZ == 1)                         //2D
               node_depth = m_msh->nod_vector[i]->Y();
            if (onZ == 2)                         //3D
               node_depth = m_msh->nod_vector[i]->Z();
            node_val = gradient_ref_depth_gradient * (gradient_ref_depth
               - node_depth) + gradient_ref_depth_value;
            this->getProcess()->SetNodeValue(
               m_msh->nod_vector[i]->GetIndex(), nidx, node_val);
         }
      }                                           //if(dis_type_name.find("GRADIENT")!=string::npos)
      //----------------------------------------------------------------------
      if (this->getProcessDistributionType() == FiniteElement::RESTART)
      {
         char line[MAX_ZEILEN];
         int no_var;
         int* var_n;
         std::vector<std::string> var_name;
         std::string var_name_string;
         std::string sdummy;
         double ddummy, dddummy;
         long ldummy;
         //....................................................................
         // File handling
         if (rfr_file_name.size() == 0)
         {
            cout << "Warning in CInitialCondition::SetDomain - no RFR file"
               << endl;
            return;
         }
         std::ifstream rfr_file;
         CGSProject* m_gsp = GSPGetMember("msh");
         std::string restart_file_name;
         restart_file_name = rfr_file_name;
         if (m_gsp)
            restart_file_name = m_gsp->path + rfr_file_name;
         //---------------------------------------------------------------------
         else                                     //Get absolut path of the file. 07.01.2009. WW
         {
            basic_string<char>::size_type indexChWin, indexChLinux;
            indexChWin = indexChLinux = 0;
            indexChWin = FileName.find_last_of('\\');
            indexChLinux = FileName.find_last_of('/');
            //
            string funfname;
            if (indexChWin == string::npos && indexChLinux == string::npos)
               funfname = rfr_file_name;
            else if (indexChWin != string::npos)
            {
               funfname = FileName.substr(0, indexChWin);
               funfname = funfname + "\\" + rfr_file_name;
            }
            else if (indexChLinux != string::npos)
            {
               funfname = FileName.substr(0, indexChLinux);
               funfname = funfname + "/" + rfr_file_name;
            }
            restart_file_name = funfname;
         }
         //-------------------------------------------------------------------
         rfr_file.open(restart_file_name.c_str(), ios::in);
         if (!rfr_file.good())
         {
            cout << "Warning in CInitialCondition::SetDomain - no RFR file"
               << endl;
            return;
         }
         rfr_file.seekg(0L, ios::beg);            // rewind?
         rfr_file.getline(line, MAX_ZEILEN);      //#0#0#0#1#100000#0#4.2.13 ###########################################
         rfr_file.getline(line, MAX_ZEILEN);      //1 1 4
         rfr_file >> no_var;                      //2 1 1
         var_n = new int[no_var];
         for (i = 0; i < no_var; i++)
         {
            rfr_file >> var_n[i];
         }
         size_t pos;
         for (i = 0; i < no_var; i++)             // HEAD, m ...
         {
            rfr_file >> var_name_string >> sdummy;
            pos = var_name_string.find_first_of(',');
            var_name_string.erase(pos, 1);
            var_name.push_back(var_name_string);
         }
         while (!rfr_file.eof())                  //TODO Implementation/Manual Correction of Keyword #STOP: TK
         {
            rfr_file >> dddummy;
            ldummy = (long) dddummy;
            //WW        cout << ldummy << endl;
            for (i = 0; i < no_var; i++)          // HEAD, m ...
            {
               rfr_file >> ddummy;
               this->getProcess()->SetNodeValue(ldummy, nidx, ddummy);
            }
         }
         //....................................................................
         delete var_n;
         var_name.clear();
      }
      //----------------------------------------------------------------------
   }
   //========================================================================
   else                                           //WW
   {
      bool quadratic = false;
      /// In case of P_U coupling monolithic scheme
      if (this->getProcess()->type == 41)         //WW Mono
      {
         if (convertPrimaryVariableToString(
            this->getProcessPrimaryVariable()).find("DISPLACEMENT")
            != string::npos)                      //Deform
            quadratic = true;
         else
            quadratic = false;
      } else if (this->getProcess()->type == 4)
      quadratic = true;
      else
         quadratic = false;
      //WW if (m_msh){
      for (k = 0; (size_t) k < SubNumber; k++)
      {
         GEOGetNodesInMaterialDomain(m_msh, subdom_index[k], nodes_vector,
            quadratic);
         if (this->getProcessDistributionType() == FiniteElement::GRADIENT)
         {
            if (k == 0)                           //TEST for DECOVALEX
            {
               for (i = 0; i < (int) nodes_vector.size(); i++)
                  this->getProcess()->SetNodeValue(nodes_vector[i], nidx,
                     subdom_ic[k]);
            }
            else
            {
               for (i = 0; i < (int) nodes_vector.size(); i++)
               {
                  if (onZ == 1)                   //2D
                     node_depth
                        = m_msh->nod_vector[nodes_vector[i]]->Y();
                  if (onZ == 2)                   //2D
                     node_depth
                        = m_msh->nod_vector[nodes_vector[i]]->Z();
                  node_val = ((gradient_ref_depth_gradient)
                     * (gradient_ref_depth - node_depth))
                     + gradient_ref_depth_value;
                  this->getProcess()->SetNodeValue(nodes_vector[i], nidx,
                     node_val);
               }
            }
         }
         /* KR not used
          else if(dis_type_name.find("FUNCTION")!=string::npos) //01.07.2008 WW
          {
          for(i=0;i<(int)nodes_vector.size();i++)
          {
          CNode *thisNode = m_msh->nod_vector[nodes_vector[i]];
          this->getProcess()->SetNodeValue(nodes_vector[i],nidx,DistributionFuntion(k, thisNode->X(),thisNode->Y(), thisNode->Z()));
          }
          }
          */
         else
         {
            for (i = 0; i < (int) nodes_vector.size(); i++)
               this->getProcess()->SetNodeValue(nodes_vector[i], nidx,
                  subdom_ic[k]);
         }
      }
      //  }
      /*  // Comment by WW
       else {
       for(k=0; k<SubNumber; k++)
       {
       GEOGetNodesInMaterialDomain(subdom_index[k], nodes_vector);
       if(dis_type_name.find("GRADIENT")!=string::npos)
       {
       if(k==0) //TEST for DECOVALEX
       {
       for(i=0;i<(int)nodes_vector.size();i++)
       SetNodeVal(nodes_vector[i],nidx, subdom_ic[k]);
      }
      else
      for(i=0;i<(int)nodes_vector.size();i++)
      {
      if(max_dim==1) //2D
      node_depth = GetNodeY(nodes_vector[i]);
      if(max_dim==2) //3D
      node_depth = GetNodeZ(nodes_vector[i]);
      node_val = ((gradient_ref_depth_gradient)*(gradient_ref_depth-node_depth))+
      gradient_ref_depth_value;
      SetNodeVal(nodes_vector[i],nidx,node_val);
      }

      }
      else
      {
      for(i=0;i<(int)nodes_vector.size();i++)
      SetNodeVal(nodes_vector[i],nidx, subdom_ic[k]);
      }

      }
      }
      */
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
09/2010 KR cleaned up code
**************************************************************************/
void ICDelete()
{
   size_t no_ic = ic_vector.size();
   for(size_t i=0;i<no_ic;i++)
      delete ic_vector[i];
   ic_vector.clear();

   no_ic = ic_group_vector.size();
   for(size_t i=0;i<no_ic;i++)
      delete ic_group_vector[i];
   ic_group_vector.clear();
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
09/2010 KR cleaned up code
**************************************************************************/
void ICGroupDelete()
{
   size_t no_ic_group = ic_group_vector.size();
   for (size_t i=0; i<no_ic_group; i++)
      delete ic_group_vector[i];
   ic_group_vector.clear();
}


/**************************************************************************
FEMLib-Method:
Task: set domain values to elements
Programing:
06/2006 MX Implementation
**************************************************************************/
void CInitialCondition::SetDomainEle(int nidx)
{
   int k;
   long i;
   double ele_val=0.0;
   std::vector<long>ele_vector;
   CFEMesh* m_msh = this->getProcess()->m_msh;
   k=0;
   CElem* m_ele = NULL;

   if(SubNumber==0)                               //only for constant values
   {
      ele_val = node_value_vector[0]->node_value;
      for(i=0;i<(long)m_msh->ele_vector.size();i++)
      {
         m_ele = m_msh->ele_vector[i];
         if(m_ele->GetMark())                     // Marked for use
         {
            this->getProcess()->SetElementValue(i,nidx, ele_val);
         }
      }
   }
   //========================================================================
   else                                           //MX
   {
      if (m_msh)
      {
         for(k=0; (size_t)k<SubNumber; k++)
         {
            for(i=0;i<(long)m_msh->ele_vector.size();i++)
            {
               m_ele = m_msh->ele_vector[i];
                                                  // Marked for use
               if(m_ele->GetMark() && m_ele->GetPatchIndex() == subdom_index[k] )
               {
                  this->getProcess()->SetElementValue(i, nidx, subdom_ic[k]);
                  ele_val = this->getProcess()->GetElementValue(i, nidx);
               }
            }
         }
      }                                           //if
   }
}


/**************************************************************************
FEMLib-Method:
07/2007 OK Implementation
**************************************************************************/
CInitialCondition* ICGet(string ic_name)
{
   CInitialCondition *m_ic = NULL;
   for(int i=0;i<(int)ic_vector.size();i++)
   {
      m_ic = ic_vector[i];
      if(m_ic->getProcessType() == convertProcessType(ic_name))
         return m_ic;
   }
   return NULL;
}
