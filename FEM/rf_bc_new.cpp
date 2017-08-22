/**************************************************************************
 FEMLib - Class: BC BoundaryConditions
 Task:
 Programing:
 02/2004 OK Implementation
 last modified
 **************************************************************************/
#include "makros.h"
// C++ STL
#include <cmath>
#include <cfloat>
#include <iostream>

// GEOLib
//#include "geo_lib.h"
//#include "geo_sfc.h"
#include "files0.h"

// GEOLIB
#include "GEOObjects.h"

// MSHLib
//#include "mshlib.h"
// FEMLib
extern void remove_white_space(std::string*);
//#include "problem.h"
#include "gs_project.h"
#include "tools.h"
//#include "rf_node.h"
#include "rf_bc_new.h"
//#include "rf_pcs.h"
//#include "rf_fct.h"
#include "rfmat_cp.h"
//#include "geo_ply.h"
// MathLib
#include "mathlib.h"

CBoundaryConditionNode::CBoundaryConditionNode()
{
   conditional = false;
}


//==========================================================================
std::list<CBoundaryCondition*> bc_list;
std::vector<std::string> bc_db_head;
std::list<CBoundaryConditionsGroup*> bc_group_list;
std::vector<CBoundaryCondition*> bc_db_vector;

/**************************************************************************
 FEMLib-Method:
 Task: BC constructor
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
CBoundaryCondition::CBoundaryCondition() :
GeoInfo (), geo_name (""), _curve_index (-1)
{
   this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
   // FCT
   conditional = false;
   time_dep_interpol = false;
   epsilon = 1e-9;                                //NW

}


/**************************************************************************
 FEMLib-Method:
 Task: BC deconstructor
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
CBoundaryCondition::~CBoundaryCondition()
{
   // DIS
   node_number_vector.clear();
   geo_node_number = -1;
   geo_node_value = 0.0;
}


const std::string& CBoundaryCondition::getGeoName () const
{
   return geo_name;
}


/**************************************************************************
 FEMLib-Method:
 Task: BC read function
 Programing:
 01/2004 OK Implementation
 09/2004 OK POINTS method
 11/2004 MX stream string
 **************************************************************************/
std::ios::pos_type CBoundaryCondition::Read(std::ifstream *bc_file,
const GEOLIB::GEOObjects& geo_obj, const std::string& unique_fname)
{
   std::string line_string;
   bool new_keyword = false;
   std::ios::pos_type position;

   std::string sub_string, strbuff;
   int ibuff;                                     //pos,
   double dbuff;                                  //WW
   std::stringstream in;

   //========================================================================
   // Schleife ueber alle Phasen bzw. Komponenten
   while (!new_keyword)
   {
      position = bc_file->tellg();
      line_string = GetLineFromFile1(bc_file);
      if (line_string.size() < 1)
         break;
      if (line_string.find("#") != std::string::npos)
      {
         new_keyword = true;
         break;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$PCS_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));
         std::string tmp;
         in >> tmp;                               // pcs_type_name;
         this->setProcessType(convertProcessType (tmp));
         in.clear();
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));
         std::string tmp;
         in >> tmp;                               // _pcs_pv_name;
         if ( this->_pcs_type == MASS_TRANSPORT )
         {
             // HS set the pointer to MCP based on component name.
             // a check whether this name is existing and unique. 
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
         }
         else
         {
             setProcess( PCSGet( this->getProcessType() ) );
             setProcessPrimaryVariable (convertPrimaryVariable (tmp));
         }
         in.clear();
      }
      //....................................................................
      // HS, this is new. later on we should stick to COMP_NAME, PRIMARY_VARIABLE support will be removed. 
      if (line_string.find("$COMP_NAME") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));
         std::string tmp;
         in >> tmp;                               // _pcs_pv_name;
         if ( this->_pcs_type == MASS_TRANSPORT )
             // HS set the pointer to MCP based on component name.
             // check whether this name is existing and unique. 
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
      }
      //....................................................................
                                                  //subkeyword found
      if (line_string.find("$GEO_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));
         in >> sub_string;                        //sub_line

         if (sub_string.compare("POINT") == 0)    //OK
         {
            in >> geo_name;                       //sub_line

            // TF 05/2010 - get the point from the point vector
            const GEOLIB::Point *pnt ((geo_obj.getPointVecObj(unique_fname))->getElementByName (geo_name));
            if (pnt == NULL)
            {
               std::cerr << "ERROR in CBoundaryCondition::Read() - Point name \"" << geo_name << "\" from geometry \"" << unique_fname << "\" not found!" << std::endl;
               exit (1);
            }

            in.clear();
            geo_type_name = "POINT";
            setGeoType (GEOLIB::POINT);
            setGeoObj (pnt);
         }
         if (sub_string.find("POLYLINE") != std::string::npos)
         {
            in >> geo_name;                       //sub_line
            in.clear();
            // 07/2010 TF substitute source code in order to use new GEOLIB
            geo_type_name = "POLYLINE";
            //				CGLPolyline* m_ply = NULL;
            //				m_ply = GEOGetPLYByName(geo_name); //CC 10/05
            //				if (!m_ply)
            //					cout << "Warning in BCRead: no PLY data" << endl;

            // TF 07/2010 - get the index from the polyline vector
            if (geo_obj.getPolylineVecObj(unique_fname) != NULL)
            {
               const GEOLIB::Polyline* ply ((geo_obj.getPolylineVecObj(unique_fname))->getElementByName (geo_name));
               if (ply == NULL)
               {
                  std::cerr << "Warning CBoundaryCondition::Read: polyline name \"" << geo_name << "\" not found!" << std::endl;
                  exit (1);
               }
               setGeoObj (ply);
            }

            setGeoType (GEOLIB::POLYLINE);
         }

         if (sub_string.find("SURFACE") != std::string::npos)
         {
            in >> geo_name;                       //sub_line
            in.clear();
            // TF 07/2010 - get the index from the surface vector
            //				if (geo_obj.getSurfaceVecObj(unique_fname) != NULL) {
            //					if (!((geo_obj.getSurfaceVecObj(unique_fname))->getElementIDByName (geo_name, _geo_obj_idx))) {
            //						std::cerr << "Warning CBoundaryCondition::Read: surface name \"" << geo_name << "\" not found!" << std::endl;
            //						exit (1);
            //					}
            //				}
            geo_type_name = "SURFACE";
            setGeoType (GEOLIB::SURFACE);
         }
         if (sub_string.find("VOLUME") != std::string::npos)
         {
            in >> geo_name;                       //sub_line
            in.clear();

            //				geo_type_name = "VOLUME";
            setGeoType (GEOLIB::VOLUME);
         }
      }
      //....................................................................
                                                  //PCH
      if (line_string.find("$DIS_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));
         in >> line_string;                       //sub_line
         _periodic = false;                       // JOD

         // Soure terms are assign to element nodes directly. 23.02.2009. WW
         if (line_string.find("DIRECT") != std::string::npos) {
			 this->setProcessDistributionType(FiniteElement::DIRECT);
            in >> fname;
            fname = FilePath + fname;
            in.clear();
         }

		 /* KR not used
         // Patch-wise constant. 19.03.2009. WW
         if (line_string.find("PATCH_C") != string::npos) {
            dis_type_name = "PATCH_C";
            in >> fname;
            fname = FilePath + fname;
            in >> geo_node_value;
            in.clear();
         }
         */
         if (line_string.find("CONSTANT") != std::string::npos)
         {
            this->setProcessDistributionType(FiniteElement::CONSTANT);
            in >> geo_node_value;                 //sub_line
            in.clear();
         }
         if (line_string.find("LINEAR") != std::string::npos)
         {
            this->setProcessDistributionType(FiniteElement::LINEAR);
            // Distribued. WW
            int nLBC;
            in >> nLBC;                           //sub_line
            in.clear();

            //        sub_string = strtok(buffer,seps);
            //        sub_string = strtok( NULL, seps );
            //        int nLBC = atoi(sub_string.c_str());
            for (int i = 0; i < nLBC; i++)
            {
               in.str(GetLineFromFile1(bc_file));
               in >> ibuff >> dbuff >> strbuff;
               in.clear();

               //           *bc_file>>ibuff>>dbuff;
               _PointsHaveDistribedBC.push_back(ibuff);
               _DistribedBC.push_back(dbuff);
               if (strbuff.size() > 0)
               {
                  _PointsFCTNames.push_back(strbuff);
                  time_dep_interpol = true;
               }
            }
            //        bc_file->ignore(MAX_ZEILE,'\n');
         }
         /* KR not used
         if (line_string.find("PERIODIC") != string::npos) { // JOD
            dis_type_name = "PERIODIC";
            _periodic = true;
            dis_type = 0;
            in >> geo_node_value;
            in.clear();
            in.str(GetLineFromFile1(bc_file));
            in >> _periode_time_length >> _periode_phase_shift; //sub_line
            in.clear();
         } // subkeyword found
         */
      }
      // Time dependent function
      //..Time dependent curve ............................................
                                                  // subkeyword found
      if (line_string.find("$TIM_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));
         in >> line_string;

         //      bc_file->getline(buffer,MAX_ZEILE);
         //      line_string = buffer;
         //      pos1 = 0;
         //      sub_string = get_sub_string(buffer,delimiter_type,pos1,&pos2);
         if (line_string.find("CURVE") != std::string::npos)
         {
            //				tim_type_name = "CURVE";
            this->setProcessDistributionType(FiniteElement::CONSTANT);
            in >> _curve_index;
            in.clear();

            //        pos1=pos2+1;
            //        sub_string = get_sub_string(buffer,"  ",pos1,&pos2);
            //		_curve_index = atoi(sub_string.c_str());
         }
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$FCT_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));
         in >> fct_name;                          //sub_line
         in.clear();

         //		bc_file->getline(line,MAX_ZEILE);
         //      line_string = line;
         //      remove_white_space(&line_string);
         //      fct_name = line_string.substr(0);
      }
      //....................................................................
      //OK4105
                                                  //subkeyword found
      if (line_string.find("$MSH_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));
         in >> sub_string;                        //sub_line
         _msh_type_name = "NODE";
         if (sub_string.find("NODE") != std::string::npos)
         {
            in >> _msh_node_number;
            in.clear();
         }
      }
      //....................................................................
      //OK4108
                                                  // subkeyword found
      if (line_string.find("$DIS_TYPE_CONDITION") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));       // CONSTANT -21500.0
         in >> line_string;
         if (line_string.find("CONSTANT") != std::string::npos)
         {
            this->setProcessDistributionType(FiniteElement::CONSTANT);
            in >> geo_node_value;
            in.clear();
         }
         in.str(GetLineFromFile1(bc_file));       // 0.0 IF HEAD > 0.04
         std::string pcs_pv_name_cond;            // 07/2010 TF temp string
         in >> node_value_cond >> line_string >> pcs_pv_name_cond
            >> line_string >> condition;
         in.clear();
         in.str(GetLineFromFile1(bc_file));       // PCS OVERLAND_FLOW
         std::string pcs_type_name_cond;
         in >> line_string >> pcs_type_name_cond;
         in.clear();
         conditional = true;
      }
      //....................................................................
                                                  // NW
      if (line_string.find("$EPSILON") != std::string::npos)
      {
         in.str(GetLineFromFile1(bc_file));
         in >> epsilon;
         in.clear();
      }
      //....................................................................
   }
   return position;
}


/**************************************************************************
 FEMLib-Method: CBoundaryCondition::Write
 02/2004 OK Implementation
 07/2007 OK LINEAR
 10/2008 OK NOD
 06/2009 OK MSH_TYPE off
 **************************************************************************/
void CBoundaryCondition::Write(std::fstream* rfd_file) const
{
   //KEYWORD
   *rfd_file << "#BOUNDARY_CONDITION" << std::endl;
   //--------------------------------------------------------------------
   //NAME+NUMBER
   *rfd_file << " $PCS_TYPE" << std::endl;
   *rfd_file << "  " << convertProcessTypeToString(getProcessType()) << std::endl;
   *rfd_file << " $PRIMARY_VARIABLE" << std::endl;
   *rfd_file << "  " << convertPrimaryVariableToString(this->getProcessPrimaryVariable()) << std::endl;
   //--------------------------------------------------------------------
   //GEO_TYPE
   *rfd_file << " $GEO_TYPE" << std::endl;
   *rfd_file << "  ";
   *rfd_file << getGeoTypeAsString() << " " << geo_name << std::endl;

   //--------------------------------------------------------------------
   /*OK4910
    //MSH_TYPE
    if(msh_node_number>0){
    *rfd_file << " $MSH_TYPE" << endl;
    *rfd_file << "  ";
    *rfd_file << "NODE" << " " << msh_node_number << endl;
    }
    */
   //--------------------------------------------------------------------
   //DIS_TYPE
   *rfd_file << " $DIS_TYPE" << std::endl;
   *rfd_file << "  ";
   *rfd_file << convertDisTypeToString(this->getProcessDistributionType());
   //switch (dis_type_name[0]) {
   //case 'C': // Constant
   if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
   {
      *rfd_file << " " << geo_node_value;
      *rfd_file << std::endl;
      //break;
   }
   //case 'L': // Linear
   else if (this->getProcessDistributionType() == FiniteElement::LINEAR)
   {
      *rfd_file << " " << _PointsHaveDistribedBC.size() << std::endl;
      for (size_t i = 0; i < _PointsHaveDistribedBC.size(); i++)
      {
         *rfd_file << "  " << _PointsHaveDistribedBC[i] << " ";
         *rfd_file << "  " << _DistribedBC[i] << std::endl;
      }
      //break;
   }

   //FCT
   if (fct_name.length() > 0)                     //OK4108
   {
      *rfd_file << " $FCT_TYPE" << std::endl;
      *rfd_file << "  ";
      *rfd_file << fct_name << std::endl;
   }
}


/**************************************************************************
 FEMLib-Method: CBoundaryCondition::Write
 Task: write function
 Programing:
 02/2004 OK Implementation
 last modification:
 **************************************************************************/
void CBoundaryCondition::WriteTecplot(std::fstream* tec_file) const
{
   long i;
   CGLPolyline* m_polyline1 = NULL;
   CGLPolyline* m_polyline2 = NULL;
   // list<CGLPolyline*>::const_iterator p;
   std::vector<CGLPolyline*>::iterator p;
   Surface *m_surface = NULL;
   long no_points = 0;
   std::vector<CTriangle*> triangle_vector;

   *tec_file << "VARIABLES = X,Y,Z,V1" << std::endl;

   if (getGeoType () == GEOLIB::SURFACE)
   {
      m_surface = GEOGetSFCByName(geo_name);      //CC
      if (m_surface)
         switch (m_surface->type)
         {
            case 2:
               p = m_surface->polyline_of_surface_vector.begin();
               while (p != m_surface->polyline_of_surface_vector.end())
               {
                  m_polyline1 = *p;
                  ++p;
                  m_polyline2 = *p;
                  break;
               }
            no_points = (long) m_polyline1->point_vector.size();
            /*
             for(i=0;i<no_points-1;i++) {
             m_triangle = new CTriangle;
             m_triangle->x[0] = m_polyline1->point_vector[i]->x;
             m_triangle->y[0] = m_polyline1->point_vector[i]->y;
             m_triangle->z[0] = m_polyline1->point_vector[i]->z;
             m_triangle->x[1] = m_polyline1->point_vector[i+1]->x;
             m_triangle->y[1] = m_polyline1->point_vector[i+1]->y;
             m_triangle->z[1] = m_polyline1->point_vector[i+1]->z;
             m_triangle->x[2] = m_polyline2->point_vector[i+1]->x;
             m_triangle->y[2] = m_polyline2->point_vector[i+1]->y;
            m_triangle->z[2] = m_polyline2->point_vector[i+1]->z;
            triangle_vector.push_back(m_triangle);
            m_triangle = new CTriangle;
            m_triangle->x[0] = m_polyline2->point_vector[i]->x;
            m_triangle->y[0] = m_polyline2->point_vector[i]->y;
            m_triangle->z[0] = m_polyline2->point_vector[i]->z;
            m_triangle->x[1] = m_polyline2->point_vector[i+1]->x;
            m_triangle->y[1] = m_polyline2->point_vector[i+1]->y;
            m_triangle->z[1] = m_polyline2->point_vector[i+1]->z;
            m_triangle->x[2] = m_polyline1->point_vector[i+1]->x;
            m_triangle->y[2] = m_polyline1->point_vector[i+1]->y;
            m_triangle->z[2] = m_polyline1->point_vector[i+1]->z;
            triangle_vector.push_back(m_triangle);
            }
            */
            break;
         }
      }

      long no_nodes = 2* no_points ;
      //long no_elements = triangle_vector.size();
      long no_elements = 2* (no_points -1);
      // Write
      *tec_file << "ZONE T = " << geo_name << ", "
         << "N = " << no_nodes << ", "
         << "E = " << no_elements << ", "
         << "F = FEPOINT" << ", " << "ET = TRIANGLE" << std::endl;
      if(m_polyline1)
         for(i=0;i<no_points;i++)
      {
         *tec_file
            << m_polyline1->point_vector[i]->x << " " << m_polyline1->point_vector[i]->y << " " << m_polyline1->point_vector[i]->z << " " \
            << geo_node_value << std::endl;
      }
      if(m_polyline2)
         for(i=0;i<no_points;i++)
      {
         *tec_file \
            << m_polyline2->point_vector[i]->x << " " << m_polyline2->point_vector[i]->y << " " << m_polyline2->point_vector[i]->z << " " \
            << geo_node_value << std::endl;
      }
      for(i=0;i<no_points-1;i++)
      {
         *tec_file \
            << i+1 << " " << i+1+1 << " " << no_points+i+1 << std::endl;
      }
      for(i=0;i<no_points-1;i++)
      {
         *tec_file \
            << no_points+i+1 << " " << no_points+i+1+1 << " " << i+1+1 << std::endl;
      }
   }


/**************************************************************************
 FEMLib-Method:
 Task: BC read function
 Programing:
 01/2004 OK Implementation
 01/2005 OK Boolean type
 01/2005 OK Destruct before read
 05/2010 TF changes due to new GEOLIB integration, some improvements
 **************************************************************************/
   bool BCRead(std::string file_base_name, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
   {
      char line[MAX_ZEILE];
      std::string line_string, bc_file_name;

      // File handling
      bc_file_name = file_base_name + BC_FILE_EXTENSION;

      std::ifstream bc_file(bc_file_name.data(), std::ios::in);
      if (!bc_file.good())
      {
         std::cout << "! Error in BCRead: No boundary conditions !" << std::endl;
         return false;
      }

      // Keyword loop
      std::cout << "BCRead ... " << std::flush;
      while (!bc_file.eof())
      {
         bc_file.getline(line, MAX_ZEILE);
         line_string = line;
         if (line_string.find("#STOP") != std::string::npos)
         {
            std::cout << "done, read " << bc_list.size() << " boundary conditions" << std::endl;
            return true;
         }
         if (line_string.find("#BOUNDARY_CONDITION") != std::string::npos)
         {
            CBoundaryCondition *bc (new CBoundaryCondition());
            std::ios::pos_type position = bc->Read (&bc_file, geo_obj, unique_name);
            bc_list.push_back(bc);
            bc_file.seekg(position, std::ios::beg);
         }                                        // keyword found
      }                                           // eof
      std::cout << "done, read " << bc_list.size() << " boundary conditions" << std::endl;
      return true;
   }


/**************************************************************************
 FEMLib-Method: BCWrite
 Task: master write function
 Programing:
 02/2004 OK Implementation
 last modification:
 **************************************************************************/

   void BCWrite(std::string base_file_name)
   {
      CBoundaryCondition *m_bc = NULL;
      std::string sub_line;
      std::string line_string;
      //========================================================================
      // File handling
      std::string bc_file_name = base_file_name + BC_FILE_EXTENSION;
      std::fstream bc_file(bc_file_name.data(), std::ios::trunc | std::ios::out);
      bc_file.setf(std::ios::scientific, std::ios::floatfield);
      bc_file.precision(12);
      //OK string tec_file_name = base_file_name + ".tec";
      //OK fstream tec_file (tec_file_name.data(),ios::trunc|ios::out);
      //OK tec_file.setf(ios::scientific,ios::floatfield);
      //OK tec_file.precision(12);
      if (!bc_file.good())
         return;
      bc_file.seekg(0L, std::ios::beg);           // rewind?
      //========================================================================
      bc_file
         << "GeoSys-BC: Boundary Conditions ------------------------------------------------\n";
      //========================================================================
      // BC list
      std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
      while (p_bc != bc_list.end())
      {
         m_bc = *p_bc;
         m_bc->Write(&bc_file);
         //OK m_bc->WriteTecplot(&tec_file);
         ++p_bc;
      }
      bc_file << "#STOP";
      bc_file.close();
      //OK tec_file.close();
   }


/**************************************************************************
 FEMLib-Method:
 01/2004 OK Implementation
 07/2007 OK V2, global function
 **************************************************************************/
//CBoundaryCondition* BCGet(const std::string &pcs_name, const std::string &geo_type_name,
//		const std::string &geo_name)
//{
//	std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
//	while (p_bc != bc_list.end()) {
//		if (((*p_bc)->pcs_type_name.compare(pcs_name) == 0)
//				&& ((*p_bc)->geo_type_name.compare(geo_type_name) == 0)
//				&& ((*p_bc)->getGeoName().compare(geo_name) == 0))
//			return *p_bc;
//		++p_bc;
//	}
//	return NULL;
//}

/**************************************************************************
 GeoSys source term function:
 02/2009 WW Implementation
 **************************************************************************/
   inline void CBoundaryCondition::DirectAssign(long ShiftInNodeVector)
   {
      std::string line_string;
      std::stringstream in;
      long n_index;
      double n_val;
      CRFProcess* m_pcs = NULL;
      CBoundaryConditionNode *m_node_value = NULL;

      m_pcs = PCSGet(convertProcessTypeToString(this->getProcessType()));

      //========================================================================
      // File handling
      std::ifstream d_file(fname.c_str(), std::ios::in);
      //if (!st_file.good()) return;

      if (!d_file.good())
      {
         std::cout << "! Error in direct node source terms: Could not find file:!\n"
            << fname << std::endl;
         abort();
      }
      // Rewind the file
      d_file.clear();
      d_file.seekg(0L, std::ios::beg);
      //========================================================================
      while (!d_file.eof())
      {
         line_string = GetLineFromFile1(&d_file);
         if (line_string.find("#STOP") != std::string::npos)
            break;

         in.str(line_string);
         in >> n_index >> n_val;
         in.clear();
         //
         m_node_value = new CBoundaryConditionNode;
         m_node_value->conditional = false;
         m_node_value->msh_node_number = n_index + ShiftInNodeVector;
         m_node_value->geo_node_number = n_index;
         m_node_value->node_value = n_val;
         m_node_value->CurveIndex = _curve_index;
         m_pcs->bc_node.push_back(this);
         m_pcs->bc_node_value.push_back(m_node_value);
      }                                           // eof
   }


/**************************************************************************
 GeoSys BC function:
 03/2009 WW Implementation
 **************************************************************************/
   inline void CBoundaryCondition::PatchAssign(long ShiftInNodeVector)
   {
      std::string line_string;
      std::string st_file_name;
      std::stringstream in;
      long n_index;
      std::vector<long> sfc_nodes;
      CRFProcess* m_pcs = NULL;
      CBoundaryConditionNode *m_node_value = NULL;

      m_pcs = PCSGet(convertProcessTypeToString(this->getProcessType()));

      Surface *m_surface = NULL;
      m_surface = GEOGetSFCByName(geo_name);
      //========================================================================
      // File handling
      std::ifstream d_file(fname.c_str(), std::ios::in);
      //if (!st_file.good()) return;

      if (!d_file.good())
      {
         std::cout << "! Error in direct node source terms: Could not find file:!\n"
            << fname << std::endl;
         abort();
      }
      // Rewind the file
      d_file.clear();
      d_file.seekg(0L, std::ios::beg);
      //========================================================================
      while (!d_file.eof())
      {
         line_string = GetLineFromFile1(&d_file);
         if (line_string.find("#STOP") != std::string::npos)
            break;

         in.str(line_string);
         in >> n_index;
         in.clear();
         sfc_nodes.push_back(n_index);
      }
      if (m_surface)
         m_pcs->m_msh->GetNODOnSFC_PLY_XY(m_surface, sfc_nodes, true);
      for (long i = 0; i < (long) sfc_nodes.size(); i++)
      {
         //
         m_node_value = new CBoundaryConditionNode;
         m_node_value->conditional = false;
         n_index = sfc_nodes[i];
         m_node_value->msh_node_number = n_index + ShiftInNodeVector;
         m_node_value->geo_node_number = n_index;
         m_node_value->node_value = geo_node_value;
         m_node_value->CurveIndex = _curve_index;
         m_pcs->bc_node.push_back(this);
         m_pcs->bc_node_value.push_back(m_node_value);
      }                                           // eof
   }


   CBoundaryConditionsGroup::CBoundaryConditionsGroup(void)
   {
      msh_node_number_subst = -1;                 //
      time_dep_bc = -1;
   }


   CBoundaryConditionsGroup::~CBoundaryConditionsGroup(void)
   {
      /*
       int group_vector_length = group_vector.size();
       int i;
       for(i=0;i<group_vector_length;i++)
       group_vector.pop_back();
       */
      //  group_vector.clear();

   }


/**************************************************************************
 FEMLib-Method: CBoundaryCondition::Set
 Task: set boundary conditions
 Programing:
 02/2004 OK Implementation
 09/2004 WW Interpolation of piecewise linear BC
 02/2005 OK MSH types
 03/2005 OK MultiMSH, PNT
 08/2005 WW Changes due to the new geometry finite element.
 12/2005 OK FCT
 04/2006 WW New storage
09/2006 WW Move linear interpolation to new MSH structure
12/2007 WW Linear distributed BC in a surface
10/2008 WW/CB SetTransientBCtoNodes
last modification:
**************************************************************************/
   void CBoundaryConditionsGroup::Set(CRFProcess* m_pcs, int ShiftInNodeVector, const std::string& this_pv_name)
   {
      //	long number_of_nodes = 0;
      long *nodes = NULL;
      std::vector<long> nodes_vector;
      std::vector<double> node_value;
      CGLPolyline *m_polyline = NULL;
      long i, j;
      CBoundaryCondition *m_bc = NULL;
      CBoundaryConditionNode *m_node_value = NULL;
      group_name = _pcs_type_name;
      bool quadratic = false;
      bool cont = false;                          //WW
      //----------------------------------------------------------------------
      if (!this_pv_name.empty())
         _pcs_pv_name = this_pv_name;
      CFEMesh* m_msh = m_pcs->m_msh;
      // Tests //OK

      if (!m_msh)
      {
         std::cout << "Warning in CBoundaryConditionsGroup::Set - no MSH data"
            << std::endl;
         //return;
      }
      if (m_msh)                                  //WW
      {
         /// In case of P_U coupling monolithic scheme
         if (m_pcs->type == 41)                   //WW Mono
         {
                                                  //Deform
            if (_pcs_pv_name.find("DISPLACEMENT") != std::string::npos || _pcs_pv_name.find("VELOCITY_DM") != std::string::npos)
               quadratic = true;
            else
               quadratic = false;
         }
         else if (m_pcs->type == 4)
            quadratic = true;
         else
            quadratic = false;
         m_pcs->m_msh->SwitchOnQuadraticNodes(quadratic);
      }

      std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
      while (p_bc != bc_list.end())
      {
         m_bc = *p_bc;
         if (m_bc->time_dep_interpol)             //WW/CB
         {
            ++p_bc;
            continue;
         }

         if(m_bc->getProcessDistributionType()==FiniteElement::DIRECT)
         {
            m_bc->DirectAssign(ShiftInNodeVector);
            continue;
         }

         //====================================================================
         //OK if(m_bc->pcs_type_name.compare(pcs_type_name)==0){ //OK/SB 4108
         if ( ( m_bc->getProcess() == m_pcs ) && (m_bc->getProcessPrimaryVariable() == convertPrimaryVariable(_pcs_pv_name)) )
         {

            cont = false;
            //------------------------------------------------------------------
            if (m_bc->getGeoType () == GEOLIB::POINT)
            {
               m_node_value = new CBoundaryConditionNode;
               // Get MSH node number
               if (m_bc->getProcessDistributionType() == FiniteElement::CONSTANT)
               {
                  // KR not used || m_bc->dis_type_name.compare("PERIODIC") == 0) { //JOD
                  // TF - tests from CC removed -> should be checked while reading data
                  m_node_value->geo_node_number = m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(m_bc->getGeoObj()));
               }

               m_node_value->conditional = cont;
               m_node_value->CurveIndex = m_bc->getCurveIndex();
               m_node_value->node_value = m_bc->geo_node_value;
               m_node_value->msh_node_number = m_node_value->geo_node_number
                  + ShiftInNodeVector;            //WW
                                                  //YD/WW
               m_node_value->pcs_pv_name = _pcs_pv_name;
               m_node_value->msh_node_number_subst = msh_node_number_subst;
               m_pcs->bc_node.push_back(m_bc);    //WW
                                                  //WW
               m_pcs->bc_node_value.push_back(m_node_value);
         }
         //------------------------------------------------------------------
         if (m_bc->getGeoType () == GEOLIB::POLYLINE)
         {
                                                  //CC
            m_polyline = GEOGetPLYByName(m_bc->geo_name);

            // 08/2010 TF get the new polyline data structure
            const GEOLIB::Polyline *ply (static_cast<const GEOLIB::Polyline*> (m_bc->getGeoObj()));

            if (m_polyline)
            {
               if (m_bc->getProcessDistributionType() == FiniteElement::CONSTANT)
               {
                  // KR not used || m_bc->dis_type_name.compare("PERIODIC") == 0) { // JOD

                  //						std::vector<long> long_nodes_vector;
                  //						m_msh->GetNODOnPLY(m_polyline, long_nodes_vector);

                  // 08/2010 TF
                  double msh_min_edge_length = m_msh->getMinEdgeLength();
                  m_msh->setMinEdgeLength (m_polyline->epsilon);
                  std::vector<size_t> my_nodes_vector;
                  m_msh->GetNODOnPLY (ply, my_nodes_vector);
                  m_msh->setMinEdgeLength (msh_min_edge_length);

                  nodes_vector.clear ();
                  for (size_t k(0); k<my_nodes_vector.size(); k++)
                     nodes_vector.push_back (my_nodes_vector[k]);

                  // for some benchmarks we need the vector entries sorted by index
                  std::vector<size_t> perm;
                  for (size_t k(0); k<my_nodes_vector.size(); k++) perm.push_back (k);
                  Quicksort<long> (nodes_vector, 0, nodes_vector.size(), perm);

                  for (i = 0; i < (long) nodes_vector.size(); i++)
                  {
                     m_node_value = new CBoundaryConditionNode();
                     m_node_value->msh_node_number = -1;
                     m_node_value->msh_node_number = nodes_vector[i]
                        + ShiftInNodeVector;
                     m_node_value->geo_node_number = nodes_vector[i];
                                                  //dis_prop[0];
                     m_node_value->node_value = m_bc->geo_node_value;
                     m_node_value->CurveIndex = m_bc->getCurveIndex();
                                                  //YD/WW
                     m_node_value->pcs_pv_name = _pcs_pv_name;
                                                  //WW
                     m_pcs->bc_node.push_back(m_bc);
                                                  //WW
                     m_pcs->bc_node_value.push_back(m_node_value);
                     //WW group_vector.push_back(m_node_value);
                     //WW bc_group_vector.push_back(m_bc); //OK
                  }
               }
               //................................................................
                                                  //WW
               if (m_bc->getProcessDistributionType() == FiniteElement::LINEAR)
               {
                  //						m_msh->GetNODOnPLY(m_polyline, nodes_vector);
                  double msh_min_edge_length = m_msh->getMinEdgeLength();
                  m_msh->setMinEdgeLength (m_polyline->epsilon);
                  std::vector<size_t> my_nodes_vector;
                  m_msh->GetNODOnPLY (ply, my_nodes_vector);
                  m_msh->setMinEdgeLength (msh_min_edge_length);

                  nodes_vector.clear ();
                  for (size_t k(0); k<my_nodes_vector.size(); k++)
                     nodes_vector.push_back (my_nodes_vector[k]);

                  // for some benchmarks we need the vector entries sorted by index
                  std::vector<size_t> perm;
                  for (size_t k(0); k<my_nodes_vector.size(); k++) perm.push_back (k);
                  Quicksort<long> (nodes_vector, 0, nodes_vector.size(), perm);

                  //            nodes = MSHGetNodesClose(&number_of_nodes, m_polyline);//CC
                  node_value.resize(nodes_vector.size());
                  // Piecewise linear distributed. WW
                  for (i = 0; i < (int) m_bc->getDistribedBC().size(); i++)
                  {
                     for (j = 0; j
                        < (int) m_polyline->point_vector.size(); j++)
                     {
                        if (m_bc->getPointsWithDistribedBC()[i]
                           == m_polyline->point_vector[j]->id)
                        {
                           if (fabs(m_bc->getDistribedBC()[i])
                              < MKleinsteZahl)
                              m_bc->getDistribedBC()[i] = 1.0e-20;
                           m_polyline->point_vector[j]->setPropert (m_bc->getDistribedBC()[i]);
                           break;
                        }
                     }
                  }
                  InterpolationAlongPolyline(m_polyline, node_value);
                  for (i = 0; i < (int) nodes_vector.size(); i++)
                  {
                     m_node_value = new CBoundaryConditionNode();
                     m_node_value->msh_node_number = -1;
                     m_node_value->msh_node_number = nodes_vector[i]
                        + ShiftInNodeVector;
                     m_node_value->geo_node_number = nodes_vector[i];
                     m_node_value->node_value = node_value[i];
                                                  //YD/WW
                     m_node_value->pcs_pv_name = _pcs_pv_name;
                     m_node_value->CurveIndex = m_bc->getCurveIndex();
                                                  //WW
                     m_pcs->bc_node.push_back(m_bc);
                                                  //WW
                     m_pcs->bc_node_value.push_back(m_node_value);
                     //WW group_vector.push_back(m_node_value);
                     //WW bc_group_vector.push_back(m_bc); //OK
                  }
                  node_value.clear();
               }
               //................................................................
               //delete(values);
               Free(nodes);
            }                                     // if(m_ply)
         }
         //------------------------------------------------------------------
         if (m_bc->getGeoType () == GEOLIB::SURFACE)
         {
                                                  //CC10/05
            Surface *m_surface = GEOGetSFCByName(m_bc->geo_name);
            if (m_surface)
            {
               //..............................................................
               if (m_msh)                         //MSH
                  m_msh->GetNODOnSFC(m_surface, nodes_vector);
               long nodes_vector_length = (long) nodes_vector.size();
               //..............................................................
               if (m_bc->getProcessDistributionType() == FiniteElement::LINEAR)
               {
                  // Linear interpolation  polygon-wise. WW
                  // list<CGLPolyline*>::const_iterator p = m_surface->polyline_of_surface_list.begin();
                  std::vector<CGLPolyline*>::iterator p =
                     m_surface->polyline_of_surface_vector.begin();
                  node_value.resize(nodes_vector_length);
                  p = m_surface->polyline_of_surface_vector.begin();
                  while (p != m_surface->polyline_of_surface_vector.end())
                  {
                     m_polyline = *p;
                     for (i = 0; i < (int) m_bc->getDistribedBC().size(); i++)
                     {
                        for (j = 0; j
                           < (int) m_polyline->point_vector.size(); j++)
                        {
                           if (m_bc->getPointsWithDistribedBC ()[i]
                              == m_polyline->point_vector[j]->id)
                           {
                              if (fabs(m_bc->getDistribedBC()[i])
                                 < MKleinsteZahl)
                                 m_bc->getDistribedBC()[i] = 1.0e-20;
                              m_polyline->point_vector[j]->setPropert (m_bc->getDistribedBC()[i]);
                              break;
                           }
                        }
                     }
                     //InterpolationAlongPolyline(m_polyline, node_value);
                     p++;
                  }
                                                  //WW
                  node_value.resize(nodes_vector_length);
                  m_bc->SurfaceInterpolation(m_pcs, nodes_vector,
                     node_value);                 //WW
               }
               //..............................................................
               for (i = 0; i < nodes_vector_length; i++)
               {
                  m_node_value = new CBoundaryConditionNode();
                  m_node_value->msh_node_number = -1;
                  m_node_value->msh_node_number = nodes_vector[i]
                     + ShiftInNodeVector;         //nodes[i];

                  //nodes[i];
                  m_node_value->geo_node_number = nodes_vector[i];
                                                  //YD/WW
                  m_node_value->pcs_pv_name = _pcs_pv_name;
                                                  //WW
                  if (m_bc->getProcessDistributionType() == FiniteElement::LINEAR)

                     m_node_value->node_value = node_value[i];
                  else
                     m_node_value->node_value = m_bc->geo_node_value;
                  m_node_value->CurveIndex = m_bc->getCurveIndex();
                                                  //OK
                  m_bc->node_number_vector = nodes_vector;
                  m_pcs->bc_node.push_back(m_bc); //WW
                                                  //WW
                  m_pcs->bc_node_value.push_back(m_node_value);
                  //WW group_vector.push_back(m_node_value);
                  //WW bc_group_vector.push_back(m_bc); //OK
               }
               node_value.clear();
            }
         }
         //------------------------------------------------------------------
         // Material domain
         //			if (m_bc->geo_type_name.find("MATERIAL_DOMAIN") == 0) {
         //				GEOGetNodesInMaterialDomain(m_msh, m_bc->_geo_type,
         //						nodes_vector, quadratic);
         //				for (i = 0; i < (long) nodes_vector.size(); i++) {
         //					m_node_value = new CBoundaryConditionNode();
         //					m_node_value->msh_node_number = -1;
         //					m_node_value->msh_node_number = nodes_vector[i]
         //							+ ShiftInNodeVector; //nodes[i];
         //					m_node_value->geo_node_number = nodes_vector[i]; //nodes[i];
         //					m_node_value->node_value = m_bc->geo_node_value;
         //					m_node_value->pcs_pv_name = pcs_pv_name; //YD/WW
         //					m_node_value->CurveIndex = m_bc->getCurveIndex();
         //					m_pcs->bc_node.push_back(m_bc); //WW
         //					m_pcs->bc_node_value.push_back(m_node_value); //WW
         //					//WW group_vector.push_back(m_node_value);
         //					//WW bc_group_vector.push_back(m_bc); //OK
         //				}
         //			}
         //------------------------------------------------------------------
         // MSH types //OK4105
         if (m_bc->getMeshTypeName().compare("NODE") == 0)
         {
            m_node_value = new CBoundaryConditionNode;
            m_node_value->msh_node_number = m_bc->getMeshNodeNumber();
            m_node_value->geo_node_number = m_bc->getMeshNodeNumber();
            m_node_value->node_value = m_bc->geo_node_value;
            m_node_value->CurveIndex = m_bc->getCurveIndex();
                                                  //YD/WW
            m_node_value->pcs_pv_name = _pcs_pv_name;
            m_pcs->bc_node.push_back(m_bc);       //WW
                                                  //WW
            m_pcs->bc_node_value.push_back(m_node_value);
            //WW group_vector.push_back(m_node_value);
            //WW bc_group_vector.push_back(m_bc); //OK
         }
         //------------------------------------------------------------------
         // FCT types //OK
         if (m_bc->fct_name.size() > 0)
         {
            //WW
            for (i = 0; i < (long) m_pcs->bc_node_value.size(); i++)
            {
               m_pcs->bc_node_value[i]->fct_name = m_bc->fct_name;
               m_pcs->bc_node_value[i]->msh_node_number_subst
                  = msh_node_number_subst;
            }
            //WW fct_name = m_bc->fct_name;
         }
         //------------------------------------------------------------------
      }                                           // PCS
      ++p_bc;
   }                                              // list
   //======================================================================
   // SetTransientBCtoNodes  10/2008 WW/CB Implementation
   p_bc = bc_list.begin();
   while (p_bc != bc_list.end())
   {
      m_bc = *p_bc;
      if (!m_bc->time_dep_interpol)               //WW/CB
      {
         ++p_bc;
         continue;
      }
      //====================================================================
      //OK if(m_bc->pcs_type_name.compare(pcs_type_name)==0){ //OK/SB 4108
      // if ( convertProcessTypeToString(m_bc->getProcessType()).compare(_pcs_type_name) == 0 ) 
      //    && (m_bc->getProcessPrimaryVariable() == convertPrimaryVariable(_pcs_pv_name)))
      if ( m_bc->getProcess() == m_pcs )
      {
         //................................................................
         if (m_bc->getGeoType () == GEOLIB::POLYLINE)
         {
                                                  //CC
            m_polyline = GEOGetPLYByName(m_bc->geo_name);

            if (m_polyline)
            {
               //................................................................
                                                  //WW
               if (m_bc->getProcessDistributionType() == FiniteElement::LINEAR)
               {
                                                  //WW
                  if (m_polyline->getType() == 100)
                     m_msh->GetNodesOnArc(m_polyline, nodes_vector);
                  else
                  {
                                                  // *** REMOVE CANDIDATE
                     m_msh->GetNODOnPLY(m_polyline, nodes_vector);

                     double msh_min_edge_length = m_msh->getMinEdgeLength();
                     m_msh->setMinEdgeLength (m_polyline->epsilon);
                     std::vector<size_t> my_nodes_vector;

                     m_msh->GetNODOnPLY (static_cast<const GEOLIB::Polyline*>(m_bc->getGeoObj()), my_nodes_vector);
                     m_msh->setMinEdgeLength (msh_min_edge_length);

                     nodes_vector.clear ();
                     for (size_t k(0); k<my_nodes_vector.size(); k++)
                        nodes_vector.push_back (my_nodes_vector[k]);

                     // for some benchmarks we need the vector entries sorted by index
                     std::vector<size_t> perm;
                     for (size_t k(0); k<my_nodes_vector.size(); k++) perm.push_back (k);
                     Quicksort<long> (nodes_vector, 0, nodes_vector.size(), perm);
                  }

                  node_value.resize(nodes_vector.size());
                  //InterpolationAlongPolyline(m_polyline, node_value);
                  m_pcs->bc_transient_index.push_back(
                     (long) m_pcs->bc_node.size());
                  for (i = 0; i < (int) nodes_vector.size(); i++)
                  {
                     m_node_value = new CBoundaryConditionNode();
                     m_node_value->msh_node_number = -1;
                     m_node_value->msh_node_number = nodes_vector[i]
                        + ShiftInNodeVector;
                     m_node_value->geo_node_number = nodes_vector[i];
                     m_node_value->node_value = 0.0;
                                                  //YD/WW
                     m_node_value->pcs_pv_name = _pcs_pv_name;
                     m_node_value->CurveIndex = m_bc->getCurveIndex();
                                                  //WW
                     m_pcs->bc_node.push_back(m_bc);
                                                  //WW
                     m_pcs->bc_node_value.push_back(m_node_value);
                  }
                  node_value.clear();
               }
               //................................................................
               //delete(values);
               Free(nodes);
            }                                     // if(m_ply)
         }
         //------------------------------------------------------------------
      }                                           // PCS
      ++p_bc;
   }                                              // list
   /* // Make the following as comment by WW
    // Test
    long no_bc = (long)m_pcs->bc_node_value.size();
    if(no_bc<1)
    cout << "Warning: no boundary conditions specified for " << pcs_type_name << endl;
    */
}


/**************************************************************************
 FEMLib-Method: CBoundaryCondition::Get
 Task: set boundary conditions
 Programing:
 02/2004 OK Implementation
 last modification:
 **************************************************************************/
CBoundaryConditionsGroup* CBoundaryConditionsGroup::Get(const std::string& pcs_name)
{
   CBoundaryConditionsGroup *m_bc_group = NULL;
   std::list<CBoundaryConditionsGroup*>::const_iterator p_bc_group =
      bc_group_list.begin();
   while (p_bc_group != bc_group_list.end())
   {
      m_bc_group = *p_bc_group;
      if (m_bc_group->group_name.compare(pcs_name) == 0)
         return m_bc_group;
      ++p_bc_group;
   }
   return NULL;
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 04/2005 OK Implementation
 last modification:
 **************************************************************************/
CBoundaryConditionsGroup* BCGetGroup(const std::string& pcs_type_name, const std::string& pcs_pv_name)
{
   std::list<CBoundaryConditionsGroup*>::const_iterator it = bc_group_list.begin();
   while (it != bc_group_list.end())
   {
      if (((*it)->getProcessTypeName().compare(pcs_type_name) == 0)
         && ((*it)->getProcessPrimaryVariableName().compare(pcs_pv_name) == 0))
         return *it;
      ++it;
   }
   return NULL;
}


/**************************************************************************
 FEMLib-Method:
 Task: write function based on CBoundaryCondition::WriteTecplot
 Programing:
 01/2005 OK Implementation
 last modification:
 **************************************************************************/
void CBoundaryConditionsGroup::WriteTecplot() const
{
   // File handling
   std::string bc_path;
   CGSProject* m_gsp = GSPGetMember("bc");
   if (m_gsp)
      bc_path = m_gsp->path;
   std::string tec_file_name;
   tec_file_name = bc_path + "BC_" + group_name + ".tec";
   std::fstream tec_file(tec_file_name.data(), std::ios::trunc | std::ios::out);
   tec_file.setf(std::ios::scientific, std::ios::floatfield);
   tec_file.precision(12);
   // Write data
   /*
    // To be moved to pcs
    long i;
    long no_nodes = (long)group_vector.size();
    tec_file << "TITLE = Scatter Plot Types" << endl;
    tec_file << "VARIABLES = X,Y,Z,V" << endl;
    for(i=0;i<no_nodes;i++){
    tec_file << GetNodeX(group_vector[i]->msh_node_number) << ", ";
    tec_file << GetNodeY(group_vector[i]->msh_node_number) << ", ";
    tec_file << GetNodeZ(group_vector[i]->msh_node_number) << ", ";
    tec_file << group_vector[i]->node_value << endl;
   }
   */
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void BCDelete()
{
   CBoundaryCondition* m_bc = NULL;
   std::list<CBoundaryCondition*>::const_iterator p = bc_list.begin();
   while (p != bc_list.end())
   {
      //bc_list.remove(*p);
      m_bc = *p;
      delete m_bc;
      ++p;
   }
   bc_list.clear();
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void BCGroupDelete()
{
   CBoundaryConditionsGroup* m_bc_group = NULL;
   std::list<CBoundaryConditionsGroup*>::const_iterator p = bc_group_list.begin();
   while (p != bc_group_list.end())
   {
      m_bc_group = *p;
      delete m_bc_group;
      //bc_group_list.remove(*p);
      ++p;
   }
   bc_group_list.clear();
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 05/2005 OK Implementation
 last modification:
 **************************************************************************/
/*
 double CBoundaryConditionsGroup::GetConditionalNODValue(int i,CBoundaryCondition* m_bc)
 {
 int j;
 int nidx;
 int no_nodes = 4;
 double m_val = 0.0;
 CRFProcess* m_pcs_cond = NULL;
 CRFProcess* m_pcs_this = NULL;
 CNode* m_nod = NULL;
 Mesh_Group::CElem* m_ele = NULL;

//----------------------------------------------------------------------
m_pcs_this = PCSGet(m_bc->pcs_type_name);
m_pcs_cond = PCSGet(m_bc->pcs_type_name_cond);
m_nod = m_pcs_this->m_msh->nod_vector[group_vector[i]->msh_node_number];

//if(!m_nod->selected)
//  return 0.0;
//MB ??? CheckMarkedElement();

m_nod->connected_elements[0];
m_ele = m_pcs_cond->m_msh->ele_vector[0];

nidx = m_pcs_cond->GetNodeValueIndex(m_bc->pcs_pv_name_cond)+1;
for(j=0;j<no_nodes;j++){
m_val += m_pcs_cond->GetNodeValue(m_ele->GetNodesNumber(m_pcs_cond->m_msh->getOrder()),nidx);
}

m_val /= no_nodes;

//..................................................................
if(m_val > m_bc->condition){
//cout << m_val << " " << m_bc->node_value_cond << endl;
return m_bc->node_value_cond;
}
else{
//cout << m_val << " " << m_bc->geo_node_value << endl;
return m_bc->geo_node_value;
}
}
*/
/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void BCGroupDelete(const std::string& pcs_type_name, const std::string& pcs_pv_name)
{
   std::list<CBoundaryConditionsGroup*>::iterator p = bc_group_list.begin();
   while (p != bc_group_list.end())
   {
      if (((*p)->getProcessTypeName().compare(pcs_type_name) == 0)
         && ((*p)->getProcessPrimaryVariableName().compare(pcs_pv_name) == 0))
      {
         delete *p;
         bc_group_list.erase(p);
         return;
      }
      ++p;
   }
}


/**************************************************************************
 FEMLib-Method:
 07/2007 OK Implementation
 **************************************************************************/
CBoundaryCondition* BCGet(const std::string& pcs_type_name)
{
   CBoundaryCondition *m_bc = NULL;

   std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
   while (p_bc != bc_list.end())
   {
      m_bc = *p_bc;
      if (m_bc->getProcessType() == convertProcessType(pcs_type_name))
         return m_bc;
      ++p_bc;
   }
   return NULL;
}


/**************************************************************************
 ROCKFLOW - Funktion:
 Programming:
 11/2007 WW Implementation
 **************************************************************************/
void CBoundaryCondition::SurfaceInterpolation(CRFProcess* m_pcs, std::vector<long>&nodes_on_sfc, std::vector<double>&node_value_vector)
{

   long i, j, k, l;

   //----------------------------------------------------------------------
   // Interpolation of polygon values to nodes_on_sfc
   int nPointsPly = 0;
   double Area1, Area2;
   //NW. Default tolerance is 1e-9 but it can be changed in a BC file.
   double Tol = this->epsilon;
   bool Passed;
   double gC[3], p1[3], p2[3], pn[3], vn[3], unit[3], NTri[3];
   //
   CGLPolyline* m_polyline = NULL;
   Surface *m_surface = NULL;
   m_surface = GEOGetSFCByName(geo_name);         //CC

   // list<CGLPolyline*>::const_iterator p = m_surface->polyline_of_surface_list.begin();
   std::vector<CGLPolyline*>::iterator p =
      m_surface->polyline_of_surface_vector.begin();

   for (j = 0; j < (long)nodes_on_sfc.size(); j++)
   {
      pn[0] = m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->X();
      pn[1] = m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->Y();
      pn[2] = m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->Z();
      node_value_vector[j] = 0.0;
      Passed = false;
      // nodes close to first polyline
      p = m_surface->polyline_of_surface_vector.begin();
      while (p != m_surface->polyline_of_surface_vector.end())
      {
         m_polyline = *p;
         // Gravity center of this polygon
         for (i = 0; i < 3; i++)
            gC[i] = 0.0;
         vn[2] = 0.0;
         nPointsPly = (int) m_polyline->point_vector.size();
         for (i = 0; i < nPointsPly; i++)
         {
            gC[0] += m_polyline->point_vector[i]->x;
            gC[1] += m_polyline->point_vector[i]->y;
            gC[2] += m_polyline->point_vector[i]->z;
            vn[2] += m_polyline->point_vector[i]->getPropert();
         }
         for (i = 0; i < 3; i++)
            gC[i] /= (double) nPointsPly;
         // BC value at center is an average of all point values of polygon
         vn[2] /= (double) nPointsPly;
         // Area of this polygon by the gravity center
         for (i = 0; i < nPointsPly; i++)
         {
            p1[0] = m_polyline->point_vector[i]->x;
            p1[1] = m_polyline->point_vector[i]->y;
            p1[2] = m_polyline->point_vector[i]->z;
            k = i + 1;
            if (i == nPointsPly - 1)
               k = 0;
            p2[0] = m_polyline->point_vector[k]->x;
            p2[1] = m_polyline->point_vector[k]->y;
            p2[2] = m_polyline->point_vector[k]->z;
            vn[0] = m_polyline->point_vector[i]->getPropert();
            vn[1] = m_polyline->point_vector[k]->getPropert();

            Area1 = fabs(ComputeDetTri(p1, gC, p2));

            Area2 = 0.0;
            // Check if pn is in the triangle by points (p1, gC, p2)
            Area2 = fabs(ComputeDetTri(p2, gC, pn));
            unit[0] = fabs(ComputeDetTri(gC, p1, pn));
            unit[1] = fabs(ComputeDetTri(p1, p2, pn));
            Area2 += unit[0] + unit[1];
            if (fabs(Area1 - Area2) < Tol)
            {
               // Interpolation within a triangle (p1,p2,gC)
               // Shape function
               for (l = 0; l < 2; l++)
                  unit[l] /= Area1;
               ShapeFunctionTri(NTri, unit);
               for (l = 0; l < 3; l++)
                  node_value_vector[j] += vn[l] * NTri[l];
               Passed = true;
               break;
            }

         }
         //
         p++;
         if (Passed)
            break;
      }                                           // while
   }                                              //j
}
