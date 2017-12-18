/**************************************************************************
 FEMLib - Object: Source Terms ST
 Task:
 Programing:
 01/2004 OK Implementation
 last modified
 **************************************************************************/

#include "makros.h"
// C++ STL
//#include <fstream>
#include <cfloat>
#include <iostream>
#include <set>
//#include <algorithm>

//#include "files0.h"
// FEM-Makros
//#include "geo_lib.h"
#include "files0.h"
#include "mathlib.h"

// GeoSys-GeoLib
#include "GEOObjects.h"
// GEO
//#include "GeoType.h"

// GeoSys-MshLib
#include "fem_ele.h"
//#include "msh_elem.h"
//#include "msh_lib.h"

#include "tools.h"                                //GetLineFromFile
/* Tools */
#ifndef NEW_EQS                                   //WW. 06.11.2008
#include "matrix.h"
#endif

// GeoSys-FEMLib
//OK_IC #include "rfsousin.h"
#include "rf_st_new.h"
#include "rf_tim_new.h"

// Math
#include "matrix_class.h"

//#include "pcs_dm.h"

// FEM
//#include "problem.h"
// For analytical source terms
#include "rf_mfp_new.h"
#include "rf_node.h"
#include "rfmat_cp.h"

// Base
//#include "StringTools.h"

// MathLib
#include "LinearInterpolation.h"

// FileIO
#include "FEMIO/GeoIO.h"

using FiniteElement::CElement;
using Mesh_Group::CElem;
using Mesh_Group::CEdge;
using Mesh_Group::CNode;
using Math_Group::vec;

#ifndef GRAVITY_CONSTANT
#define GRAVITY_CONSTANT 9.81
#endif

//==========================================================================
std::vector<CSourceTerm*> st_vector;
std::list<CSourceTermGroup*> st_group_list;
std::vector<std::string> analytical_processes;
std::vector<std::string> analytical_processes_polylines;
std::vector<NODE_HISTORY*> node_history_vector;   //CMCD
/**************************************************************************
 FEMLib-Method:
 Task: ST constructor
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
CSourceTerm::CSourceTerm() :
	ProcessInfo(), GeoInfo(), _coupled (false), _sub_dom_idx(-1),
	start_pos_in_st(-1), GIS_shape_head(NULL)
                                                  // 07.06.2010, 03.2010. WW
{
   CurveIndex = -1;
   //KR critical_depth = false;
   //	COUPLING_SWITCH = false;
   geo_node_value = 0.0;
   nodes = NULL;                                  //OK
   analytical = false;                            //CMCD
   //  display_mode = false; //OK
}


/**************************************************************************
 FEMLib-Method:
 Task: BC deconstructor
 Programing:
 04/2004 OK Implementation
 **************************************************************************/
CSourceTerm::~CSourceTerm()
{
   DeleteHistoryNodeMemory();
   //    dis_file_name.clear();
   node_number_vector.clear();
   node_value_vector.clear();
   node_renumber_vector.clear();
   PointsHaveDistribedBC.clear();
   DistribedBC.clear();
   element_st_vector.clear();
   //WW----------22.02.2007-------------------
   // TF 06/2010
   size_t size(normal2surface.size());
   for (size_t i = 0; i < size; i++)
      delete normal2surface[i];
   size = pnt_parameter_vector.size();
   for (size_t i = 0; i < size; i++)
      delete pnt_parameter_vector[i];
   if(GIS_shape_head)                             // 07.06.2010. WW
   {
      delete [] GIS_shape_head;
      GIS_shape_head = NULL;
   }
   //WW---------------------------------------
}


const std::string& CSourceTerm::getGeoName() const
{
   return geo_name;
}


double CSourceTerm::getCoupLeakance () const
{
   return _coup_leakance;
}


/**************************************************************************
 FEMLib-Method:
 Task: ST read function
 Programing:
 01/2004 OK Implementation
 11/2004 MB neues read Konzept
 02/2005 MB River condition
 03/2005 WW Node force released by excavation
 11/2005 CMCD Analytical source for matrix
 04/2006 OK CPL
 04/2006 OK MSH_TYPE
06/2010 TF modification of the signature, added geo_obj and unique_name
**************************************************************************/
std::ios::pos_type CSourceTerm::Read(std::ifstream *st_file,
const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
   char line[MAX_ZEILE];
   std::string line_string, sub_string;
   bool new_keyword = false;

   std::stringstream in;

                                                  // JOD 4.10.01
   channel = 0, node_averaging = 0, no_surface_water_pressure = 0;
   std::ios::pos_type position;

   // read loop
   while (!new_keyword)
   {
      position = st_file->tellg();
      if (!GetLineFromFile(line, st_file))
         break;
      line_string = line;
      if (line_string.find("#") != std::string::npos)
      {
         new_keyword = true;
         break;
      }
      remove_white_space(&line_string);           //OK

      /* search for keywords */
                                                  // subkeyword found
      if (line_string.find("$PCS_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(st_file));
         std::string tmp;
         in >> tmp;
         setProcessType (convertProcessType (tmp));
         in.clear();
         continue;
      }

                                                  // subkeyword found
      if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos)
      {
         in.str(GetLineFromFile1(st_file));
         std::string tmp;
         in >> tmp;
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
                 DisplayErrorMsg("Error: In reading ST file, the input component names are not found in MCP file!!!");
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

      if (line_string.find("$COMP_NAME") != std::string::npos)
      {
         in.str(GetLineFromFile1(st_file));
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
             DisplayErrorMsg("Error: In reading ST file, the input component names are not found in MCP file!!!");
             exit(1);
         }
         in.clear();
         continue;
      }

      if (line_string.find("$GEO_TYPE") != std::string::npos)
      {
         ReadGeoType(st_file, geo_obj, unique_name);
         continue;
      }

                                                  //05.09.2008 WW
      if (line_string.find("$DIS_TYPE") != std::string::npos)
      {
         //10.04.2008. WW  if(line_string.compare("$DIS_TYPE")==0) {
         if (line_string.find("CONDITION") != std::string::npos)
         {
            _coupled = true;
            ReadDistributionType(st_file);
            in.str(GetLineFromFile1(st_file));    //PCS type
            in >> line_string >> pcs_type_name_cond;
            in.clear();
            in.str(GetLineFromFile1(st_file));    //
            in >> pcs_pv_name_cond;
            in.clear();
            in.str(GetLineFromFile1(st_file));
            in >> _coup_leakance >> rill_height;
            in.clear();
         }                                        //05.09.2008 WW
         else
         {
            ReadDistributionType(st_file);
            continue;
         }
      }

      //		if (line_string.find("$COUPLING_SWITCH") != std::string::npos) {
      //			COUPLING_SWITCH = true; // switch case
      //			in.str(GetLineFromFile1(st_file));
      //			in >> rainfall >> rainfall_duration;
      //			in.clear();
      //		}

      if (line_string.find("$NODE_AVERAGING") != std::string::npos)
      {
         in.clear();
         node_averaging = 1;
         continue;
      }
                                                  // JOD 4.10.01
      if (line_string.find("$NEGLECT_SURFACE_WATER_PRESSURE") != std::string::npos)
      {
         in.clear();
         no_surface_water_pressure = 1;
         continue;
      }
      if (line_string.find("$CHANNEL") != std::string::npos)
      {
         in.clear();
         in.str(GetLineFromFile1(st_file));
         in >> channel_width;
         channel = 1;
         continue;
      }

      if (line_string.find("$TIM_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(st_file));
         in >> tim_type_name;
         if (tim_type_name.find("CURVE") != std::string::npos)
         {
        	 //				dis_type = 0;
            in >> CurveIndex;
         }
         in.clear();
         continue;
      }

      if (line_string.find("$FCT_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(st_file));
         in >> fct_name;                          //sub_line
                                                  //WW
         if (fct_name.find("METHOD") != std::string::npos)
            in >> fct_method;
         in.clear();
      }

      if (line_string.find("$MSH_TYPE") != std::string::npos)
      {
         in.str(GetLineFromFile1(st_file));
         std::string sub_string;
         in >> sub_string;                        //sub_line
         msh_type_name = "NODE";
         if (sub_string.find("NODE") != std::string::npos)
         {
            in >> msh_node_number;
            in.clear();
         }
         continue;
      }
   }                                              // end !new_keyword
   return position;
}


/**************************************************************************
 FEMLib-Method:
 Task: for CSourceTerm::Read
 Programing:
 11/2007 JOD Implementation
 02/2009 WW  Add a functionality to directly assign source terms to element nodes.
 **************************************************************************/
void CSourceTerm::ReadDistributionType(std::ifstream *st_file)
{
   std::stringstream in;
   // 03.2010 WW
   std::string aline;
   std::stringstream ss;
   int abuff, nLBC = 0;
   double bbuff;

   std::string dis_type_name;
   in.str(GetLineFromFile1(st_file));
   in >> dis_type_name;

   this->setProcessDistributionType (FiniteElement::convertDisType(dis_type_name));

   if (dis_type_name.compare (convertDisTypeToString (this->getProcessDistributionType())) != 0)
   {
      std::cerr << "Error in CSourceTerm::ReadDistributionType (): dist_type_name #" << dis_type_name << "#, new: " << convertDisTypeToString (this->getProcessDistributionType()) << std::endl;
      exit (1);
   }

   if (   this->getProcessDistributionType() == FiniteElement::CONSTANT
       || this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN
       || this->getProcessDistributionType() == FiniteElement::CONSTANT_GEO      )
   {
      in >> geo_node_value;
      in.clear();
   }

   //	if (dis_type_name.find("ANALYTICAL") != std::string::npos) {
   if (this->getProcessDistributionType() == FiniteElement::ANALYTICAL)
   {
      in >> analytical_material_group;            //Which material group is it being applied to
      in >> analytical_diffusion;                 //D value
      in >> analytical_porosity;                  //n value of matrix
      in >> analytical_tortousity;                //t value of matrix
      in >> analytical_linear_sorption_Kd;        //Linear sorption coefficient
      in >> analytical_matrix_density;            //Density of solid
      in >> number_of_terms;                      //no timesteps to consider in solution
      in >> resolution;                           //every nth term will be considered
      in >> factor;                               //to convert temperature to energy
      analytical = true;
      analytical_processes.push_back(convertPrimaryVariableToString (getProcessPrimaryVariable()));
      //		if (geo_type_name.compare("POLYLINE") == 0)
      if (this->getGeoType() == GEOLIB::POLYLINE)
         analytical_processes_polylines.push_back(geo_name);
      in.clear();
   }

   if (this->getProcessDistributionType() == FiniteElement::LINEAR || this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
   {
      in >> nLBC;
      in.clear();
      for (int i = 0; i < nLBC; i++)
      {
         in.str(GetLineFromFile1(st_file));
         in >> abuff >> bbuff;
         in.clear();
         PointsHaveDistribedBC.push_back(abuff);
         DistribedBC.push_back(bbuff);
      }

      //      Read LINENODES AND VALUES......
      in.clear();
   }

   if (this->getProcessDistributionType() == FiniteElement::CRITICALDEPTH)
   {
      //KR critical_depth = true;
      in >> geo_node_value;
      in.clear();
      in.str(GetLineFromFile1(st_file));
      in >> rill_height;
      in.clear();
      //		dis_type = 6;
   }

   if (this->getProcessDistributionType() == FiniteElement::NORMALDEPTH)
   {
      dis_type_name = "NORMALDEPTH";
      in >> geo_node_value;
      in.clear();
      in.str(GetLineFromFile1(st_file));
      in >> normaldepth_slope >> rill_height;
      in.clear();
   }


   if (this->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
   {
      dis_type_name = "GREEN_AMPT";
      in >> geo_node_value;
      in.clear();
      in.str(GetLineFromFile1(st_file));
      in >> sorptivity >> constant >> rainfall >> moistureDeficit;
      in.clear();
   }
   // Soure terms are assign to element nodes directly. 23.02.2009. WW
   if(getProcessDistributionType() == FiniteElement::DIRECT)
   {
      dis_type_name = "DIRECT";
      in >> fname;
      fname = FilePath+fname;
      in.clear();
   }
   // Soure terms from precipitation are assign to element nodes directly.03.2010. WW
   if(dis_type_name.find("PRECIPITATION")!=std::string::npos)
   {
      dis_type_name = "PRECIPITATION";
      in >> fname;
      fname = FilePath+fname;
      std::ifstream ins(fname.c_str());
      if(!ins.good())
      {
         std::cout<<"Could not find file "<<fname<<std::endl;
         exit(0);
      }
      double timess;
      GIS_shape_head = new double[6];             // 07.06.2010. WW
      for(int i=0; i<6; i++)
      {
         getline(ins, aline);
         ss.str(aline);
         ss>> aline >> GIS_shape_head[i];
         ss.clear();

      }
      while(!ins.eof())
      {
         getline(ins, aline);
         if(aline.find("#STOP")!=std::string::npos)
            break;
         ss.str(aline);
         ss>> timess >> aline;
         precip_times.push_back(timess);
         precip_files.push_back(aline);

         ss.clear();
      }
      in.clear();
   }
}


/**************************************************************************
 FEMLib-Method:
 Task: for CSourceTerm::Read
 Programing:
 11/2007 JOD Implementation
 06/2010 TF modification of the signature, added geo_obj and unique_name
 **************************************************************************/
void CSourceTerm::ReadGeoType(std::ifstream *st_file,
const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
   FileIO::GeoIO geo_io;
   geo_io.readGeoInfo(this, *st_file, geo_name, geo_obj, unique_name);

   if (getProcessPrimaryVariable() == EXCAVATION) //WW
   {
      std::stringstream strstr;
      strstr.str(GetLineFromFile1(st_file));
      //size_t tmp_geo_type;
      std::string sub_string;
      strstr >> sub_string >> _sub_dom_idx;
      strstr.clear();
   }
}


/**************************************************************************
 FEMLib-Method:
 Task: ST read function
 Programing:
 01/2004 OK Implementation
 06/2010 TF modification of the signature, added geo_obj and unique_name
 **************************************************************************/
bool STRead(const std::string &file_base_name,
const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
   char line[MAX_ZEILE];
   std::string line_string, st_file_name;
   std::ios::pos_type position;

   // File handling
   st_file_name = file_base_name + ST_FILE_EXTENSION;
   std::ifstream st_file(st_file_name.data(), std::ios::in);

   if (!st_file.good())
   {
      std::cout << "! Warning in STRead: No source terms !" << std::endl;
      return false;
   }

   // Keyword loop
   std::cout << "STRead" << std::endl;
   while (!st_file.eof())
   {
      st_file.getline(line, MAX_ZEILE);
      line_string = line;
                                                  //Code included to make dynamic memory for analytical solution
      if (line_string.find("#STOP") != std::string::npos)
      {
         size_t no_source_terms(st_vector.size());
         size_t no_an_sol = 0, number_of_terms = 0;
                                                  //Need to find the memory size limits for the anal. solution.
         for (size_t i = 0; i < no_source_terms; i++)
         {
            if (st_vector[i]->isAnalytical())
            {
               no_an_sol++;
               number_of_terms = std::max(st_vector[i]->getNumberOfTerms(), number_of_terms);
            }
         }
         if (no_an_sol > 0)
         {
            for (size_t i = 0; i < no_source_terms; i++)
            {
               st_vector[i]->setNumberOfAnalyticalSolutions (no_an_sol);
               st_vector[i]->setMaxNumberOfTerms (number_of_terms);
            }
         }

         return true;
      }
      //----------------------------------------------------------------------
                                                  // keyword found
      if (line_string.find("#SOURCE_TERM") != std::string::npos)
      {
         CSourceTerm *st(new CSourceTerm());
         std::ios::pos_type pos (st_file.tellg());
         position = st->Read(&st_file, geo_obj, unique_name);
         if (pos != position)
         {
            st_vector.push_back(st);
         }
         else
         {
            std::cerr << "WARNING: in STRead: could not read source term" << std::endl;
            delete st;
         }
         st_file.seekg(position, std::ios::beg);
      }                                           // keyword found
   }                                              // eof
   return true;
}


/**************************************************************************
 FEMLib-Method:
 Task: ST to right-hand-side vector
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
void ST2RHS(std::string pcs_function, double* rhs_vector)
{
   int j;
   long i;
   long node_number_vector_length;
   CSourceTerm *m_st = NULL;
   long no_st = (long) st_vector.size();
   for (j = 0; j < no_st; j++)
   {
      m_st = st_vector[j];
      if (pcs_function.compare(convertPrimaryVariableToString (m_st->getProcessPrimaryVariable())) == 0)
      {
         node_number_vector_length = (long) m_st->node_number_vector.size();
         for (i = 0; i < node_number_vector_length; i++)
         {
            rhs_vector[m_st->node_number_vector[i]]
               = m_st->node_value_vector[i];
         }
      }
   }
}


/**************************************************************************
 FEMLib-Method:
 Task: ST to mesh nodes
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
//void ST2NOD()
//{
//	CGLVolume *m_volume = NULL;
//	double st_value;
//
//	CGLPolyline *ply (NULL);
//	size_t points_along_polyline;
//
//	// Nodes
//	size_t no_st (st_vector.size());
//	for (size_t j = 0; j < no_st; j++) {
//		CSourceTerm *m_st = st_vector[j];
//		switch (m_st->getGeoType()) {
//		case GS_POINT:
//			m_st->node_number_vector.push_back(m_st->getGeoObjIdx());
//			break;
//		case GS_POLYLINE:
//			CGLPolyline *ply = GEOGetPLYByName(m_st->geo_prop_name);
//			points_along_polyline = ply->point_vector.size();
//			for (size_t i = 0; i < points_along_polyline; i++) {
//				m_st->node_number_vector.push_back(ply->point_vector[i]->id);
//			}
//			break;
//		case GS_SURFACE:
//			Surface *m_surface (GEOGetSFCByName(m_st->geo_prop_name));
//			long points_in_surface;
//			long* nodes (NULL);
//			nodes = GetPointsIn(m_surface, &points_in_surface);
//			//MB patch areas
//			for (size_t i = 0; i < static_cast<size_t>(points_in_surface); i++) {
//				m_st->node_number_vector.push_back(nodes[i]);
//			}
//			delete [] nodes;
//			break;
//		case GS_VOLUME:
//			m_volume = GEOGetVOL(m_st->geo_prop_name);
//			break;
//		default:
//			break;
//		} // switch
//	} // while
//	//========================================================================
//	size_t st_point_number;
//	size_t geo_point_number;
//	// Values
//	for (size_t j = 0; j < no_st; j++) {
//		CSourceTerm *m_st = st_vector[j];
//		switch (m_st->dis_type) {
//		case CONSTANT:
//			st_point_number = m_st->node_number_vector.size();
//			for (size_t i = 0; i < st_point_number; i++) {
//				m_st->node_value_vector.push_back(m_st->dis_prop[0]);
//			}
//			break;
//		case LINEAR: // for polylines
//			CGLPolyline *ply = GEOGetPLYByName(m_st->geo_prop_name);//CC
//
//			geo_point_number = ply->point_vector.size();
//			//if(!(st_point_number==geo_point_number)) Warning !
//			for (size_t i = 0; i < geo_point_number; i++) {
//				st_value = ply->point_vector[i]->propert; // i.e. node property
//				m_st->node_value_vector.push_back(st_value);
//			}
//			break;
//		} // switch
//	} // while
//}

/**************************************************************************
 FEMLib-Method: STWrite
 Task: master write function
 Programing:
 04/2004 OK Implementation
 last modification:
 05/2010 TF
 **************************************************************************/
void STWrite(std::string base_file_name)
{
   // File handling
   std::string st_file_name = base_file_name + ST_FILE_EXTENSION;
   std::fstream st_file(st_file_name.data(), std::ios::trunc | std::ios::out);
   st_file.setf(std::ios::scientific, std::ios::floatfield);
   st_file.precision(12);
   if (!st_file.good())
      return;

   st_file
      << "GeoSys-ST: Source Terms ------------------------------------------------\n";

   size_t no_st(st_vector.size());
   for (size_t j = 0; j < no_st; j++)
   {
      st_vector[j]->Write(&st_file);
   }
   st_file << "#STOP";
   st_file.close();
}


/**************************************************************************
 FEMLib-Method:
 Task: write function
 Programing:
 02/2004 OK Implementation
 04/2005 OK PRIMARY_VARIABLE
 06/2005 OK RIVER
 last modification:
 **************************************************************************/
void CSourceTerm::Write(std::fstream* st_file)
{
   //KEYWORD
   *st_file << "#SOURCE_TERM" << std::endl;
   //--------------------------------------------------------------------
   //NAME+NUMBER
   *st_file << " $PCS_TYPE" << std::endl;
   *st_file << "  ";
   //	*st_file << pcs_type_name << std::endl;
   *st_file << convertProcessTypeToString (getProcessType()) << std::endl;
   *st_file << " $PRIMARY_VARIABLE" << std::endl;
   *st_file << "  ";
   *st_file << convertPrimaryVariableToString (getProcessPrimaryVariable()) << std::endl;
   //--------------------------------------------------------------------
   //GEO_TYPE
   *st_file << " $GEO_TYPE" << std::endl;
   *st_file << "  ";
   *st_file << getGeoTypeAsString() << " " << geo_name << std::endl;
   //--------------------------------------------------------------------
   // TIM_TYPE
   if (tim_type_name.size() > 0)                  //OK
   {
      *st_file << " $TIM_TYPE" << std::endl;
      *st_file << "  ";
      *st_file << tim_type_name << std::endl;
   }
   //--------------------------------------------------------------------
   //DIS_TYPE
   *st_file << " $DIS_TYPE" << std::endl;
   *st_file << "  ";
   *st_file << convertDisTypeToString(this->getProcessDistributionType());
   switch (this->getProcessDistributionType())
   {
      case FiniteElement::CONSTANT:
         *st_file << " " << geo_node_value;
         *st_file << std::endl;
         break;
      case  FiniteElement::CONSTANT_NEUMANN:
         *st_file << " " << geo_node_value;
         *st_file << std::endl;
         break;
      case FiniteElement::LINEAR:
         *st_file << " " << (int) PointsHaveDistribedBC.size() << std::endl;
         for (long i = 0; i < (long) PointsHaveDistribedBC.size(); i++)
         {
            *st_file << "  " << PointsHaveDistribedBC[i] << " ";
            *st_file << "  " << DistribedBC[i] << std::endl;
         }
         break;
      case  FiniteElement::LINEAR_NEUMANN:
         *st_file << " " << PointsHaveDistribedBC.size() << std::endl;
         for (size_t i = 0; i < PointsHaveDistribedBC.size(); i++)
         {
            *st_file << "  " << PointsHaveDistribedBC[i] << " ";
            *st_file << "  " << DistribedBC[i] << std::endl;
         }
         break;
      default:
         std::cerr << "this distributition type is not handled in CSourceTerm::Write" << std::endl;
   }
   //--------------------------------------------------------------------
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 04/2004 OK Implementation
 11/2007 JOD Reaktivation
 last modification:
 **************************************************************************/
//void CSourceTerm::SetDISType()
//{
////	if (this->getProcessDistributionType() == CONSTANT)
////		dis_type = 1;
////////	if (dis_type_name.compare("CONSTANT_GEO") == 0)
////////		dis_type = 12; //SB flux is distributed along polyline. To do. 4.10.06
////	if (this->getProcessDistributionType() == LINEAR)
////		dis_type = 2;
////	if (this->getProcessDistributionType() == CONSTANT_NEUMANN)
////		dis_type = 3;
////	if (this->getProcessDistributionType() == LINEAR_NEUMANN)
////		dis_type = 4;
////	if (this->getProcessDistributionType() == RIVER) {
////		dis_type = 5;
////	}
//
////	if (this->getProcessDistributionType() == CRITICALDEPTH)
////		dis_type = 6;
////	if (this->getProcessDistributionType() == SYSTEM_DEPENDENT)
////		dis_type = 7; //YD
////	if (this->getProcessDistributionType() == NORMALDEPTH)
////		dis_type = 8; //JOD MB
////	if (this->getProcessDistributionType() == ANALYTICAL)
////		dis_type = 9;//CMCD 02 2006
//////	if (dis_type_name.compare("PHILIP") == 0)
////////	if (this->getProcessDistributionType() == PHILIP)
//////		dis_type = 10; // JOD
////	if (this->getProcessDistributionType() == GREEN_AMPT)
////		dis_type = 11; // JOD
//	// if(dis_type_name.compare("CONSTANT")==0) dis_type = 0;
//	// if(dis_type_name.compare("LINEAR")  ==0) dis_type = 1;
//}

/**************************************************************************
 FEMLib-Method:
 Task: set ST group member
 Programing:
 02/2004 OK Implementation
 09/2004 WW Face integration of Neumann boundary condition for all element type
 09/2004 WW Interpolation for piece-wise linear distributed source term or BC
 03/2005 OK LINE sources
 02/2005 MB River condition, CriticalDepth
 08/2006 WW Re-implementing edge,face and domain integration versatile for all element types
 04/2006 OK MSH types
02/2009 WW Direct assign node source terms
**************************************************************************/
void CSourceTermGroup::Set(CRFProcess* m_pcs, const int ShiftInNodeVector,
		std::string this_pv_name)
{

   if (this_pv_name.size() != 0)                  //WW
      pcs_pv_name = this_pv_name;
   m_msh = FEMGet(pcs_type_name);                 //m_pcs->m_msh; //

   if (m_msh)                                     //WW
   {
      /// In case of P_U coupling monolithic scheme
      if (m_pcs->type == 41)                      //WW Mono
      {
                                                  //Deform
         if (pcs_pv_name.find("DISPLACEMENT") != std::string::npos)
            m_pcs->m_msh->SwitchOnQuadraticNodes(true);
         else
            m_pcs->m_msh->SwitchOnQuadraticNodes(false);
      } else if (m_pcs->type == 4)
      m_pcs->m_msh->SwitchOnQuadraticNodes(true);
      else
         m_pcs->m_msh->SwitchOnQuadraticNodes(false);
      //====================================================================
      long no_st = (long) st_vector.size();
      for (long i = 0; i < no_st; i++)
      {
         CSourceTerm *source_term (st_vector[i]);

         // 07.01.2011. WW
         if(source_term->getProcessDistributionType()==FiniteElement::PRECIPITATION)
            continue;

         //-- 07.04.2011--23.02.3009. WW
         if(source_term->getProcessDistributionType()==FiniteElement::DIRECT)
         {
            source_term->DirectAssign(ShiftInNodeVector);
            continue;
         }



         if (source_term->isCoupled())
            m_msh_cond = FEMGet(source_term->pcs_type_name_cond);

         if (source_term->getProcessType() == MASS_TRANSPORT)
             if ( cp_vec[cp_name_2_idx[convertPrimaryVariableToString(source_term->getProcessPrimaryVariable())]]->getProcess() != m_pcs )
                 continue;

         if ((convertProcessTypeToString (source_term->getProcessType ()).compare(pcs_type_name) == 0)
            && (convertPrimaryVariableToString(source_term->getProcessPrimaryVariable()).compare(pcs_pv_name) == 0))
         // if ( source_term->getProcess() == m_pcs )
         {
             source_term->setProcess (m_pcs);      // HS: 01.09.2009
             if (source_term->getGeoType() == GEOLIB::POINT)
                 SetPNT(m_pcs, source_term, ShiftInNodeVector);
             if (source_term->getGeoType () == GEOLIB::POLYLINE)
                 SetPLY(source_term, ShiftInNodeVector);
             if (source_term->getGeoType () == GEOLIB::SURFACE)
                 SetSFC(source_term, ShiftInNodeVector);
             if (source_term->getGeoType () == GEOLIB::GEODOMAIN)
                 SetDMN(source_term, ShiftInNodeVector);
             if (source_term->fct_name.size() > 0)
                 fct_name = source_term->fct_name;
         }                                        // end pcs name & pv
      }                                           // end st loop
   }                                              // end msh
   else
      std::cout << "Warning in CSourceTermGroup::Set - no MSH data" << std::endl;

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 08/2004 OK Implementation
 last modification:
 **************************************************************************/
//void STCreateFromPNT() {
//	size_t no_pnt_properties(pnt_properties_vector.size());
//
//	CSourceTerm *m_st = NULL;
//	CGLPoint *m_pnt = NULL;
//
//	for (size_t i = 0; i < no_pnt_properties; i++) {
//		m_pnt = pnt_properties_vector[i];
//		m_st = new CSourceTerm();
//		// PCS ---------------------------------------------------------------
//		m_st->pcs_pv_name = "PRESSURE1"; // ToDo
//		// GEO ---------------------------------------------------------------
//		m_st->geo_type_name = "POINT";
//		m_st->setGeoType (GEOLIB::POINT);
//		m_st->geo_node_number = m_pnt->id;
//		m_st->setGeoObjIdx(m_st->geo_node_number); // 05/2010 TF
//
//		// DIS ---------------------------------------------------------------
//		m_st->dis_type_name = "CONSTANT";
//		m_st->dis_type = 0;
//		m_st->geo_node_value = m_pnt->value;
//		// TIM ---------------------------------------------------------------
//		m_st->tim_type_name = m_pnt->name;
//		// VEC ---------------------------------------------------------------
//		st_vector.push_back(m_st);
//	}
//	pnt_properties_vector.clear();
//}

/**************************************************************************
 ROCKFLOW - Funktion: FaceIntegration
 Task: Translate distributed Neumann boundary condition /source term on edges
 found on a polyline to nodes value for all kinds of element
 Programming:
 07/2005 WW Re-Implementation
 12/2005 WW Axismmetry
 **************************************************************************/
//void CSourceTerm::EdgeIntegration(CFEMesh* msh, vector<long>&nodes_on_ply,
//		vector<double>&node_value_vector) {
//	long i, j, k, l;
//	long this_number_of_nodes;
//	int elemsCnode;
//	int nedges, ii;
//	vec<CNode*> e_nodes(3);
//	vec<CEdge*> e_edges(12);
//
//	double Jac = 0.0;
//	double Weight = 0.0;
//	double eta = 0.0;
//	double v1, v2, radius = 0.0;
//	double Shfct[3];
//	bool Const = false;
//
//	if (dis_type_name.find("CONSTANT") != std::string::npos)
//		Const = true;
//
//	//CFEMesh* msh = m_pcs->m_msh;
//	//CFEMesh* msh;  // JOD
//	//msh = FEMGet(pcs_type_name);
//	CElem* elem = NULL;
//	CEdge* edge = NULL;
//	CNode* node = NULL;
//
//	int nSize = (long) msh->nod_vector.size();
//	this_number_of_nodes = (long) nodes_on_ply.size();
//	vector<long> G2L(nSize);
//	vector<double> NVal(this_number_of_nodes);
//
//	// Unmakr edges.
//	for (i = 0; i < (long) msh->edge_vector.size(); i++)
//		msh->edge_vector[i]->SetMark(false);
//	for (i = 0; i < nSize; i++) {
//		msh->nod_vector[i]->SetMark(false);
//		G2L[i] = -1;
//	}
//
//	// Search edges on polyline
//	for (i = 0; i < this_number_of_nodes; i++) {
//		NVal[i] = 0.0;
//		k = nodes_on_ply[i];
//		G2L[k] = i;
//		node = msh->nod_vector[k];
//		elemsCnode = (int) node->connected_elements.size();
//		for (j = 0; j < elemsCnode; j++) {
//			l = msh->nod_vector[k]->connected_elements[j];
//			elem = msh->ele_vector[l];
//			nedges = elem->GetEdgesNumber();
//			elem->GetEdges(e_edges);
//			for (ii = 0; ii < nedges; ii++) {
//				edge = e_edges[ii];
//				if (edge->GetMark())
//					continue;
//				edge->GetNodes(e_nodes);
//				// Edge A
//				if (*node == *e_nodes[0])
//					e_nodes[0]->SetMark(true);
//				// Edge B
//				if (*node == *e_nodes[1])
//					e_nodes[1]->SetMark(true);
//				if (msh->getOrder()) // Quadratic
//				{
//					if (*node == *e_nodes[2])
//						e_nodes[2]->SetMark(true);
//				}
//				if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark()) {
//					if (msh->getOrder()) {
//						if (e_nodes[2]->GetMark())
//							edge->SetMark(true);
//					} else
//						edge->SetMark(true);
//				}
//
//			}// e_edges
//		}
//	}
//
//	for (i = 0; i < (long) msh->edge_vector.size(); i++) {
//		edge = msh->edge_vector[i];
//		if (!edge->GetMark())
//			continue;
//		edge->GetNodes(e_nodes);
//		if (msh->getOrder()) // Quad
//		{
//			if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark()
//					&& e_nodes[2]->GetMark()) {
//				Jac = 0.5 * edge->Length();
//				v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
//				v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
//				if (Const && (!msh->isAxisymmetry())) {
//					NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1 / 3.0;
//					NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1 / 3.0;
//					NVal[G2L[e_nodes[2]->GetIndex()]] += 4.0 * Jac * v1 / 3.0;
//
//				} else {
//					for (k = 0; k < 3; k++) // Three nodes
//					{
//						// Numerical integration
//						for (l = 0; l < 3; l++) // Gauss points
//						{
//							Weight = Jac * MXPGaussFkt(3, l);
//							eta = MXPGaussPkt(3, l);
//							ShapeFunctionLineHQ(Shfct, &eta);
//							//Axisymmetical problem
//							if (msh->isAxisymmetry()) {
//								radius = 0.0;
//								for (ii = 0; ii < 3; ii++)
//									radius += Shfct[ii] * e_nodes[ii]->X();
//								Weight *= radius; //2.0*pai*radius;
//							}
//							NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
//									+ eta * (v2 - v1)) * Shfct[k] * Weight;
//						}
//					}
//				}
//			}
//		} else // Linear
//		{
//			if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark()) {
//				Jac = 0.5 * edge->Length();
//				v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
//				v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
//				if (!msh->isAxisymmetry()) {
//					if (Const) {
//						NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1;
//						NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1;
//					} else {
//						NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * (2.0 * v1
//								+ v2) / 3.0;
//						NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * (v1 + 2.0
//								* v2) / 3.0;
//					}
//				} else // Axisymmetry
//				{
//
//					for (k = 0; k < 2; k++) // Three nodes
//					{
//						// Numerical integration
//						for (l = 0; l < 3; l++) // Gauss points
//						{
//							Weight = Jac * MXPGaussFkt(3, l);
//							eta = MXPGaussPkt(3, l);
//							ShapeFunctionLine(Shfct, &eta);
//							//Axisymmetical problem
//							if (msh->isAxisymmetry()) {
//								radius = 0.0;
//								for (ii = 0; ii < 2; ii++)
//									radius += Shfct[ii] * e_nodes[ii]->X();
//								Weight *= radius; //2.0*pai*radius;
//							}
//							NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
//									+ eta * (v2 - v1)) * Shfct[k] * Weight;
//						}
//					}
//				}// End of is (!axi)
//			}
//		}
//	}
//	for (i = 0; i < this_number_of_nodes; i++)
//		node_value_vector[i] = NVal[i];
//	for (i = 0; i < (long) msh->edge_vector.size(); i++)
//		msh->edge_vector[i]->SetMark(true);
//	for (i = 0; i < nSize; i++)
//		msh->nod_vector[i]->SetMark(true);
//	NVal.clear();
//	G2L.clear();
//	e_nodes.resize(0);
//	e_edges.resize(0);
//}

void CSourceTerm::EdgeIntegration(CFEMesh* msh, const std::vector<long>&nodes_on_ply,
std::vector<double>&node_value_vector) const
{
   long i, j, k, l;
   long this_number_of_nodes;
   int elemsCnode;
   int nedges, ii;
   vec<CNode*> e_nodes(3);
   vec<CEdge*> e_edges(12);

   double Jac = 0.0;
   double Weight = 0.0;
   double eta = 0.0;
   double v1, v2, radius = 0.0;
   double Shfct[3];

   bool Const = false;
   if (this->getProcessDistributionType() == FiniteElement::CONSTANT || this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
      Const = true;

   //CFEMesh* msh = m_pcs->m_msh;
   //CFEMesh* msh;  // JOD
   //msh = FEMGet(pcs_type_name);
   CElem* elem = NULL;
   CEdge* edge = NULL;
   CNode* node = NULL;

   int nSize = (long) msh->nod_vector.size();
   this_number_of_nodes = (long) nodes_on_ply.size();
   std::vector<long> G2L(nSize);
   std::vector<double> NVal(this_number_of_nodes);

   // Unmakr edges.
   for (i = 0; i < (long) msh->edge_vector.size(); i++)
      msh->edge_vector[i]->SetMark(false);
   for (i = 0; i < nSize; i++)
   {
      msh->nod_vector[i]->SetMark(false);
      G2L[i] = -1;
   }

   // Search edges on polyline
   for (i = 0; i < this_number_of_nodes; i++)
   {
      NVal[i] = 0.0;
      k = nodes_on_ply[i];
      G2L[k] = i;
      node = msh->nod_vector[k];
      elemsCnode = (int) node->connected_elements.size();
      for (j = 0; j < elemsCnode; j++)
      {
         l = msh->nod_vector[k]->connected_elements[j];
         elem = msh->ele_vector[l];
         nedges = elem->GetEdgesNumber();
         elem->GetEdges(e_edges);
         for (ii = 0; ii < nedges; ii++)
         {
            edge = e_edges[ii];
            if (edge->GetMark())
               continue;
            edge->GetNodes(e_nodes);
            // Edge A
            if (*node == *e_nodes[0])
               e_nodes[0]->SetMark(true);
            // Edge B
            if (*node == *e_nodes[1])
               e_nodes[1]->SetMark(true);
            if (msh->getOrder())                  // Quadratic
            {
               if (*node == *e_nodes[2])
                  e_nodes[2]->SetMark(true);
            }
            if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark())
            {
               if (msh->getOrder())
               {
                  if (e_nodes[2]->GetMark())
                     edge->SetMark(true);
               } else
               edge->SetMark(true);
            }

         }                                        // e_edges
      }
   }

   for (i = 0; i < (long) msh->edge_vector.size(); i++)
   {
      edge = msh->edge_vector[i];
      if (!edge->GetMark())
         continue;
      edge->GetNodes(e_nodes);
      if (msh->getOrder())                        // Quad
      {
         if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark()
            && e_nodes[2]->GetMark())
         {
            Jac = 0.5 * edge->Length();
            v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
            v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
            if (Const && (!msh->isAxisymmetry()))
            {
               NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1 / 3.0;
               NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1 / 3.0;
               NVal[G2L[e_nodes[2]->GetIndex()]] += 4.0 * Jac * v1 / 3.0;

            }
            else
            {
               for (k = 0; k < 3; k++)            // Three nodes
               {
                  // Numerical integration
                  for (l = 0; l < 3; l++)         // Gauss points
                  {
                     Weight = Jac * MXPGaussFkt(3, l);
                     eta = MXPGaussPkt(3, l);
                     ShapeFunctionLineHQ(Shfct, &eta);
                     //Axisymmetical problem
                     if (msh->isAxisymmetry())
                     {
                        radius = 0.0;
                        for (ii = 0; ii < 3; ii++)
                           radius += Shfct[ii] * e_nodes[ii]->X();
                        Weight *= radius;         //2.0*pai*radius;
                     }
                     NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
                        + eta * (v2 - v1)) * Shfct[k] * Weight;
                  }
               }
            }
         }
      } else                                      // Linear
      {
         if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark())
         {
            Jac = 0.5 * edge->Length();
            v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
            v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
            if (!msh->isAxisymmetry())
            {
               if (Const)
               {
                  NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1;
                  NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1;
               }
               else
               {
                  NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * (2.0 * v1
                     + v2) / 3.0;
                  NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * (v1 + 2.0
                     * v2) / 3.0;
               }
            } else                                // Axisymmetry
            {

               for (k = 0; k < 2; k++)            // Three nodes
               {
                  // Numerical integration
                  for (l = 0; l < 3; l++)         // Gauss points
                  {
                     Weight = Jac * MXPGaussFkt(3, l);
                     eta = MXPGaussPkt(3, l);
                     ShapeFunctionLine(Shfct, &eta);
                     //Axisymmetical problem
                     if (msh->isAxisymmetry())
                     {
                        radius = 0.0;
                        for (ii = 0; ii < 2; ii++)
                           radius += Shfct[ii] * e_nodes[ii]->X();
                        Weight *= radius;         //2.0*pai*radius;
                     }
                     NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
                        + eta * (v2 - v1)) * Shfct[k] * Weight;
                  }
               }
            }                                     // End of is (!axi)
         }
      }
   }
   for (i = 0; i < this_number_of_nodes; i++)
      node_value_vector[i] = NVal[i];
   for (i = 0; i < (long) msh->edge_vector.size(); i++)
      msh->edge_vector[i]->SetMark(true);
   for (i = 0; i < nSize; i++)
      msh->nod_vector[i]->SetMark(true);
   NVal.clear();
   G2L.clear();
   e_nodes.resize(0);
   e_edges.resize(0);
}


/**************************************************************************
 ROCKFLOW - Funktion: FaceIntegration
 Task: Translate distributed Neumann boundary condition /source term on faces
 found on a surface into nodes value for all kinds of element
 Programming:
 08/2005 WW Re-Implementation
 11/2005 WW/OK Layer optimization
 01/2010 NW improvement of efficiency to search faces
 **************************************************************************/

void CSourceTerm::FaceIntegration(CFEMesh* msh, std::vector<long>&nodes_on_sfc,
		std::vector<double>&node_value_vector)
{
   if (!msh)
   {
      std::cout
         << "Warning in CSourceTerm::FaceIntegration: no MSH data, function doesn't function";
      return;
   }

   long i, j, k, l;
   long this_number_of_nodes;
   int nfaces, nfn;
   int nodesFace[8];
   double nodesFVal[8];

   bool Const = false;
   if (this->getProcessDistributionType() == FiniteElement::CONSTANT || this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
      //	if (dis_type_name.find("CONSTANT") != std::string::npos)
      Const = true;
   //----------------------------------------------------------------------
   // Interpolation of polygon values to nodes_on_sfc
   if (!Const)                                    // Get node BC by interpolation with surface
   {
      int nPointsPly = 0;
      double Area1, Area2;
      double Tol = 1.0e-9;
      bool Passed;
      const int Size = (int) nodes_on_sfc.size();
      double gC[3], p1[3], p2[3], pn[3], vn[3], unit[3], NTri[3];

      CGLPolyline* m_polyline = NULL;
      Surface *m_surface = NULL;
      m_surface = GEOGetSFCByName(geo_name);      //CC

      // list<CGLPolyline*>::const_iterator p = m_surface->polyline_of_surface_list.begin();
      std::vector<CGLPolyline*>::iterator p =
         m_surface->polyline_of_surface_vector.begin();

      for (j = 0; j < Size; j++)
      {
         pn[0] = msh->nod_vector[nodes_on_sfc[j]]->X();
         pn[1] = msh->nod_vector[nodes_on_sfc[j]]->Y();
         pn[2] = msh->nod_vector[nodes_on_sfc[j]]->Z();
         node_value_vector[j] = 0.0;
         Passed = false;
         // nodes close to first polyline
         p = m_surface->polyline_of_surface_vector.begin();
         while (p != m_surface->polyline_of_surface_vector.end())
         {
            m_polyline = *p;
            // Grativity center of this polygon
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

            // Area of this polygon by the grativity center
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
                  // Intopolation whin triangle (p1,p2,gC)
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
         }                                        // while
      }                                           //j
   }

   int Axisymm = 1;                               // ani-axisymmetry
   //CFEMesh* msh = m_pcs->m_msh;
   if (msh->isAxisymmetry())
      Axisymm = -1;                               // Axisymmetry is true
   CElem* elem = NULL;
   CElem* face = new CElem(1);
   CElement* fem = new CElement(Axisymm * msh->GetCoordinateFlag());
   CNode* e_node = NULL;
   CElem* e_nei = NULL;
   //vec<CNode*> e_nodes(20);
   // vec<CElem*> e_neis(6);

   face->SetFace();
   this_number_of_nodes = (long) nodes_on_sfc.size();
   int nSize = (long) msh->nod_vector.size();
   std::vector<long> G2L(nSize);
   std::vector<double> NVal(this_number_of_nodes);

   for (i = 0; i < nSize; i++)
   {
      msh->nod_vector[i]->SetMark(false);
      G2L[i] = -1;
   }

   for (i = 0; i < this_number_of_nodes; i++)
   {
      NVal[i] = 0.0;
      k = nodes_on_sfc[i];
      G2L[k] = i;
   }

   //----------------------------------------------------------------------
   // NW 15.01.2010
   // 1) search element faces on the surface
   // 2) face integration

   //init
   for (i = 0; i < (long) msh->ele_vector.size(); i++)
   {
      msh->ele_vector[i]->selected = 0;           //TODO can use a new variable
   }
   std::set<long> set_nodes_on_sfc;               //unique set of node id on the surface
   for (i = 0; i < (long) nodes_on_sfc.size(); i++)
   {
      set_nodes_on_sfc.insert(nodes_on_sfc[i]);
   }

   //filtering elements: elements should have nodes on the surface
   //Notice: node-elements relation has to be constructed beforehand
   std::vector<long> vec_possible_elements;
   for (i = 0; i < this_number_of_nodes; i++)
   {
      k = nodes_on_sfc[i];
      for (j = 0; j < (long) msh->nod_vector[k]->connected_elements.size(); j++)
      {
         l = msh->nod_vector[k]->connected_elements[j];
         if (msh->ele_vector[l]->selected == 0)
            vec_possible_elements.push_back(l);
         msh->ele_vector[l]->selected += 1;       // remember how many nodes of an element are on the surface
      }
   }
   //search elements & face integration
   int count;
   double fac = 1.0;
   for (i = 0; i < (long) vec_possible_elements.size(); i++)
   {
      elem = msh->ele_vector[vec_possible_elements[i]];
      if (!elem->GetMark())
         continue;
      nfaces = elem->GetFacesNumber();
      elem->SetOrder(msh->getOrder());
      for (j = 0; j < nfaces; j++)
      {
         e_nei = elem->GetNeighbor(j);
         nfn = elem->GetElementFaceNodes(j, nodesFace);
         //1st check
         if (elem->selected < nfn)
            continue;
         //2nd check: if all nodes of the face are on the surface
         count = 0;
         for (k = 0; k < nfn; k++)
         {
            e_node = elem->GetNode(nodesFace[k]);
            if (set_nodes_on_sfc.count(e_node->GetIndex()) > 0)
            {
               count++;
            }
         }
         if (count != nfn)
            continue;
         // face integration
         for (k = 0; k < nfn; k++)
         {
            e_node = elem->GetNode(nodesFace[k]);
            nodesFVal[k] = node_value_vector[G2L[e_node->GetIndex()]];
         }
         fac = 1.0;
                                                  // Not a surface face
         if (elem->GetDimension() == e_nei->GetDimension())
            fac = 0.5;
         face->SetFace(elem, j);
         face->SetOrder(msh->getOrder());
         face->ComputeVolume();
         fem->setOrder(msh->getOrder() + 1);
         fem->ConfigElement(face, true);
         fem->FaceIntegration(nodesFVal);
         for (k = 0; k < nfn; k++)
         {
            e_node = elem->GetNode(nodesFace[k]);
            NVal[G2L[e_node->GetIndex()]] += fac * nodesFVal[k];
         }
      }
   }

   /*
    //----------------------------------------------------------------------
    int count;
    double fac=1.0;
    for (i = 0; i < (long)msh->ele_vector.size(); i++)
    {
    elem = msh->ele_vector[i];
    if(!elem->GetMark()) continue;
    nfaces = elem->GetFacesNumber();
    elem->SetOrder(msh->getOrder());
    for(j=0; j<nfaces; j++)
   {
   e_nei =  elem->GetNeighbor(j);
   nfn = elem->GetElementFaceNodes(j, nodesFace);
   count=0;
   for(k=0; k<nfn; k++)
   {
   e_node = elem->GetNode(nodesFace[k]);
   for (l = 0; l <this_number_of_nodes; l++)
   {
   if(*e_node==*msh->nod_vector[nodes_on_sfc[l]])
   {
   count++;
   break;
   }
   }
   }
   if(count!=nfn) continue;
   for(k=0; k<nfn; k++)
   {
   e_node = elem->GetNode(nodesFace[k]);
   nodesFVal[k] = node_value_vector[G2L[e_node->GetIndex()]];
   }
   fac = 1.0;
   if(elem->GetDimension()==e_nei->GetDimension()) // Not a surface face
   fac = 0.5;
   face->SetFace(elem, j);
   face->SetOrder(msh->getOrder());
   face->ComputeVolume();
   fem->setOrder(msh->getOrder()+1);
   fem->ConfigElement(face, true);
   fem->FaceIntegration(nodesFVal);
   for(k=0; k<nfn; k++)
   {
   e_node = elem->GetNode(nodesFace[k]);
   NVal[G2L[e_node->GetIndex()]] += fac*nodesFVal[k];
   }
   }
   }
   */

   for (i = 0; i < this_number_of_nodes; i++)
      node_value_vector[i] = NVal[i];
   for (i = 0; i < nSize; i++)
      msh->nod_vector[i]->SetMark(true);

   NVal.clear();
   G2L.clear();
   delete fem;
   delete face;
}


/**************************************************************************
 ROCKFLOW - Funktion: DomainIntegration
 Task:  Translate distributed source term within elements into nodes value
 for all kinds of element
 Programming:
 08/2005 WW Re-Implementation
 09/2010 TF re structured some things
 **************************************************************************/
void CSourceTerm::DomainIntegration(CFEMesh* msh, const std::vector<long>&nodes_in_dom,
std::vector<double>&node_value_vector) const
{
   double nodesFVal[8];

   int Axisymm = 1;                               // ani-axisymmetry
   if (msh->isAxisymmetry())
      Axisymm = -1;                               // Axisymmetry is true
   CElement* fem = new CElement(Axisymm * msh->GetCoordinateFlag());
   vec<CNode*> e_nodes(20);

   const size_t this_number_of_nodes (nodes_in_dom.size());
   const size_t nSize (msh->nod_vector.size());
   std::vector<long> G2L(nSize);
   std::vector<double> NVal(this_number_of_nodes);

   for (size_t i = 0; i < nSize; i++)
   {
      msh->nod_vector[i]->SetMark(false);
      G2L[i] = -1;
   }

   for (size_t i = 0; i < this_number_of_nodes; i++)
   {
      NVal[i] = 0.0;
      G2L[nodes_in_dom[i]] = i;
   }

   size_t count = 0;
   for (size_t i = 0; i < msh->ele_vector.size(); i++)
   {
      CElem* elem (msh->ele_vector[i]);
      if (!elem->GetMark())
         continue;
      elem->GetNodes(e_nodes);
      size_t nn = elem->GetNodesNumber(msh->getOrder());
      count = 0;
      for (size_t j = 0; j < nn; j++)
      {
         for (size_t k = 0; k < this_number_of_nodes; k++)
         {
            if (*e_nodes[j] == *msh->nod_vector[nodes_in_dom[k]])
            {
               count++;
               break;
            }
         }
      }
      if (count != nn)
         continue;
      for (size_t j = 0; j < nn; j++)
         nodesFVal[j] = node_value_vector[G2L[e_nodes[j]->GetIndex()]];
      fem->ConfigElement(elem, true);
      fem->setOrder(msh->getOrder() + 1);
      fem->FaceIntegration(nodesFVal);
      for (size_t j = 0; j < nn; j++)
         NVal[G2L[e_nodes[j]->GetIndex()]] += nodesFVal[j];
   }

   for (size_t i = 0; i < this_number_of_nodes; i++)
      node_value_vector[i] = NVal[i];
   for (size_t i = 0; i < nSize; i++)
      msh->nod_vector[i]->SetMark(true);

   NVal.clear();
   G2L.clear();
   e_nodes.resize(0);
   delete fem;
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void STDelete()
{
   long i;
   int no_st = (int) st_vector.size();
   for (i = 0; i < no_st; i++)
   {
      delete st_vector[i];
   }
   st_vector.clear();
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void STGroupsDelete()
{
   CSourceTermGroup* m_st_group = NULL;
   std::list<CSourceTermGroup*>::const_iterator p = st_group_list.begin();
   while (p != st_group_list.end())
   {
      m_st_group = *p;
      delete m_st_group;
      //st_group_list.remove(*p);
      ++p;
   }
   st_group_list.clear();
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 03/2005 OK Implementation
 last modification:
 **************************************************************************/
// 05/2010 TF not needed for compilation
//void STCreateFromLIN(vector<CGLLine*>lin_properties_vector)
//{
//  long i;
//  CSourceTerm *m_st = NULL;
//  CGLLine *m_lin = NULL;
//  long lin_properties_vector_size = (long)lin_properties_vector.size();
//  for(i=0;i<lin_properties_vector_size;i++){
//    m_st = new CSourceTerm();
//    m_lin = lin_properties_vector[i];
//    m_st->pcs_pv_name = "PRESSURE1"; // ToDo
//    m_st->geo_type_name = "LINE";
//    m_st->setGeoName (m_lin->name);
//    m_st->geo_id = m_lin->gli_line_id;
//    m_st->dis_type_name = "CONSTANT_NEUMANN";
//    m_st->geo_node_value = m_lin->value;
//    m_st->tim_type_name = m_lin->name;
//    st_vector.push_back(m_st);
//  }
//}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 03/2005 OK Implementation
 last modification:
 05/2010 TF restructured a little bit
 **************************************************************************/
CSourceTerm* STGet(std::string geo_name)
{
   size_t st_vector_size(st_vector.size());
   for (size_t i = 0; i < st_vector_size; i++)
   {
      if (st_vector[i]->getGeoName().compare(geo_name) == 0)
         return st_vector[i];
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
CSourceTermGroup* STGetGroup(std::string pcs_type_name, std::string pcs_pv_name)
{
   CSourceTermGroup *m_st_group = NULL;
   std::list<CSourceTermGroup*>::const_iterator p_st_group = st_group_list.begin();
   while (p_st_group != st_group_list.end())
   {
      m_st_group = *p_st_group;
      if (m_st_group->pcs_type_name.compare(pcs_type_name) == 0
         && (m_st_group->pcs_pv_name.compare(pcs_pv_name) == 0))
         return m_st_group;
      ++p_st_group;
   }
   return NULL;
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2006 MB Implementation
 **************************************************************************/
                                                  //WW
double GetConditionalNODValue(CSourceTerm* st, CNodeValue* cnodev)
{
   int nidx;
   double value_cond = 0.0;
   double NodeReachLength;
   CRFProcess* pcs_cond (PCSGet(st->pcs_type_name_cond));
   long node_cond;

   //WW  node_cond = group_vector[i]->msh_node_number_conditional;
                                                  //WW
   node_cond = cnodev->msh_node_number_conditional;
   nidx = pcs_cond->GetNodeValueIndex(st->pcs_pv_name_cond) + 1;
   value_cond = pcs_cond->GetNodeValue(node_cond, nidx);

   if (st->pcs_pv_name_cond.find("FLUX") != std::string::npos)
   {
      //WW    NodeReachLength = group_vector[i]->node_value;
      NodeReachLength = cnodev->node_value;       //WW
      value_cond = value_cond * NodeReachLength;
   }

   return value_cond;
}


/**************************************************************************
 FEMLib-Method:
 Task: Calculate Philips (1957) two-term infiltration flux term
 calculated separately for each node
 Programing:
 05/2007 JOD Implementation
 09/2010 KR cleaned up code
 **************************************************************************/
void GetPhilipNODValue(double &value, const CSourceTerm* m_st)
{
   double infiltration = m_st->constant + m_st->sorptivity / sqrt(aktuelle_zeit);
   infiltration = std::min(m_st->rainfall, infiltration);
   value = infiltration * (-value);
}


/**************************************************************************
 FEMLib-Method:
 Task: Calculate Green-Ampt infiltration flux term
 for homogeneous soil, includes water depth
 writes cumulative infiltration in COUPLINGFLUX
 solution is sensitive to time step
 infiltration is estimated with first order derivative
 of cumulative infiltration
 Programing:
 05/2007 JOD Implementation

**************************************************************************/
void GetGreenAmptNODValue(double &value, CSourceTerm* m_st, long msh_node)
{

   double F, Fiter, Fold, infiltration;
   double conductivity, suction, Theta, wdepth;
   double a, b;
   CFEMesh* m_msh = NULL;
   CRFProcess* m_pcs_this = NULL;
   m_pcs_this = PCSGet(convertProcessTypeToString (m_st->getProcessType()));
   m_msh = m_pcs_this->m_msh;

   double area = value;

   wdepth = std::max(0., m_pcs_this->GetNodeValue(msh_node,
      m_pcs_this->GetNodeValueIndex("HEAD") + 1)
      - m_msh->nod_vector[msh_node]->Z());
   conductivity = m_st->constant;
   suction = m_st->sorptivity + wdepth;           // water depth included
   Theta = m_st->moistureDeficit * suction;

   Fold = m_pcs_this->GetNodeValue(msh_node, m_pcs_this->GetNodeValueIndex(
      "COUPLING"));
   F = Fold;

   do                                             // Newton iteration loop
   {
      Fiter = F;
      if (Fiter == 0)                             // avoids a = 0
         Fiter = 1.e-5;
      a = 1 - Theta / (Fiter + Theta);

      b = Fold - Fiter + Theta * log((Fiter + Theta) / (Fold + Theta))
         + conductivity * dt;                     // dt = timeStep
      F = Fiter + b / a;
   } while (fabs(F - Fiter) > 1.e-10);

   infiltration = (F - Fold) / dt;

   if (infiltration > m_st->rainfall)             // + wdepth / timestep )   // compare with available water
      infiltration = m_st->rainfall;              //  +  wdepth / timestep ;

   F = infiltration * dt + Fold;

   m_pcs_this->SetNodeValue(msh_node,
                                                  // output is cumulative infiltration
      m_pcs_this->GetNodeValueIndex("COUPLING") + 1, F);

   infiltration *= -area;
   value = infiltration;

}


/**************************************************************************
 FEMLib-Method:
 Task: Calculate coupling flux
 Mutual coupling of overland, Richards and groundwater flow
 Programing:
 01/2007 JOD Implementation
 09/2010 KR cleaned up code
 **************************************************************************/
#ifndef NEW_EQS                                   //WW. 06.11.2008
void GetCouplingNODValue(double &value, CSourceTerm* st, CNodeValue* cnodev)
{
   //	if (st->COUPLING_SWITCH == true ||
   //		st->getProcessType() == GROUNDWATER_FLOW ||
   //		st->getProcessType() == RICHARDS_FLOW ||
   //		st->getProcessType() == OVERLAND_FLOW)
   //		GetCouplingNODValuePicard(value, st, cnodev);
   //	else
   //		cout << "Error in GetCouplingNODValue";

   if (st->isCoupled() ||
      st->getProcessType() == GROUNDWATER_FLOW ||
      st->getProcessType() == RICHARDS_FLOW ||
      st->getProcessType() == OVERLAND_FLOW)
      GetCouplingNODValuePicard(value, st, cnodev);
   else
      std::cout << "Error in GetCouplingNODValue";
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task: Calculate coupling flux for GetCouplingNODValue
 for the soil or groundwater flow case
 prerequisite: fluid data in mfp_vector[0] (2 times)
 Programing:
 01/2007 JOD Implementation
 10/2008 JOD overland node shifting for soil columns, averaging automatically 4.7.10
 **************************************************************************/
#ifndef NEW_EQS                                   //WW. 06.11.2008
void GetCouplingNODValuePicard(double &value, CSourceTerm* m_st,
CNodeValue* cnodev)
{

   double relPerm, condArea, gamma;
   CRFProcess* m_pcs_cond = NULL;
   CRFProcess* m_pcs_this = NULL;
   double h_this, h_cond, z_this, z_cond, h_cond_shifted, help;

                                                  // only one phase
   gamma = mfp_vector[0]->Density() * GRAVITY_CONSTANT;
   m_pcs_this = PCSGet(convertProcessTypeToString(m_st->getProcessType()));
   m_pcs_cond = PCSGet(m_st->pcs_type_name_cond);
   GetCouplingFieldVariables(&h_this, &h_cond, &h_cond_shifted, &help,
      &z_this, &z_cond, m_st, cnodev);            // z_cond shifted for soil columns
   ////////   coupling factor
   if (h_this > h_cond_shifted)
      relPerm = GetRelativeCouplingPermeability(m_pcs_this, h_this, m_st,
         cnodev->msh_node_number);                // groundwater or richards
   else
      relPerm = GetRelativeCouplingPermeability(m_pcs_cond, h_cond, m_st,
         cnodev->msh_node_number_conditional);    // overland

   condArea = value * relPerm;                    // leakance in relPerm
   if (m_st->channel)                             // wetted perimeter, pcs_cond has to be overland flow
      condArea *= m_st->channel_width + (h_cond - z_cond);
   ///////   source term & matrix
                                                  // hydrostatic overland pressure
   value = CalcCouplingValue(condArea, h_this, h_cond_shifted, z_cond, m_st);

   if (m_st->getProcessType() == RICHARDS_FLOW)
      condArea /= gamma;

   if (m_st->getProcessType() == GROUNDWATER_FLOW && h_this < z_cond
      && m_st->pcs_type_name_cond == "OVERLAND_FLOW")
      condArea = 0;

                                                  // groundwater / soil pressure
   MXInc(cnodev->msh_node_number, cnodev->msh_node_number, condArea);
   /////  output
   m_pcs_this->SetNodeValue(cnodev->msh_node_number,
                                                  // update coupling variable for error estimation
      m_pcs_this->GetNodeValueIndex("COUPLING") + 1, h_this);

   if (m_st->no_surface_water_pressure)           // 4.10.01
      value = 0;
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task: Calculate coupling flux for GetCouplingNODValue
 for the overland flow case
 prerequisite: fluid data in mfp_vector[0]
 Programing:
 01/2007 JOD Implementation
 10/2008 JOD node shifting for soil columns 4.7.10
 **************************************************************************/
#ifndef NEW_EQS                                   //WW. 06.11.2008
void GetCouplingNODValueNewton(double &value, CSourceTerm* m_st,
CNodeValue* cnodev)
{
   double relPerm, area, condArea;
   CRFProcess* m_pcs_cond = NULL;
   CRFProcess* m_pcs_this = NULL;
   double h_this, h_cond, z_this, z_cond, h_this_shifted, h_this_averaged;
   double epsilon = 1.e-7, value_jacobi, h_this_epsilon = 0.0,
      relPerm_epsilon, condArea_epsilon;          //OK411 epsilon as in pcs->assembleParabolicEquationNewton

   m_pcs_this = PCSGet(convertProcessTypeToString(m_st->getProcessType()));
   m_pcs_cond = PCSGet(m_st->pcs_type_name_cond);
   GetCouplingFieldVariables(&h_this, &h_cond, &h_this_shifted,
                                                  // z_this shifted for soil columns
      &h_this_averaged, &z_this, &z_cond, m_st, cnodev);
   ///// factor
   if (h_this_shifted > h_cond)
      relPerm = GetRelativeCouplingPermeability(m_pcs_this, h_this, m_st,
         cnodev->msh_node_number);                // overland
   else
      relPerm = GetRelativeCouplingPermeability(m_pcs_cond, h_cond, m_st,
         cnodev->msh_node_number_conditional);    // richards or groundwater

   area = value;                                  // for condArea_epsilon
   condArea = value * relPerm;                    // leakance in relPerm
   if (m_st->channel)                             // wetted perimeter
      condArea *= m_st->channel_width + (h_this - z_this);
   //////  source term
   value = CalcCouplingValue(condArea, h_this_averaged, h_cond, z_cond, m_st);
   /////   epsilon shift for jacobian
   if (m_st->node_averaging)
   {
      for (long i = 0; i < (long) cnodev->msh_node_numbers_averaging.size(); i++)
         if (cnodev->msh_node_numbers_averaging[i]
         == cnodev->msh_node_number)
            h_this_epsilon = h_this_averaged + epsilon
               * cnodev->msh_node_weights_averaging[i];
   } else
   h_this_epsilon = h_this + epsilon;
   /////  factor
   if (h_this_shifted + epsilon > h_cond)
      relPerm_epsilon = GetRelativeCouplingPermeability(m_pcs_this, h_this
         + epsilon, m_st, cnodev->msh_node_number);
   else
      relPerm_epsilon = GetRelativeCouplingPermeability(m_pcs_cond, h_cond,
         m_st, cnodev->msh_node_number_conditional);// richards or groundwater

   condArea_epsilon = area * relPerm_epsilon;
   if (m_st->channel)                             // wetted perimeter
      condArea_epsilon *= m_st->channel_width + (h_this_epsilon - z_this);

   if (m_st->no_surface_water_pressure)           // 4.10.01
      h_this_epsilon = z_this;
   /////  jacobian
   value_jacobi = -condArea_epsilon * (h_cond - h_this_epsilon) + value;
   MXInc(cnodev->msh_node_number, cnodev->msh_node_number, value_jacobi
      / epsilon);
   /////  output
   m_pcs_this->SetNodeValue(cnodev->msh_node_number,
                                                  // coupling flux (m/s)
      m_pcs_this->GetNodeValueIndex("COUPLING") + 1, -value / area);

}
#endif
/**************************************************************************
 FEMLib-Method:
 Task: Calculate relative coupling permeability for GetCouplingNODValue(***).
 Programing: prerequisites: phase 0 in mfp
 06/2007 JOD Implementation
 10/2008 JOD include leakance 4.7.10
 **************************************************************************/
double GetRelativeCouplingPermeability(const CRFProcess* pcs, double head, const CSourceTerm* source_term, long msh_node)
{
   double relPerm;
   double z = pcs->m_msh->nod_vector[msh_node]->Z();

   //	if (pcs->pcs_type_name == "OVERLAND_FLOW") {
   if (pcs->getProcessType() == OVERLAND_FLOW)
   {
      double sat = (head - z) / std::max(1.e-6, source_term->rill_height);
      if (sat > 1)
         relPerm = 1;
      else if (sat < 0)
         relPerm = 0;
      else
         relPerm = pow(sat, 2* (1 - sat));
   }
   //	else if( pcs->pcs_type_name == "GROUNDWATER_FLOW")
   else if( pcs->getProcessType() == GROUNDWATER_FLOW)
      relPerm = 1;
   //	else if( pcs->pcs_type_name == "RICHARDS_FLOW") {
   else if( pcs->getProcessType() == RICHARDS_FLOW)
   {
      /*CElem *m_ele = NULL;
       long msh_ele;
       int group;
       double gamma =  mfp_vector[0]->Density() * GRAVITY_CONSTANT; // only one phase

       msh_ele = pcs->m_msh->nod_vector[msh_node]->connected_elements[0];
       m_ele = pcs->m_msh->ele_vector[msh_ele];
       group = pcs->m_msh->ele_vector[msh_ele]->GetPatchIndex();

       sat = mmp_vector[group]->SaturationCapillaryPressureFunction( -(head - z) * gamma, 0);
       relPerm = mmp_vector[group]->PermeabilitySaturationFunction(sat,0);*/

      relPerm = 1;                                // JOD 4.10.01
   }
   else
   {
      std::cout <<"!!!!! Coupling flux upwinding not implemented for this case !!!!!" << std::endl;
      relPerm = 1;
   }

   return relPerm * source_term->getCoupLeakance();
}


/**************************************************************************
 FEMLib-Method:
 Task: Coupling of overland and soil flow by using water depth as soil boundary
 condition and flux term as overland source term according to
 Morita and Yen, J. Hydr. Eng. 184, 2002
 Programing: prerequisites: constant precipitation with assigned duration,
 phase = 0 in mfp, soil data in mmp_vetor[1] !!!!!
 06/2007 JOD Implementation
 **************************************************************************/
#ifndef NEW_EQS                                   //WW. 06.11.2008
void GetCouplingNODValueMixed(double& value, CSourceTerm* m_st,
CNodeValue* cnodev)
{

   double cond1, cond0, pressure1, pressure0, bc_value, depth, gamma, sat,
      area;
   double leakance, deltaZ;
   int phase = 0;                                 // RESTRICTION for mfp !!!!!!!

   CElem *m_ele = NULL;
   long msh_ele;
   int group, nidx;
   CRFProcess* m_pcs_cond = NULL;
   CRFProcess* m_pcs_this = NULL;
   m_pcs_this = PCSGet(convertProcessTypeToString (m_st->getProcessType()));
   m_pcs_cond = PCSGet(m_st->pcs_type_name_cond);

   area = value;
   leakance = m_st->getCoupLeakance();
   deltaZ = m_st->rill_height;
                                                  // phase  = 0 !!!!
   gamma = mfp_vector[0]->Density() * GRAVITY_CONSTANT;
   long msh_node_2nd;
   double x_this = m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->X();
   double y_this = m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->Y();
   double z_this = m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->Z();

   msh_node_2nd = -1;                             //WW

   cond0 = leakance * deltaZ;
   cond1 = cond0;

   if (m_st->getProcessType () == OVERLAND_FLOW)
   {

      ///// get number of second mesh node, provisional implementation
      double epsilon = 1.e-5;
      double
         x_cond =
         m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->X();
      double
         y_cond =
         m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->Y();
      double
         z_cond =
         m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->Z();

      for (int i = 0; i < (long) m_pcs_cond->m_msh->nod_vector.size(); i++)
      {
         if (m_pcs_cond->m_msh->nod_vector[i]->X() - x_cond < epsilon)
            if (m_pcs_cond->m_msh->nod_vector[i]->Y() - y_cond < epsilon)
               if (m_pcs_cond->m_msh->nod_vector[i]->Z() - (z_cond
            - deltaZ) < epsilon)
                  msh_node_2nd = i;
      }
      //////////////////////////

      nidx = m_pcs_cond->GetNodeValueIndex("PRESSURE1") + 1;

      pressure0 = m_pcs_cond->GetNodeValue(
         cnodev->msh_node_number_conditional, nidx);
      pressure1 = m_pcs_cond->GetNodeValue(msh_node_2nd, nidx);

                                                  // only one phase
      double gamma = mfp_vector[phase]->Density() * GRAVITY_CONSTANT;

      msh_ele
         = m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->connected_elements[0];
      m_ele = m_pcs_cond->m_msh->ele_vector[msh_ele];
      group = m_pcs_cond->m_msh->ele_vector[msh_ele]->GetPatchIndex();

      //sat = mmp_vector[group]->SaturationCapillaryPressureFunction( -pressure0, 0);
      //cond0 *=  mmp_vector[group]->PermeabilitySaturationFunction(sat,0);

      sat = mmp_vector[group]->SaturationCapillaryPressureFunction(
         -pressure1, phase);
      cond1 *= mmp_vector[group]->PermeabilitySaturationFunction(sat, phase);
      // use of relative permeability for second node (absolute perm. for top node !!!!)

      value = (pressure1 - pressure0 - deltaZ * gamma) * (cond0 + cond1)
         / (2* deltaZ * gamma);

      m_pcs_this->SetNodeValue(cnodev->msh_node_number,
         m_pcs_this->GetNodeValueIndex("COUPLING") + 1, -value);

      value *= area;
   }                                              // end overland
   else                                           // Richards
   {

      ///// get number of second mesh node, provisional implementation
      double epsilon = 1.e-5;
      for (int i = 0; i < (long) m_pcs_this->m_msh->nod_vector.size(); i++)
      {
         if (m_pcs_this->m_msh->nod_vector[i]->X() - x_this < epsilon)
            if (m_pcs_this->m_msh->nod_vector[i]->Y() - y_this < epsilon)
               if (m_pcs_this->m_msh->nod_vector[i]->Z() - (z_this
            - deltaZ) < epsilon)
                  msh_node_2nd = i;
      }
      //////////////////////////

      double inf_cap, supplyRate, rainfall;
      long
         bc_eqs_index =
         m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->GetEquationIndex();
      double
         z_cond =
         m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->Z();
      depth = std::max(0., m_pcs_cond->GetNodeValue(
         cnodev->msh_node_number_conditional,
         m_pcs_cond->GetNodeValueIndex("HEAD") + 1) - z_cond);

      nidx = m_pcs_this->GetNodeValueIndex("PRESSURE1") + 1;
      pressure0 = m_pcs_this->GetNodeValue(cnodev->msh_node_number, nidx);
      pressure1 = m_pcs_this->GetNodeValue(msh_node_2nd, nidx);

      msh_ele
         = m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->connected_elements[0];
      m_ele = m_pcs_this->m_msh->ele_vector[msh_ele];
      group = m_pcs_this->m_msh->ele_vector[msh_ele]->GetPatchIndex();

      //sat = mmp_vector[group]->SaturationCapillaryPressureFunction( -pressure0, phase);
      //cond0 *=  mmp_vector[group]->PermeabilitySaturationFunction(sat,phase);

      sat = mmp_vector[group]->SaturationCapillaryPressureFunction(
         -pressure1, phase);
      cond1 *= mmp_vector[group]->PermeabilitySaturationFunction(sat, phase);
      // use of relative permeability for second node (absolute perm. for top node !!!!)

      // calculate infiltration capacity
      if (aktuelle_zeit < m_st->rainfall_duration)
         rainfall = m_st->rainfall;
      else
         rainfall = 0;

      inf_cap = (depth + deltaZ - pressure1 / gamma) * (cond0 + cond1) / (2
         * deltaZ);
      supplyRate = m_st->rainfall;                //+ (depth ) / dt; // dt = timeStep

      m_pcs_this->SetNodeValue(cnodev->msh_node_number,
                                                  // update coupling variable for error estimation
         m_pcs_this->GetNodeValueIndex("COUPLING") + 1, inf_cap);

      if (inf_cap > supplyRate)
         bc_value = pressure1 - deltaZ * gamma + gamma * supplyRate * deltaZ
            * 2 / (cond0 + cond1);
      else
         bc_value = pressure1 - deltaZ * gamma + gamma * inf_cap * deltaZ
            * 2 / (cond0 + cond1);
      // bc_value = supplyRate * gamma * dt;

      MXRandbed(bc_eqs_index, bc_value, m_pcs_this->eqs->b);
      value = 0;

   }                                              // end Richards

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2006 MB Implementation
 02/2006 WW Change argument
 09/2010 KR cleaned up code
 09/2010 TF commented out method
 **************************************************************************/
//double CSourceTermGroup::GetRiverNODValue(int i,CSourceTerm* m_st, long msh_node) //WW
//void GetRiverNODValue(double &value, CNodeValue* cnodev, const CSourceTerm* m_st) //WW
//{
//	double h;
//	double paraA(cnodev->node_parameterA); //HRiver
//	double paraB(cnodev->node_parameterB); //KRiverBed
//	double paraC(cnodev->node_parameterC); //WRiverBed
//	double paraD(cnodev->node_parameterD); //TRiverBed
//	double paraE(cnodev->node_parameterE); //BRiverBed
//	double NodeReachLength(value);
//	CRFProcess* m_pcs_this(NULL);
//
//	m_pcs_this = PCSGet(convertProcessTypeToString (m_st->getProcessType()));
//	//_________________________________________________________________
//	//paraA jetzt aus dem Prozess Overland Flow
//	if (m_st->isCoupled()) {
//
//		CRFProcess* m_pcs_cond = NULL;
//		m_pcs_cond = PCSGet(m_st->pcs_type_name_cond);
//
//		int nidx = m_pcs_cond->GetNodeValueIndex(m_st->pcs_pv_name_cond) + 1;
//		//WW    long node_cond = group_vector[i]->msh_node_number_conditional;
//		long node_cond = cnodev->msh_node_number_conditional; //WW
//		paraA = m_pcs_cond->GetNodeValue(node_cond, nidx);
//	}
//
//	double RiverConductance = paraB * paraC * NodeReachLength / (paraD - paraE);
//	int nidx1 = m_pcs_this->GetNodeValueIndex("HEAD") + 1;
//	h = m_pcs_this->GetNodeValue(cnodev->msh_node_number, nidx1);
//
//	if (h > paraD) { //HAquiver > BRiverBed
//		//q = (RiverConductance * HRiver)   -  (RiverConductance * HAquifer)
//		value = RiverConductance * paraA;
//		MXInc(cnodev->msh_node_number, cnodev->msh_node_number,	RiverConductance);
//	}
//	if (h < paraD) { //HAquiver < BRiverBed
//		//q = (RiverConductance * HRiver)   - (RiverConductance * BRiverBed)
//		value = RiverConductance * (paraA - paraD);
//	}
//	if (h == paraE)
//		value = 0.;
//	//_________________________________________________________________
//	//Safe Flux values
//	int nidxFLUX = m_pcs_this->GetNodeValueIndex("FLUX") + 1;
//	double flux = value / NodeReachLength; //fluxes in m^2/s !!
//	m_pcs_this->SetNodeValue(cnodev->msh_node_number, nidxFLUX, flux);
//}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2006 MB Implementation
 02/2006 WW Change argument
 **************************************************************************/
//double CSourceTermGroup::GetCriticalDepthNODValue(CNodeValue* cnodev,CSourceTerm* m_st, long msh_node)
void GetCriticalDepthNODValue(double &value, CSourceTerm* m_st, long msh_node)
{
   double value_jacobi;
   double width, flowdepth, flowdepth3, flowdepth3_epsilon;
   long msh_ele;
   double epsilon = 1.e-7;                        // like in pcs->assembleParabolicEquationNewton

   CRFProcess* m_pcs_this = NULL;
   m_pcs_this = PCSGet(convertProcessTypeToString(m_st->getProcessType()));
   long nidx1 = m_pcs_this->GetNodeValueIndex("HEAD") + 1;
   flowdepth = m_pcs_this->GetNodeValue(msh_node, nidx1)
      - m_pcs_this->m_msh->nod_vector[msh_node]->Z() - m_st->rill_height;

   if (flowdepth < 0.0)
   {
      value = 0;
      m_pcs_this->SetNodeValue(msh_node,
         m_pcs_this->GetNodeValueIndex("FLUX") + 0, -value);
   }
   else
   {
      flowdepth3 = pow(flowdepth, 3.);
      flowdepth3_epsilon = pow(flowdepth + epsilon, 3.);
      width = value;
      if (m_pcs_this->m_msh->GetMaxElementDim() == 1)
      {
         msh_ele
            = m_pcs_this->m_msh->nod_vector[msh_node]->connected_elements[0];
         int group = m_pcs_this->m_msh->ele_vector[msh_ele]->GetPatchIndex();
         width = mmp_vector[group]->overland_width;
      }

      value = -sqrt(GRAVITY_CONSTANT * flowdepth3) * width;

      value_jacobi = sqrt(GRAVITY_CONSTANT * flowdepth3_epsilon) * width
         + value;
                                                  // write source term into jacobi
      MXInc(msh_node, msh_node, value_jacobi / epsilon);

      m_pcs_this->SetNodeValue(msh_node,
         m_pcs_this->GetNodeValueIndex("FLUX") + 0, -value);

   }
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2006 MB JOD Implementation
 06/2007 JOD 2D case with slope in st-file
 **************************************************************************/
void GetNormalDepthNODValue(double &value, CSourceTerm* st, long msh_node)
{
   //WW  int AnzNodes = 0;
   //WW  double Haverage = 0;
   CRFProcess* pcs_this (PCSGet(st->getProcessType()));
   CFEMesh* mesh (pcs_this->m_msh);

   double value_for_jacobi, S_0;
   double epsilon = 1.e-7;                        // pcs->assembleParabolicEquationNewton !!!!!!!!!

   long msh_ele = mesh->nod_vector[msh_node]->connected_elements[0];
   CElem *m_ele = mesh->ele_vector[msh_ele];
   int group = mesh->ele_vector[msh_ele]->GetPatchIndex();
   double width = mmp_vector[group]->overland_width;
   double fric_coef = mmp_vector[group]->friction_coefficient;
   double slope_exp = mmp_vector[group]->friction_exp_slope;
   double depth_exp = mmp_vector[group]->friction_exp_depth;
   if (st->getNormalDepthSlope() == -1)           // TF: WARNING comparison of double via ==
   {
      if (mesh->GetMaxElementDim() > 1)
         std::cout << "!!!!! give slope for NORMAL DEPTH in st-file !!!!!"
            << std::endl;

      double elementlength = sqrt(pow(m_ele->GetNode(1)->X()
         - m_ele->GetNode(0)->X(), 2.) + pow(m_ele->GetNode(1)->Y()
         - m_ele->GetNode(0)->Y(), 2.) + pow(m_ele->GetNode(1)->Z()
         - m_ele->GetNode(0)->Z(), 2.));
      S_0 = (m_ele->GetNode(1)->Z() - m_ele->GetNode(0)->Z()) / elementlength;
      if (S_0 < 0)
         S_0 = -S_0;
   } else
   S_0 = st->getNormalDepthSlope();

   double flowdepth = pcs_this->GetNodeValue(msh_node, 1)
      - mesh->nod_vector[msh_node]->Z() - st->rill_height;
   double flowdepth_epsilon = flowdepth + epsilon;
   if (flowdepth < 0.0)
   {
      flowdepth = 0.0;
      flowdepth_epsilon = 0.0;
   }

   double temp = width * fric_coef * pow(S_0, slope_exp);
   if (mmp_vector[group]->channel == 1)
   {
      value = -pow(flowdepth * width / (2 * flowdepth + width), depth_exp)
         * flowdepth * temp;
      value_for_jacobi = pow(flowdepth_epsilon * width / (2
         * flowdepth_epsilon + width), depth_exp) * flowdepth_epsilon
         * temp + value;
   }
   else
   {
      value = -pow(flowdepth, depth_exp + 1) * temp;
      value_for_jacobi = pow(flowdepth_epsilon, depth_exp + 1) * temp + value;
   }

                                                  // write source term into jacobi
   MXInc(msh_node, msh_node, value_for_jacobi / epsilon);
   pcs_this->SetNodeValue(msh_node, pcs_this->GetNodeValueIndex("FLUX")
      + 0, -value);
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
void GetNODValue(double& value, CNodeValue* cnodev, CSourceTerm* st)
{

#ifndef NEW_EQS                                //WW. 06.11.2008
   if (st->isCoupled())
      GetCouplingNODValue(value, st, cnodev);
   else if (st->isAnalytical())
   {
      //WW      m_st_group->m_msh = m_msh;
                                                  //WW
      value = st->GetAnalyticalSolution(cnodev->msh_node_number);
      //WW         value = m_st_group->GetAnalyticalSolution(m_st,msh_node,(string)function_name[j]);
   }

   //	if (cnodev->node_distype == 5) // River Condition
   //	if (cnodev->getProcessDistributionType() == RIVER)
   //		GetRiverNODValue(value, cnodev, st); //MB
   //	if (cnodev->node_distype == 6) // CriticalDepth Condition
                                                  // CriticalDepth Condition
   if (cnodev->getProcessDistributionType() == FiniteElement::CRITICALDEPTH)
                                                  //MB
      GetCriticalDepthNODValue(value, st, cnodev->msh_node_number);
   //	if (cnodev->node_distype == 8) // NormalDepth Condition JOD
   if (cnodev->getProcessDistributionType() == FiniteElement::NORMALDEPTH)
                                                  //MB
      GetNormalDepthNODValue(value, st, cnodev->msh_node_number);
#endif
   //	if (cnodev->node_distype == 10) // Philip infiltration JOD
   //		GetPhilipNODValue(value, st);
   //	if (cnodev->node_distype == 11) // Green_Ampt infiltration JOD
   if (cnodev->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
      GetGreenAmptNODValue(value, st, cnodev->msh_node_number);

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void STGroupDelete(std::string pcs_type_name, std::string pcs_pv_name)
{
   CSourceTermGroup* m_st_group = NULL;
   std::list<CSourceTermGroup*>::const_iterator p = st_group_list.begin();
   while (p != st_group_list.end())
   {
      m_st_group = *p;
      if ((m_st_group->pcs_type_name.compare(pcs_type_name) == 0)
         && (m_st_group->pcs_pv_name.compare(pcs_pv_name) == 0))
      {
         delete m_st_group;
         st_group_list.remove(m_st_group);
         return;
      }
      ++p;
   }
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 modification:
 05/2010
 **************************************************************************/
void CSourceTermGroup::SetPNT(CRFProcess* pcs, CSourceTerm* st,
const int ShiftInNodeVector)
{
   CNodeValue *nod_val (new CNodeValue());

   // TF removed some checks - check validity of data while reading data

   nod_val->msh_node_number = m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(st->getGeoObj())) + ShiftInNodeVector;

   nod_val->CurveIndex = st->CurveIndex;
                                                  //WW
   nod_val->geo_node_number = nod_val->msh_node_number - ShiftInNodeVector;
   nod_val->node_value = st->geo_node_value;
   nod_val->tim_type_name = st->tim_type_name;

   if (st->getProcessDistributionType() == FiniteElement::CRITICALDEPTH)
   {
      //	if (st->dis_type_name.compare("CRITICALDEPTH") == 0) {
      nod_val->setProcessDistributionType (st->getProcessDistributionType());
      nod_val->node_area = 1.0;
   }

   if (st->getProcessDistributionType() == FiniteElement::NORMALDEPTH)
   {
      nod_val->setProcessDistributionType (st->getProcessDistributionType());
      nod_val->node_area = 1.0;
   }

   //	if (st->dis_type_name.compare("PHILIP") == 0) { // JOD
   //		nod_val->node_distype = 10;
   //		nod_val->node_area = 1.0;
   //	}

   if (st->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
   {
      nod_val->setProcessDistributionType (st->getProcessDistributionType());
      nod_val->node_area = 1.0;
   }

   if (st->getProcessDistributionType() == FiniteElement::SYSTEM_DEPENDENT)
   {
      nod_val->setProcessDistributionType (st->getProcessDistributionType());
      pcs->compute_domain_face_normal = true;     //WW
      CElem* elem = NULL;
      CNode* cnode = NULL;                        //WW
      for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
      {
         elem = m_msh->ele_vector[i];
         if (!elem->GetMark())
            continue;
         int nn = elem->GetNodesNumber(m_msh->getOrder());
         for (long j = 0; j < nn; j++)
         {
            cnode = elem->GetNode(j);             //WW
            if (cnode->GetIndex() == (size_t)st->geo_node_number)
               st->element_st_vector.push_back(i);
         }
      }
   }
   //WW        group_vector.push_back(m_node_value);
   //WW        st_group_vector.push_back(st); //OK
   pcs->st_node_value.push_back(nod_val);         //WW
   pcs->st_node.push_back(st);                    //WW
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
void CSourceTermGroup::SetLIN(CRFProcess* m_pcs, CSourceTerm* m_st,
const int ShiftInNodeVector)
{
   (void)m_pcs;
   (void)m_st;
   (void)ShiftInNodeVector;
   /*OK411
    long number_of_nodes;
    vector<long>lin_nod_vector;
    vector<double>lin_nod_val_vector;
    CGLLine* m_lin = NULL;
    CGLPolyline* m_ply = NULL;
    long *nodes = NULL;
    m_lin = m_lin->GEOGetLine(m_st->geo_id);

    if(m_lin){
    double* coordinates;
   m_ply = new CGLPolyline;
   m_ply->point_vector.push_back(m_lin->m_point1);
   m_ply->point_vector.push_back(m_lin->m_point2);
   nodes = MSHGetNodesClose(&number_of_nodes, m_ply);//CC
   lin_nod_val_vector.resize(number_of_nodes);
   for(long i = 0; i < number_of_nodes; i++){
   lin_nod_val_vector[i] =  m_st->geo_node_value / number_of_nodes;
   coordinates = new double[3];
   coordinates[0] = GetNodeX(nodes[i]);
   coordinates[1] = GetNodeY(nodes[i]);
   coordinates[2] = GetNodeZ(nodes[i]);
   m_lin->nodes_coor_vector.push_back(coordinates);
   }
   //InterpolationAlongPolyline(m_polyline,node_value_vector);
   for(long i=0; i < number_of_nodes; i++){
   CNodeValue* m_nod_val = NULL;
   m_nod_val = new CNodeValue();
   m_nod_val->msh_node_number = -1;
   m_nod_val->msh_node_number = nodes[i]+ShiftInNodeVector;
   m_nod_val->geo_node_number = nodes[i];
   m_nod_val->node_value = lin_nod_val_vector[i];
   m_nod_val->CurveIndex = m_st->CurveIndex;
   //WW        group_vector.push_back(m_node_value);
   //WW        st_group_vector.push_back(m_st); //OK
   m_pcs->st_node_value.push_back(m_nod_val);  //WW
   m_pcs->st_node.push_back(m_st); //WW
   }
   lin_nod_val_vector.clear();
   m_ply->point_vector.clear();
   delete m_ply;
   }
   else
   cout << "Warning - CSourceTermGroup::Set: LIN not found" << endl;
   */
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing
 07/2005 OK Implementation based on CSourceTermGroup::Set
 modified
 07/2010 TF substituted GEOGetPLYByName
 **************************************************************************/
void CSourceTermGroup::SetPLY(CSourceTerm* st, int ShiftInNodeVector)
{
   std::vector<long> ply_nod_vector;
   std::vector<long> ply_nod_vector_cond;
   std::vector<double> ply_nod_val_vector;
   CGLPolyline* m_ply = NULL;

   m_ply = GEOGetPLYByName(st->geo_name);

   if (m_ply)
   {
      SetPolylineNodeVector(m_ply, ply_nod_vector);
      if (st->_coupled)
         SetPolylineNodeVectorConditional(st, m_ply, ply_nod_vector,
            ply_nod_vector_cond);
      SetPolylineNodeValueVector(st, m_ply, ply_nod_vector,
         ply_nod_vector_cond, ply_nod_val_vector);
      st->SetNodeValues(ply_nod_vector, ply_nod_vector_cond,
         ply_nod_val_vector, ShiftInNodeVector);
   }                                              // end polyline

   ////	std::vector<size_t> ply_nod_vector;
   ////	std::vector<size_t> ply_nod_vector_cond;
   ////	std::vector<double> ply_nod_val_vector;
   //
   //	std::vector<long> long_ply_nod_vector;
   //	std::vector<long> long_ply_nod_vector_cond;
   //	std::vector<double> ply_nod_val_vector1;
   //
   //	CGLPolyline* old_ply;
   //	old_ply = GEOGetPLYByName(st->geo_name);
   //
   ////	double min_edge_length (m_msh->getMinEdgeLength());
   ////	m_msh->setMinEdgeLength (old_ply->epsilon);
   ////	m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(st->getGeoObj()), ply_nod_vector);
   ////	std::sort (ply_nod_vector.begin(), ply_nod_vector.end());
   ////	m_msh->setMinEdgeLength (min_edge_length);
   //
   //	if (old_ply) {
   //		SetPolylineNodeVector(old_ply, long_ply_nod_vector);  // this method has side effects!!!
   //
   //////		long_ply_nod_vector.clear();
   //////		for (size_t k(0); k<ply_nod_vector.size(); k++) long_ply_nod_vector.push_back(ply_nod_vector[k]);
   ////
   ////		for (size_t k(0); k<long_ply_nod_vector.size(); k++) {
   ////			if (long_ply_nod_vector[k] != ply_nod_vector[k]) {
   ////				std::cout << "points of old polyline: " << std::endl;
   ////				for (size_t k(0); k<old_ply->point_vector.size(); k++)
   ////					std::cout << "id: " << old_ply->point_vector[k]->id << ", propert: " << old_ply->point_vector[k]->propert << std::endl;
   ////				std::cout << std::endl;
   ////				std::cout << "CSourceTermGroup::SetPLY for polyline " << old_ply << " (eps = " << old_ply->epsilon << ") - long_ply_nod_vector: " << std::endl;
   ////				for (size_t k(0); k<long_ply_nod_vector.size(); k++) std::cout << long_ply_nod_vector[k] << " " << std::flush;
   ////				std::cout << std::endl;
   ////
   ////				std::cout << "points of new polyline: " << std::endl;
   ////				for (size_t k(0); k<static_cast<const GEOLIB::Polyline*>(st->getGeoObj())->getNumberOfPoints(); k++)
   ////					std::cout << "id: " << static_cast<const GEOLIB::Polyline*>(st->getGeoObj())->getPointID(k) << " ";
   ////				std::cout << std::endl;
   ////				std::cout << "CSourceTermGroup::SetPLY for polyline " <<static_cast<const GEOLIB::Polyline*>(st->getGeoObj()) << " - ply_nod_vector: " << std::endl;
   ////				for (size_t k(0); k<ply_nod_vector.size(); k++) std::cout << ply_nod_vector[k] << " " << std::flush;
   ////				std::cout << std::endl;
   ////
   ////				exit (1);
   ////			}
   ////		}
   //
   //		if (st->isCoupled()) {
   //			SetPolylineNodeVectorConditional(st, old_ply, long_ply_nod_vector,
   //					long_ply_nod_vector_cond);
   //////			long_ply_nod_vector.clear();
   ////			SetPolylineNodeVectorConditional(st, ply_nod_vector, ply_nod_vector_cond);
   //////			for (size_t k(0); k<ply_nod_vector.size(); k++) long_ply_nod_vector.push_back(ply_nod_vector[k]);
   ////
   ////			for (size_t k(0); k<long_ply_nod_vector_cond.size(); k++) {
   ////				if (long_ply_nod_vector_cond[k] != ply_nod_vector_cond[k]) {
   ////					std::cout << "CSourceTermGroup::SetPLY long_ply_nod_vector_cond: " << std::endl;
   ////					for (size_t k(0); k<long_ply_nod_vector_cond.size(); k++) std::cout << long_ply_nod_vector_cond[k] << " " << std::flush;
   ////					std::cout << std::endl;
   ////					std::cout << "CSourceTermGroup::SetPLY ply_nod_vector_cond: " << std::endl;
   ////					for (size_t k(0); k<ply_nod_vector_cond.size(); k++) std::cout << ply_nod_vector_cond[k] << " " << std::flush;
   ////					std::cout << std::endl;
   ////
   ////					exit (1);
   ////				}
   ////			}
   //		}
   //
   //		SetPolylineNodeValueVector(st, old_ply, long_ply_nod_vector,
   //				long_ply_nod_vector_cond, ply_nod_val_vector1);
   //
   ////		const Mesh_Group::MeshNodesAlongPolyline& msh_nodes_along_polyline (m_msh->GetMeshNodesAlongPolyline (static_cast<const GEOLIB::Polyline*>(st->getGeoObj())));
   ////
   ////		for (size_t k(0); k<long_ply_nod_vector.size(); k++) {
   ////			if (long_ply_nod_vector[k] != ply_nod_vector[k]) {
   ////				std::cout << "points of old polyline: " << std::endl;
   ////				for (size_t k(0); k<old_ply->point_vector.size(); k++)
   ////					std::cout << old_ply->point_vector[k]->id << " ";
   ////				std::cout << std::endl;
   ////				std::cout << "CSourceTermGroup::SetPLY for polyline " << old_ply << " (eps = " << old_ply->epsilon << ") - long_ply_nod_vector: " << std::endl;
   ////				for (size_t k(0); k<long_ply_nod_vector.size(); k++) std::cout << long_ply_nod_vector[k] << " " << std::flush;
   ////				std::cout << std::endl;
   ////
   ////				std::cout << "points of new polyline: " << std::endl;
   ////				for (size_t k(0); k<static_cast<const GEOLIB::Polyline*>(st->getGeoObj())->getNumberOfPoints(); k++)
   ////					std::cout << static_cast<const GEOLIB::Polyline*>(st->getGeoObj())->getPointID(k) << " ";
   ////				std::cout << std::endl;
   ////				std::cout << "CSourceTermGroup::SetPLY for polyline " <<static_cast<const GEOLIB::Polyline*>(st->getGeoObj()) << " - ply_nod_vector: " << std::endl;
   ////				for (size_t k(0); k<ply_nod_vector.size(); k++) std::cout << ply_nod_vector[k] << " " << std::flush;
   ////				std::cout << std::endl;
   ////
   ////				exit (1);
   ////			}
   ////		}
   ////
   //////		SetPolylineNodeValueVector(st, msh_nodes_along_polyline, long_ply_nod_vector_cond, ply_nod_val_vector);
   ////
   //////		std::cout << "CSourceTermGroup::SetPLY ply_nod_val_vector1 and ply_nod_val_vector: " << std::endl;
   //////		if (ply_nod_val_vector1.size() != ply_nod_val_vector.size()) {
   //////			std::cout << ply_nod_val_vector1.size() << " != " << ply_nod_val_vector.size() << std::endl;
   //////			exit (1);
   //////		}
   //////		for (size_t k(0); k<ply_nod_val_vector1.size(); k++)
   //////			std::cout << ply_nod_val_vector1[k] << "\t" << ply_nod_val_vector[k] << std::endl;
   //
   //		st->SetNodeValues(long_ply_nod_vector, long_ply_nod_vector_cond,
   //				ply_nod_val_vector1, ShiftInNodeVector);
   //	} // end polyline
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
void CSourceTermGroup::SetDMN(CSourceTerm *m_st, const int ShiftInNodeVector)
{
   long number_of_nodes;
   std::vector<long> dmn_nod_vector;
   std::vector<double> dmn_nod_val_vector;
   std::vector<long> dmn_nod_vector_cond;

   GEOGetNodesInMaterialDomain(m_msh, m_st->analytical_material_group,
      dmn_nod_vector, false);
   number_of_nodes = (long) dmn_nod_vector.size();
   dmn_nod_val_vector.resize(number_of_nodes);

   for (long i = 0; i < number_of_nodes; i++)
      dmn_nod_val_vector[i] = 0;

   m_st->SetNodeValues(dmn_nod_vector, dmn_nod_vector_cond,
      dmn_nod_val_vector, ShiftInNodeVector);

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2008 JOD Implementation
 **************************************************************************/
void CSourceTermGroup::SetCOL(CSourceTerm *m_st, const int ShiftInNodeVector)
{
   long number_of_nodes;
   std::vector<long> col_nod_vector;
   std::vector<double> col_nod_val_vector;
   std::vector<long> col_nod_vector_cond;

   long i = 0;
   if (m_st->geo_name == "BOTTOM")
      i = m_msh->getNumberOfMeshLayers() - 1;

   while (i < (long) m_msh->nod_vector.size())
   {
      col_nod_vector.push_back(i);
      i += m_msh->getNumberOfMeshLayers();
   }
   number_of_nodes = (long) col_nod_vector.size();
   col_nod_val_vector.resize(number_of_nodes);

   for (long i = 0; i < number_of_nodes; i++)
      col_nod_val_vector[i] = 1;

   m_st->SetSurfaceNodeVectorConditional(col_nod_vector, col_nod_vector_cond);

   m_st->SetNodeValues(col_nod_vector, col_nod_vector_cond,
      col_nod_val_vector, ShiftInNodeVector);

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 09/2010 WW  For the case that nodes are directly given
 **************************************************************************/
void CSourceTermGroup::SetSFC(CSourceTerm* m_st, const int ShiftInNodeVector)
{
   std::vector<long> sfc_nod_vector;
   std::vector<long> sfc_nod_vector_cond;
   std::vector<double> sfc_nod_val_vector;
   Surface* m_sfc = NULL;

   m_sfc = GEOGetSFCByName(m_st->geo_name);       //CC

   if (m_sfc)
   {

      SetSurfaceNodeVector(m_sfc, sfc_nod_vector);
      if (m_st->isCoupled())
         m_st->SetSurfaceNodeVectorConditional(sfc_nod_vector,
            sfc_nod_vector_cond);
      //		m_st->SetDISType();
      SetSurfaceNodeValueVector(m_st, m_sfc, sfc_nod_vector,
         sfc_nod_val_vector);
      m_st->SetNodeValues(sfc_nod_vector, sfc_nod_vector_cond,
         sfc_nod_val_vector, ShiftInNodeVector);

   }                                              // end surface

   // 09/2010 TF geo type NODES does not exist
   //    /// If nodes are directly given. WW
   //    if(m_st->geo_type_name.find("NODES")!=std::string::npos)
   //    {
   //       string t_str = FilePath+m_st->nodes_file;
   //       ifstream ins(t_str.c_str());
   //       if(!ins.good())
   //       {
   //          cout<<" File "<<t_str<<" needed by .st is not found"<<endl;
   //          exit(1);
   //       }
   //
   //       long node_i;
   //       double val;
   //       while(!ins.eof())
   //       {
   //          stringstream ss;
   //          getline(ins, t_str);
   //          if(t_str.find("#STOP")!=std::string::npos)
   //            break;
   //
   //          ss.str(t_str);
   //          ss>> node_i >> val;
   //          sfc_nod_vector.push_back(node_i);
   //          sfc_nod_val_vector.push_back(val);
   //          ss.clear();
   //       }
   //
   //       m_st->FaceIntegration(m_msh, sfc_nod_vector, sfc_nod_val_vector);
   //
   //	   m_st->SetNodeValues(sfc_nod_vector, sfc_nod_vector_cond,
   //				sfc_nod_val_vector, ShiftInNodeVector);
   //    }

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
void CSourceTerm::SetNOD()
{

   std::vector<long> nod_vector;
   std::vector<long> nod_vector_cond;
   std::vector<double> nod_val_vector;
   int ShiftInNodeVector;

   nod_vector.push_back(msh_node_number);
   nod_vector_cond.push_back(msh_node_number);
   nod_val_vector.push_back(geo_node_value);

   /*nod_vector[0] = msh_node_number;
    nod_vector_cond[0] = msh_node_number;
    nod_val_vector[0] =geo_node_value;*/
   ShiftInNodeVector = 0;

   SetNodeValues(nod_vector, nod_vector_cond, nod_val_vector,
      ShiftInNodeVector);

}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::SetPolylineNodeVector(CGLPolyline* m_ply,
std::vector<long>&ply_nod_vector)
{
   if (m_ply->getType() == 100)                   // WW
      m_msh->GetNodesOnArc(m_ply, ply_nod_vector);
   else if (m_ply->getType() == 3)                // JOD
      m_msh->GetNODOnPLY_XY(m_ply, ply_nod_vector);
   else
      m_msh->GetNODOnPLY(m_ply, ply_nod_vector);
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::SetSurfaceNodeVector(Surface* m_sfc,
		std::vector<long>&sfc_nod_vector)
{
   m_msh->GetNODOnSFC(m_sfc, sfc_nod_vector);
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::SetPolylineNodeVectorConditional(CSourceTerm* st,
CGLPolyline* m_ply, std::vector<long>& ply_nod_vector,
std::vector<long>& ply_nod_vector_cond)
{
   size_t assembled_mesh_node, number_of_nodes;

   if (st->node_averaging)
   {
      if (m_msh_cond)
      {
         if (pcs_type_name == "RICHARDS_FLOW")
         {
            //				m_msh_cond->GetNODOnPLY(m_ply, ply_nod_vector_cond);
            m_msh_cond->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(st->getGeoObj()), ply_nod_vector_cond);
            number_of_nodes = ply_nod_vector_cond.size();
            assembled_mesh_node = ply_nod_vector[0];
            ply_nod_vector.resize(number_of_nodes);
            for (size_t i = 0; i < number_of_nodes; i++)
               ply_nod_vector[i] = assembled_mesh_node;
         }                                        // end richards
         else if (pcs_type_name == "OVERLAND_FLOW" || pcs_type_name
            == "GROUNDWATER_FLOW")                // JOD 4.10.01
         {
            number_of_nodes = ply_nod_vector.size();
            //				m_msh_cond->GetNODOnPLY(m_ply, ply_nod_vector_cond);
            m_msh_cond->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(st->getGeoObj()), ply_nod_vector_cond);
            assembled_mesh_node = ply_nod_vector_cond[0];
            ply_nod_vector_cond.resize(number_of_nodes);
            for (size_t i = 0; i < number_of_nodes; i++)
               ply_nod_vector_cond[i] = assembled_mesh_node;
         }                                        // end overland, groundwater
         else
            std::cout
               << "Warning in CSourceTermGroup::SetPolylineNodeVectorConditional - no area assembly for this process"
               << std::endl;
      }                                           // end mesh_cond
      else
         std::cout << "Warning in CSourceTermGroup::SetPLY - no MSH_COND data"
            << std::endl;
   }                                              // end area_assembly
   else
   {
      number_of_nodes = ply_nod_vector.size();
      ply_nod_vector_cond.resize(number_of_nodes);
      st->SetNOD2MSHNOD(ply_nod_vector, ply_nod_vector_cond);
   }                                              // end !area_assembly
}


// 09/2010 TF
void CSourceTermGroup::SetPolylineNodeVectorConditional(CSourceTerm* st,
std::vector<size_t>& ply_nod_vector,
std::vector<size_t>& ply_nod_vector_cond)
{
   size_t assembled_mesh_node, number_of_nodes;

   if (st->node_averaging)
   {
      if (m_msh_cond)
      {
         if (pcs_type_name == "RICHARDS_FLOW")
         {
            m_msh_cond->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(st->getGeoObj()), ply_nod_vector_cond);
            number_of_nodes = ply_nod_vector_cond.size();
            assembled_mesh_node = ply_nod_vector[0];
            ply_nod_vector.resize(number_of_nodes);
            for (size_t i = 0; i < number_of_nodes; i++)
               ply_nod_vector[i] = assembled_mesh_node;
         }                                        // end richards
         else if (pcs_type_name == "OVERLAND_FLOW"
                                                  // JOD 4.10.01
            || pcs_type_name == "GROUNDWATER_FLOW")
         {
            number_of_nodes = ply_nod_vector.size();
            m_msh_cond->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(st->getGeoObj()), ply_nod_vector_cond);
            assembled_mesh_node = ply_nod_vector_cond[0];
            ply_nod_vector_cond.resize(number_of_nodes);
            for (size_t i = 0; i < number_of_nodes; i++)
               ply_nod_vector_cond[i] = assembled_mesh_node;
         }                                        // end overland, groundwater
         else
            std::cout << "Warning in CSourceTermGroup::SetPolylineNodeVectorConditional - no area assembly for this process"
               << std::endl;
      }                                           // end mesh_cond
      else
         std::cout << "Warning in CSourceTermGroup::SetPLY - no MSH_COND data" << std::endl;
   }                                              // end area_assembly
   else
   {
      number_of_nodes = ply_nod_vector.size();
      ply_nod_vector_cond.resize(number_of_nodes);
      st->SetNOD2MSHNOD(ply_nod_vector, ply_nod_vector_cond);
   }                                              // end !area_assembly
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
//void CSourceTerm::SetPolylineRiverNodeValueVectors(CGLPolyline* m_ply,
//		int number_of_nodes) {
//
//	node_value_vectorA.resize(number_of_nodes); // HRiver
//	node_value_vectorB.resize(number_of_nodes); // KRiverBed
//	node_value_vectorC.resize(number_of_nodes); // WRiverBed
//	node_value_vectorD.resize(number_of_nodes); // TRiverBed
//	node_value_vectorE.resize(number_of_nodes); // BRiverBed
//
//
//	InterpolatePolylineRiverNodeValueVector(m_ply, DistribedBC,
//			node_value_vectorA);// HRiver
//	InterpolatePolylineRiverNodeValueVector(m_ply, DistBC_KRiverBed,
//			node_value_vectorB); // KRiverBed
//	InterpolatePolylineRiverNodeValueVector(m_ply, DistBC_WRiverBed,
//			node_value_vectorC);// WRiver
//	InterpolatePolylineRiverNodeValueVector(m_ply, DistBC_TRiverBed,
//			node_value_vectorD); // TRiverBed
//	InterpolatePolylineRiverNodeValueVector(m_ply, DistBC_BRiverBed,
//			node_value_vectorE);// BRiver
//}

//void CSourceTerm::SetPolylineRiverNodeValueVectors(const Mesh_Group::MeshNodesAlongPolyline& msh_nodes_along_polyline,
//		int number_of_nodes)
//{
//	node_value_vectorA.resize(number_of_nodes); // HRiver
//	node_value_vectorB.resize(number_of_nodes); // KRiverBed
//	node_value_vectorC.resize(number_of_nodes); // WRiverBed
//	node_value_vectorD.resize(number_of_nodes); // TRiverBed
//	node_value_vectorE.resize(number_of_nodes); // BRiverBed
//
//	InterpolatePolylineRiverNodeValueVector(msh_nodes_along_polyline, DistribedBC,
//			node_value_vectorA);// HRiver
//	InterpolatePolylineRiverNodeValueVector(msh_nodes_along_polyline, DistBC_KRiverBed,
//			node_value_vectorB); // KRiverBed
//	InterpolatePolylineRiverNodeValueVector(msh_nodes_along_polyline, DistBC_WRiverBed,
//			node_value_vectorC);// WRiver
//	InterpolatePolylineRiverNodeValueVector(msh_nodes_along_polyline, DistBC_TRiverBed,
//			node_value_vectorD); // TRiverBed
//	InterpolatePolylineRiverNodeValueVector(msh_nodes_along_polyline, DistBC_BRiverBed,
//			node_value_vectorE);// BRiver
//}

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTerm::InterpolatePolylineNodeValueVector(CGLPolyline* m_ply,
		std::vector<double>& Distribed, std::vector<double>& ply_nod_vector)
{
   for (long k = 0; k < (long) DistribedBC.size(); k++)
   {
      for (long l = 0; l < (long) m_ply->point_vector.size(); l++)
      {
         if (PointsHaveDistribedBC[k] == m_ply->point_vector[l]->id)
         {
            if (fabs(DistribedBC[k]) < MKleinsteZahl)
               DistribedBC[k] = 1.0e-20;
            m_ply->point_vector[l]->setPropert (Distribed[k]);
            break;
         }
      }
   }
   InterpolationAlongPolyline(m_ply, ply_nod_vector);
}


//void CSourceTerm::InterpolatePolylineRiverNodeValueVector(
//		const Mesh_Group::MeshNodesAlongPolyline& msh_nodes_along_polyline,
//		std::vector<double>& distribution, std::vector<double>& ply_nod_vector) const
//{
//	MATHLIB::LinearInterpolation ((msh_nodes_along_polyline.getPolyline())->getLengthVec(), distribution, msh_nodes_along_polyline.getDistOfProjNodeFromPlyStart(), ply_nod_vector);
//}

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::SetPolylineNodeValueVector(CSourceTerm* st,
CGLPolyline* ply, const std::vector<long>& ply_nod_vector,
std::vector<long>& ply_nod_vector_cond, std::vector<double>& ply_nod_val_vector)
{
   long number_of_nodes = (long) ply_nod_vector.size();
   ply_nod_val_vector.resize(number_of_nodes);

   if (st->getProcessDistributionType() == FiniteElement::LINEAR || st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
      st->InterpolatePolylineNodeValueVector(ply, st->DistribedBC,
         ply_nod_val_vector);
                                                  //System Dependented YD
   else if (st->getProcessDistributionType() == FiniteElement::SYSTEM_DEPENDENT)
   {
      CRFProcess* m_pcs = NULL;
      m_pcs = PCSGet(pcs_type_name);
      m_pcs->compute_domain_face_normal = true;   //WW
      long no_face = (long) m_msh->face_vector.size();
      for (long i = 0; i < no_face; i++)
      {
         int node_on_line = 0;
         int no_vertex = m_msh->face_vector[i]->GetVertexNumber();
         for (long jj = 0; jj < no_vertex; jj++)
         {
            for (long kk = 0; kk < number_of_nodes; kk++)
            {
               if (ply_nod_vector[kk]
                  == m_msh->face_vector[i]->GetNodeIndex(jj))
                  node_on_line++;
            }                                     // end nodes
         }                                        // end vertices
         if (node_on_line == 2)
            st->element_st_vector.push_back(
               m_msh->face_vector[i]->GetOwner()->GetIndex());
      }                                           // end faces
   }                                              // end system dependent
   else                                           //WW
   {
      for (long i = 0; i < number_of_nodes; i++)
      {
         ply_nod_val_vector[i] = st->geo_node_value;
         //			if (st->dis_type == 12)
         if (st->getProcessDistributionType() == FiniteElement::CONSTANT_GEO)
            ply_nod_val_vector[i] = st->geo_node_value
               / (double) number_of_nodes;        // distribute flow to nodes along polyline. To do.. 4.10.06
      }
   }
                                                  /*|| st->getProcessDistributionType() == PHILIP */
   if (st->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN || st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN
      || st->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
   {
      if (m_msh->GetMaxElementDim() == 1)         // 1D  //WW MB
         st->DomainIntegration(m_msh, ply_nod_vector, ply_nod_val_vector);
      else
         st->EdgeIntegration(m_msh, ply_nod_vector, ply_nod_val_vector);
   }

   //	if (st->getProcessDistributionType() == RIVER) // river
   //		st->SetPolylineRiverNodeValueVectors(ply, number_of_nodes);

   //	if (st->dis_type == 5 || st->dis_type == 6 || st->dis_type == 8
   //				|| st->dis_type == 9) {
   if (                                           //st->getProcessDistributionType()  == RIVER
      /*||*/ st->getProcessDistributionType()  == FiniteElement::CRITICALDEPTH
      || st->getProcessDistributionType()  == FiniteElement::NORMALDEPTH
      || st->getProcessDistributionType()  == FiniteElement::ANALYTICAL)
   {
      st->node_value_vectorArea.resize(number_of_nodes);
      for (long i = 0; i < number_of_nodes; i++)
         st->node_value_vectorArea[i] = 1.0;      //Element width !
      st->EdgeIntegration(m_msh, ply_nod_vector,
         st->node_value_vectorArea);
   }

   if (st->isCoupled() && st->node_averaging)
      AreaAssembly(st, ply_nod_vector_cond, ply_nod_val_vector);
}


// 09/2010 TF
void CSourceTermGroup::SetPolylineNodeValueVector(CSourceTerm* st,
const Mesh_Group::MeshNodesAlongPolyline& msh_nodes_along_polyline,
const std::vector<long>& ply_nod_vector_cond,
std::vector<double>& ply_nod_val_vector) const
{
   const std::vector<size_t> ply_nod_vector(msh_nodes_along_polyline.getNodeIDs());
   size_t number_of_nodes(ply_nod_vector.size());
   ply_nod_val_vector.resize(number_of_nodes);

   // horrible !!! since template parameter difference in std::vector (size_t vs. long)
   std::vector<long> ply_nod_vector_long;
   for (size_t k(0); k < ply_nod_vector.size(); k++)
      ply_nod_vector_long.push_back(ply_nod_vector[k]);

   FiniteElement::DistributionType distype (st->getProcessDistributionType());

                                                  // linear
   if (distype == FiniteElement::LINEAR|| distype == FiniteElement::LINEAR_NEUMANN)
   {
      //		st->InterpolatePolylineRiverNodeValueVector(msh_nodes_along_polyline,
      //				st->DistribedBC, ply_nod_val_vector);
      MATHLIB::LinearInterpolation (
         (msh_nodes_along_polyline.getPolyline())->getLengthVec(),
         st->getDistribution(),
         msh_nodes_along_polyline.getDistOfProjNodeFromPlyStart(),
         ply_nod_val_vector);
   }                                              //System Dependented YD
   else if (distype == FiniteElement::SYSTEM_DEPENDENT)
   {
      CRFProcess* m_pcs(PCSGet(pcs_type_name));
      m_pcs->compute_domain_face_normal = true;   //WW
      long no_face = (long) m_msh->face_vector.size();
      for (long i = 0; i < no_face; i++)
      {
         int node_on_line = 0;
         int no_vertex = m_msh->face_vector[i]->GetVertexNumber();
         for (long jj = 0; jj < no_vertex; jj++)
         {
            for (size_t kk = 0; kk < number_of_nodes; kk++)
            {
               if (ply_nod_vector[kk]
                  == static_cast<size_t> (m_msh->face_vector[i]->GetNodeIndex(
                  jj)))
                  node_on_line++;
            }                                     // end nodes
         }                                        // end vertices
         if (node_on_line == 2)
            st->element_st_vector.push_back(
               m_msh->face_vector[i]->GetOwner()->GetIndex());
      }                                           // end faces
   }                                              // end system dependent
   else                                           //WW
   {
      for (size_t i = 0; i < number_of_nodes; i++)
      {
         ply_nod_val_vector[i] = st->geo_node_value;
      }
   }
   if (distype == FiniteElement::CONSTANT_NEUMANN || distype == FiniteElement::LINEAR_NEUMANN || distype == FiniteElement::GREEN_AMPT)
   {
      if (m_msh->GetMaxElementDim() == 1)         // 1D  //WW MB
         st->DomainIntegration(m_msh, ply_nod_vector_long,
            ply_nod_val_vector);
      else
         st->EdgeIntegration(m_msh, ply_nod_vector_long, ply_nod_val_vector);
   }

   //	if (distype == RIVER) // river
   //		st->SetPolylineRiverNodeValueVectors(msh_nodes_along_polyline, number_of_nodes);

   if (/*distype == RIVER ||*/ distype == FiniteElement::CRITICALDEPTH || distype == FiniteElement::NORMALDEPTH || distype == FiniteElement::ANALYTICAL) {
   st->node_value_vectorArea.resize(number_of_nodes);
   for (size_t i = 0; i < number_of_nodes; i++)
      st->node_value_vectorArea[i] = 1.0;         //Element width !
   st->EdgeIntegration(m_msh, ply_nod_vector_long,
      st->node_value_vectorArea);
}


if (st->isCoupled() && st->node_averaging)
   AreaAssembly(st, ply_nod_vector_cond, ply_nod_val_vector);
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::AreaAssembly(const CSourceTerm* const st,
const std::vector<long>& ply_nod_vector_cond,
std::vector<double>& ply_nod_val_vector) const
{
   if (pcs_type_name == "RICHARDS_FLOW")
   {
      if (m_msh_cond->GetMaxElementDim() == 1)    // 1D  //WW MB
         st->DomainIntegration(m_msh_cond, ply_nod_vector_cond,
            ply_nod_val_vector);
      else
         st->EdgeIntegration(m_msh_cond, ply_nod_vector_cond,
            ply_nod_val_vector);
      double sum_node_value = 0;
      for (size_t i = 0; i < ply_nod_val_vector.size(); i++)
         sum_node_value += ply_nod_val_vector[i];
      for (size_t  i = 0; i < ply_nod_val_vector.size(); i++)
         ply_nod_val_vector[i] /= sum_node_value;
   }
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::SetSurfaceNodeValueVector(CSourceTerm* st,
		Surface* m_sfc, std::vector<long>&sfc_nod_vector,
		std::vector<double>&sfc_nod_val_vector)
{
   // CRFProcess* m_pcs = NULL;
   // m_pcs = PCSGet(pcs_type_name);
   long number_of_nodes = (long) sfc_nod_vector.size();
   sfc_nod_val_vector.resize(number_of_nodes);

   for (long i = 0; i < number_of_nodes; i++)
      sfc_nod_val_vector[i] = st->geo_node_value;
   // KR & TF - case not used
   //	if (m_st->dis_type == 12) //To do. 4.10.06
   //		for (long i = 0; i < number_of_nodes; i++)
   //			sfc_nod_val_vector[i] = m_st->geo_node_value
   //					/ (double) number_of_nodes;

   //	if (st->dis_type == 2 || st->dis_type == 4) { // Piecewise linear distributed, polygon-wise WW
   if (st->getProcessDistributionType() == FiniteElement::LINEAR || st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
   {
      CGLPolyline* m_ply = NULL;
      std::vector<CGLPolyline*>::iterator p =
         m_sfc->polyline_of_surface_vector.begin();
      p = m_sfc->polyline_of_surface_vector.begin();
      while (p != m_sfc->polyline_of_surface_vector.end())
      {
         m_ply = *p;
         for (long k = 0; k < (long) st->DistribedBC.size(); k++)
         {
            for (long l = 0; l < (long) m_ply->point_vector.size(); l++)
            {
               if (st->PointsHaveDistribedBC[k]
                  == m_ply->point_vector[l]->id)
               {
                  if (fabs(st->DistribedBC[k]) < MKleinsteZahl)
                     st->DistribedBC[k] = 1.0e-20;
                  m_ply->point_vector[l]->setPropert (st->DistribedBC[k]);
                  break;
               }                                  // end l
            }                                     // end k
         }                                        // end polyline
         // InterpolationAlongPolyline(m_polyline, node_value_vector);
         p++;
      }                                           // end while
   }                                              // end linear

   // neumann, Green-Ampt, Philip
   //	if (st->dis_type == 3 || st->dis_type == 4 || st->dis_type == 10
   //				|| st->dis_type == 11) {
                                                  /*|| st->getProcessDistributionType() == PHILIP */
   if (st->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN || st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN
      || st->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
   {
      if (m_msh->GetMaxElementDim() == 2)         // For all meshes with 1-D or 2-D elements
         st->DomainIntegration(m_msh, sfc_nod_vector, sfc_nod_val_vector);
      else if (m_msh->GetMaxElementDim() == 3)    // For all meshes with 3-D elements
         st->FaceIntegration(m_msh, sfc_nod_vector, sfc_nod_val_vector);
   }                                              // end neumann

}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTerm::SetSurfaceNodeVectorConditional(std::vector<long>&sfc_nod_vector,
std::vector<long>&sfc_nod_vector_cond)
{
   long number_of_nodes;
   number_of_nodes = (long) sfc_nod_vector.size();

   sfc_nod_vector_cond.resize(number_of_nodes);
   SetNOD2MSHNOD(sfc_nod_vector, sfc_nod_vector_cond);
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTerm::SetNodeValues(const std::vector<long>& nodes, const std::vector<long>& nodes_cond,
const std::vector<double>&node_values, int ShiftInNodeVector)
{
   CNodeValue *m_nod_val = NULL;
   long number_of_nodes = (long) nodes.size();

   for (long i = 0; i < number_of_nodes; i++)
   {
      m_nod_val = new CNodeValue();
      m_nod_val->msh_node_number = nodes[i] + ShiftInNodeVector;
      m_nod_val->geo_node_number = nodes[i];
      m_nod_val->setProcessDistributionType (getProcessDistributionType());
      m_nod_val->node_value = node_values[i];
      m_nod_val->CurveIndex = CurveIndex;
      if (_coupled)                               // JOD 4.7.10
      {
         m_nod_val->msh_node_number_conditional = nodes_cond[i];
                                                  // JOD 4.10.01
         if ((getProcessType() == OVERLAND_FLOW || getProcessType() == GROUNDWATER_FLOW) && node_averaging)
         {
            double weights = 0;
            for (long j = 0; j < number_of_nodes; j++)
            {
               m_nod_val->msh_node_numbers_averaging.push_back(nodes[j]);
               m_nod_val->msh_node_weights_averaging.push_back(
                  node_values[j]);
               weights += node_values[j];
            }
            for (long j = 0; j < number_of_nodes; j++)
               m_nod_val->msh_node_weights_averaging[j] /= weights;
         }
      }
      //WW        group_vector.push_back(m_node_value);
      //WW        st_group_vector.push_back(m_st); //OK
      //		if (getProcessDistributionType() == RIVER) {
      //			m_nod_val->node_value = node_value_vectorArea[i];
      //			m_nod_val->node_parameterA = node_value_vectorA[i];
      //			m_nod_val->node_parameterB = node_value_vectorB[i];
      //			m_nod_val->node_parameterC = node_value_vectorC[i];
      //			m_nod_val->node_parameterD = node_value_vectorD[i];
      //			m_nod_val->node_parameterE = node_value_vectorE[i];
      //		}
      //		if (dis_type == 6 || dis_type == 8 || dis_type == 9) // critical depth, normal depth, analytical
      if (getProcessDistributionType() == FiniteElement::CRITICALDEPTH
         || getProcessDistributionType() == FiniteElement::NORMALDEPTH
         || getProcessDistributionType() == FiniteElement::ANALYTICAL)
      {
         m_nod_val->node_value = node_value_vectorArea[i];
                                                  //CMCD bugfix on 4.9.06
         m_nod_val->node_area = node_value_vectorArea[i];
      }
      _pcs->st_node_value.push_back(m_nod_val);   //WW
      _pcs->st_node.push_back(this);              //WW
   }                                              // end nodes

}


// 09/2010 TF
void CSourceTerm::SetNodeValues(const std::vector<size_t>& nodes, const std::vector<size_t>& nodes_cond,
		const std::vector<double>& node_values, int ShiftInNodeVector)
{
   size_t number_of_nodes (nodes.size());

   for (size_t i = 0; i < number_of_nodes; i++)
   {
      CNodeValue *m_nod_val = new CNodeValue();
      m_nod_val->msh_node_number = nodes[i] + ShiftInNodeVector;
      m_nod_val->geo_node_number = nodes[i];
      m_nod_val->setProcessDistributionType (getProcessDistributionType());
      m_nod_val->node_value = node_values[i];
      m_nod_val->CurveIndex = CurveIndex;

      if (_coupled)                               // JOD 4.7.10
      {
         m_nod_val->msh_node_number_conditional = nodes_cond[i];
         if ((getProcessType() == OVERLAND_FLOW
            || getProcessType() == GROUNDWATER_FLOW)
            && node_averaging)                    // JOD 4.10.01
         {
            double weights = 0;
            for (size_t j = 0; j < number_of_nodes; j++)
            {
               m_nod_val->msh_node_numbers_averaging.push_back(nodes[j]);
               m_nod_val->msh_node_weights_averaging.push_back(
                  node_values[j]);
               weights += node_values[j];
            }
            for (size_t j = 0; j < number_of_nodes; j++)
               m_nod_val->msh_node_weights_averaging[j] /= weights;
         }
      }

      //		if (getProcessDistributionType() == RIVER) {
      //			m_nod_val->node_value = node_value_vectorArea[i];
      //			m_nod_val->node_parameterA = node_value_vectorA[i];
      //			m_nod_val->node_parameterB = node_value_vectorB[i];
      //			m_nod_val->node_parameterC = node_value_vectorC[i];
      //			m_nod_val->node_parameterD = node_value_vectorD[i];
      //			m_nod_val->node_parameterE = node_value_vectorE[i];
      //		}
      if (getProcessDistributionType() == FiniteElement::CRITICALDEPTH
         || getProcessDistributionType() == FiniteElement::NORMALDEPTH
         || getProcessDistributionType() == FiniteElement::ANALYTICAL)
      {
         m_nod_val->node_value = node_value_vectorArea[i];
                                                  //CMCD bugfix on 4.9.06
         m_nod_val->node_area = node_value_vectorArea[i];
      }
      _pcs->st_node_value.push_back(m_nod_val);   //WW
      _pcs->st_node.push_back(this);              //WW
   }                                              // end nodes
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2005 MB
 last modification:
 **************************************************************************/
//void CSourceTerm::SetNOD2MSHNOD(vector<long>&nodes,
//		vector<long>&conditional_nodes)
//{
//	CGLPoint* m_pnt = NULL;
//	long number;
//	CFEMesh* m_msh_cond = NULL;
//	CFEMesh* m_msh_this = NULL;
//
//	m_msh_cond = FEMGet(pcs_type_name_cond);
//	m_msh_this = FEMGet(convertProcessTypeToString(m_st->getProcessType()));
//	m_pnt = new CGLPoint;
//
//	for (long i = 0; i < (long) nodes.size(); i++) {
//		m_pnt->x = m_msh_this->nod_vector[nodes[i]]->X();
//		m_pnt->y = m_msh_this->nod_vector[nodes[i]]->Y();
//		m_pnt->z = m_msh_this->nod_vector[nodes[i]]->Z();
//
//		number = m_msh_cond->GetNODOnPNT(m_pnt);
//		conditional_nodes[i] = number;
//
//	}
//
//	delete m_pnt;
//
//}

void CSourceTerm::SetNOD2MSHNOD(std::vector<long>& nodes,
std::vector<long>& conditional_nodes)
{
   CFEMesh* m_msh_cond(FEMGet(pcs_type_name_cond));
   CFEMesh* m_msh_this(FEMGet(convertProcessTypeToString(getProcessType())));

   GEOLIB::Point pnt;
   for (size_t i = 0; i < nodes.size(); i++)
   {
      pnt[0] = m_msh_this->nod_vector[nodes[i]]->X();
      pnt[1] = m_msh_this->nod_vector[nodes[i]]->Y();
      pnt[2] = m_msh_this->nod_vector[nodes[i]]->Z();

      conditional_nodes[i] = m_msh_cond->GetNODOnPNT(&pnt);
   }
}


void CSourceTerm::SetNOD2MSHNOD(const std::vector<size_t>& nodes,
std::vector<size_t>& conditional_nodes) const
{
   CFEMesh* m_msh_cond(FEMGet(pcs_type_name_cond));
   CFEMesh* m_msh_this(FEMGet(convertProcessTypeToString(getProcessType())));

   GEOLIB::Point pnt;
   for (size_t i = 0; i < nodes.size(); i++)
   {
      pnt[0] = m_msh_this->nod_vector[nodes[i]]->X();
      pnt[1] = m_msh_this->nod_vector[nodes[i]]->Y();
      pnt[2] = m_msh_this->nod_vector[nodes[i]]->Z();

      conditional_nodes[i] = m_msh_cond->GetNODOnPNT(&pnt);
   }
}


/**************************************************************************
 GeoSys source term function:
 02/2009 WW Implementation
 **************************************************************************/
inline void CSourceTerm::DirectAssign(long ShiftInNodeVector)
{
   std::string line_string;
   std::string st_file_name;
   std::stringstream in;
   long n_index;
   double n_val;
   CRFProcess* m_pcs = NULL;
   CNodeValue *m_nod_val = NULL;
   m_pcs = PCSGet(convertProcessTypeToString(getProcessType()));

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
      m_nod_val = new CNodeValue();
      m_nod_val->msh_node_number = n_index + ShiftInNodeVector;
      m_nod_val->geo_node_number = n_index;
      m_nod_val->setProcessDistributionType (getProcessDistributionType());
      m_nod_val->node_value = n_val;
      m_nod_val->CurveIndex = CurveIndex;
      m_pcs->st_node_value.push_back(m_nod_val);
      m_pcs->st_node.push_back(this);
      //
   }                                              // eof
}


/**************************************************************************
GeoSys source term function:
03/2010 WW Implementation
To read GIS raster files as source term.
**************************************************************************/
std::string CSourceTerm::DirectAssign_Precipitation(double current_time)
{
   int i, size;
   double stepA, stepB, tim_val;

   long l, nbc_node, n_index, osize = 0;
   double n_val;

   std::string fileA, fileB;

   CRFProcess* m_pcs = NULL;
   CNodeValue *m_nod_val = NULL;
                                                  //PCSGet(pcs_type_name);
   m_pcs = PCSGet(convertProcessTypeToString (getProcessType()));

   if(start_pos_in_st<0)
      osize = (long)m_pcs->st_node.size();

   size = (int)precip_times.size();
   stepA = 0.;
   stepB = 0.;
   if(current_time < m_pcs->Tim->time_start||fabs(current_time - m_pcs->Tim->time_start)<DBL_MIN)
   {
      fileA = precip_files[0];
      stepB = -1.;
   }
   else if(current_time > precip_times[size-1]||fabs(current_time - precip_times[size-1])<DBL_MIN)
   {
      fileA = precip_files[size-1];
      stepB = -1.;
   }
   else
   {
      double step_b = DBL_MAX;
      double step_f = DBL_MAX;

      for(i=0; i<size; i++)
      {
         tim_val = precip_times[i];
         if(current_time>tim_val)
         {
            if((current_time-tim_val)<step_b)
            {
               step_b = current_time-tim_val;
               stepA = tim_val;
               fileA = precip_files[i];
            }
         }
         else
         {

            if((tim_val-current_time)<step_f)
            {
               step_f = tim_val-current_time;
               stepB = tim_val;
               fileB = precip_files[i];
            }
         }
      }
      if(fabs(stepA-current_time)<DBL_MIN)
         stepB = -1.;
      if(fabs(current_time-stepB)<DBL_MIN)
      {
         fileA = fileB;
         stepB = -1.;
      }
      if(fabs(stepA-stepB)<DBL_MIN)
         stepB = -1.;
   }

   fileA = FilePath + fileA;
   std::ifstream bin_sA(fileA.c_str(), std::ios::binary);
   std::ifstream bin_sB;

   if(!bin_sA.good())
   {
      std::cout<<"Could not find file "<< fileA<<std::endl;
      exit(0);
   }
   else
   {
      std::cout<<"Successfully read file: "<< fileA<<std::endl; //written by Miao Jing. 03.2017
   }

   bin_sA.setf(std::ios::scientific,std::ios::floatfield);
   bin_sA.precision(14);

   /*
   if(stepB>0.)
   {
      fileB = FilePath + fileB;
      bin_sB.open(fileB.c_str(), ios::binary);
      if(!bin_sB.good())
      {
          cout<<"Could not find file "<< fileB<<endl;
          exit(0);
      }
      bin_sB.setf(ios::scientific,ios::floatfield);
   bin_sB.precision(14);
   bin_sB.read((char*)(&nbc_node), sizeof(nbc_node));
   }
   */

   bin_sA.read((char*)(&nbc_node), sizeof(nbc_node));

   double valA;                                   //, valB = 0.;
   for(l=0; l<nbc_node; l++)
   {

      bin_sA.read((char*)(&n_index), sizeof(n_index));
      bin_sA.read((char*)(&valA), sizeof(valA));
      /*
      if( stepB>0.)
      {
         bin_sB.read((char*)(&n_index), sizeof(n_index));
         bin_sB.read((char*)(&valB), sizeof(valB));
         n_val = valA + (current_time - stepA)*(valB-valA)/(stepB-stepA);
      }
      else
      */
      n_val = valA;

      //
      if(start_pos_in_st<0)
      {
         m_nod_val = new CNodeValue();
         m_pcs->st_node_value.push_back(m_nod_val);
         m_pcs->st_node.push_back(this);
      }
      else
         m_nod_val = m_pcs->st_node_value[l+start_pos_in_st];

      m_nod_val->msh_node_number = n_index ;
      m_nod_val->geo_node_number = n_index;
                                                  //node_distype = dis_type;
      m_nod_val->setProcessDistributionType (getProcessDistributionType());
      m_nod_val->node_value = n_val;
      m_nod_val->CurveIndex = CurveIndex;


	//  std::cout << "New node " << l << std::endl;
      //
   }                                              //

   if(start_pos_in_st<0)
      start_pos_in_st = osize;

   bin_sA.close();
   //bin_sB.close();

   return fileA;
}


/**************************************************************************
 FEMLib-Method:
 Task: Analytical diffusion in matrix. Replaces matrix. See paper to be issued.
 Programing:
 11/2005 CMCD Implementation
 04/2006 Moved from CSourceTermGroup and changed the arguments
 last modification:
 04/2006 CMCD Updated
 **************************************************************************/
                                                  // , CSourceTerm *m_st)
double CSourceTerm::GetAnalyticalSolution(long location)
{
   int idx, n;
   int size, process_no;
   long i;
   long step, no_terms_included;
   double value, source, gradient, ref_value = 0.0;
   double timevalue;
   double fac = 1.0;
   double temp_time, temp_value;
   double pi = 3.1415926;
   double D = this->analytical_diffusion;
   double ne = this->analytical_porosity;
   double tort = this->analytical_tortousity;
   double Kd = this->analytical_linear_sorption_Kd;
   double rho = this->analytical_matrix_density;
   double Dtrans = (D * ne) / ((ne + Kd * rho) * tort);
   double Dsteady = D * ne / tort;
   double t0, tn, tnn, val1, val2, node_area;
   double tvol, vol, flux_area, tflux_area;
   double mass_solute_present, mass_to_remove;
   //WW  bool out = false;
   //WW  int dimension = this->analytical_material_group;
   std::string process;
   CRFProcess* m_pcs = NULL;
   m_pcs = PCSGet(convertProcessTypeToString(this->getProcessType()), convertPrimaryVariableToString(this->getProcessPrimaryVariable()));
   CFEMesh* m_msh = m_pcs->m_msh;                 //WW
   CElem* Ele = NULL;
   long node_number = location;                   //WW m_pcs->st_node_value[location]->msh_node_number;
   CNode* Node = m_msh->nod_vector[node_number];
   double area = m_pcs->st_node_value[location]->node_area;
   std::vector<double> time_history;
   std::vector<double> value_history;
   //Initialise
   time_history.clear();
   value_history.clear();
   t0 = tn = tnn = source = gradient = val1 = val2 = node_area = flux_area
      = 0.0;
   idx = m_pcs->GetNodeValueIndex(convertPrimaryVariableToString(this->getProcessPrimaryVariable()));
   value = m_pcs->GetNodeValue(node_number, idx);
   if (value < MKleinsteZahl)
      value = 0.0;
   timevalue = aktuelle_zeit;
   step = aktueller_zeitschritt;
   if (step < 10)
      step = 10;
   size = (int) analytical_processes.size();
   process_no = 0;
   //Domain or Polyline
   for (i = 0; i < size; i++)
   {
      if (analytical_processes[i] == convertPrimaryVariableToString(this->getProcessPrimaryVariable()))
      {
         if (this->getGeoType () == GEOLIB::POLYLINE)
         {
            if (this->getGeoName().compare(
               analytical_processes_polylines[i]) == 0)
               process_no = i;
         }
         //			if (this->geo_type_name.compare("DOMAIN") == 0)
         if (this->getGeoType () == GEOLIB::GEODOMAIN)
            process_no = i;
      }
   }
   //Identify process
   if (process_no == 1)
   {
      process_no = process_no;
   }
   process_no *= 2;                               //first column time, second column value, hence two columns per process;

   //If time step require new calculation of source term then start
   if ((aktueller_zeitschritt - 1) % this->resolution == 0)
   {

      //Save data in a vector attached to the nodes
      this->SetNodePastValue(node_number, process_no, 0, timevalue);
      this->SetNodePastValue(node_number, process_no + 1, 0, value);

      //Recall historical data
      ref_value = this->GetNodePastValueReference(node_number, (process_no
         / 2));
      if ((size_t)step > this->number_of_terms)
         no_terms_included = this->number_of_terms;
      else
         no_terms_included = step;
      for (i = 0; i < no_terms_included; i++)
      {
         temp_time = this->GetNodePastValue(node_number, process_no, i);
         temp_value
            = (this->GetNodePastValue(node_number, process_no + 1, i)
            - ref_value);
         time_history.push_back(temp_time);
         value_history.push_back(temp_value);
      }

      //Calculate individual terms and sum for gradient
      for (i = no_terms_included - 1; i > 0; i--)
      {
         t0 = time_history[0];
         if (i == no_terms_included - 1)
            tn = (t0 - time_history[i]) + (time_history[i - 1]
               - time_history[i]);
         else
            tn = t0 - time_history[i + 1];
         tnn = t0 - time_history[i];
         val1 = 1 / (sqrt(pi * Dtrans * tn));
         val2 = 1 / (sqrt(pi * Dtrans * tnn));
         gradient += ((val2 - val1) * value_history[i]);
      }
      tn = t0 - time_history[1];
      tnn = 0;
      val1 = 1 / (sqrt(pi * Dtrans * tn));
      gradient -= (val1 * value_history[0]);

      //Area calculations
      mass_solute_present = 1.e99;                //Initially set value very high

      //Area for lines, triangles and quads in domain.
      //  if (area < 0) {//set in set source terms function, domain area = -1 to start with
      if (area < DBL_MIN)                         // HS 04.2008
      {
         tvol = 0.0;
         tflux_area = 0.0;
         for (i = 0; i < (int) Node->connected_elements.size(); i++)
         {
            Ele = m_msh->ele_vector[Node->connected_elements[i]];
            vol = Ele->GetVolume();               //Assuming 1m thickness
            flux_area = Ele->GetFluxArea();       //Real thickness for a fracture
            n = Ele->GetVertexNumber();
            tvol += (vol / n);
            tflux_area += (flux_area / n);
         }
         node_area = tvol * 2.;                   //* 2 because the diffusion is in two direction perpendicular to the fracture
         mass_solute_present = tflux_area * tvol * value;
      }
      //Area for polylines
      else
         node_area = area;

      //factor for conversion to energy for temperature if necessary
      fac = this->factor;
      source = gradient * Dsteady * fac * node_area;
      mass_to_remove = fabs(source) * dt;
      if (mass_to_remove > mass_solute_present)
      {
         source *= (mass_solute_present / mass_to_remove);
      }
      this->SetNodeLastValue(node_number, (process_no / 2), source);

   }                                              // If new source term calculation not required
   else
      source = this->GetNodeLastValue(node_number, (process_no / 2));
   return source;
}


/**************************************************************************
 FEMLib-Method:
 Task: master write function
 Programing:
 11/2005 CMCD Implementation, functions to access and write the time history data
 of nodes
 last modification:
 **************************************************************************/
void CSourceTerm::SetNodePastValue(long n, int idx, int pos, double value)
{
   bool endstepreached = false;
   pos = pos;                                     //WW

   //Check whether this is the first call
   CRFProcess* m_pcs = NULL;
   m_pcs = PCSGet(convertProcessTypeToString(getProcessType()), convertPrimaryVariableToString (getProcessPrimaryVariable()));
   if (!m_pcs)                                    //OK
   {
      std::cout << "Warning in SetNodePastValue - no PCS data" << std::endl;
      return;
   }

   size_t size1 (m_pcs->nod_val_vector.size());

   //Create memory for the values
   if (node_history_vector.empty())
   {
      //WW     int number_of_terms = max_no_terms;
      for (size_t k = 0; k < size1; k++)
      {
         NODE_HISTORY *nh = new NODE_HISTORY;
         CreateHistoryNodeMemory(nh);
         node_history_vector.push_back(nh);
      }
      for (size_t k = 0; k < size1; k++)
      {
         for (size_t j = 0; j < _no_an_sol; j++)
            node_history_vector[k]->value_reference.push_back(-1.0);
      }
   }

   //Store the first set of values as reference values
   int flipflop = idx % 2;
   if (aktueller_zeitschritt == 1)
   {
      //if (size2 == idx)
      if (flipflop == 1)
         node_history_vector[n]->value_reference[(idx - 1) / 2] = value;
   }
   //  size2 = (int) node_history_vector[n]->value_reference.size();
   size_t steps = 10;
   if (aktueller_zeitschritt > steps)
      steps = aktueller_zeitschritt;
   if (_max_no_terms >= steps)
      steps = aktueller_zeitschritt;
   else
   {
      steps = _max_no_terms;
      endstepreached = true;
   }
   //Enter the value and push the other values back
   if (!endstepreached)
   {
      for (size_t k = steps - 1; k > 0; k--)
         node_history_vector[n]->value_store[idx][k]
            = node_history_vector[n]->value_store[idx][k - 1];
      node_history_vector[n]->value_store[idx][0] = value;
   }
   double cutvalue = 0.0;
   double nextvalue = 0.0;
   long no_steps_past_cutof = 0;
   if (endstepreached)
   {
      no_steps_past_cutof = aktueller_zeitschritt - _max_no_terms;
      cutvalue = node_history_vector[n]->value_store[idx][steps - 1];
      nextvalue = node_history_vector[n]->value_store[idx][steps - 2];
      node_history_vector[n]->value_store[idx][steps - 1] = (cutvalue
         * (double) (no_steps_past_cutof) + nextvalue)
         / ((double) no_steps_past_cutof + 1);
      for (size_t k = steps - 2; k > 0; k--)
         node_history_vector[n]->value_store[idx][k]
            = node_history_vector[n]->value_store[idx][k - 1];
      node_history_vector[n]->value_store[idx][0] = value;
   }
}


void CSourceTerm::SetNodeLastValue(long n, int idx, double value)
{
   size_t size (node_history_vector[n]->last_source_value.size());
   if (size == 0)
   {
      for (size_t i = 0; i < _no_an_sol; i++)
         node_history_vector[n]->last_source_value.push_back(0);
   }
   node_history_vector[n]->last_source_value[idx] = value;
}


double CSourceTerm::GetNodeLastValue(long n, int idx)
{
   double value = 0.0;
   value = node_history_vector[n]->last_source_value[idx];
   return value;
}


double CSourceTerm::GetNodePastValue(long n, int idx, int pos)
{
   double value;
   value = node_history_vector[n]->value_store[idx][pos];
   return value;
}


double CSourceTerm::GetNodePastValueReference(long n, int idx)
{
   double value;
   value = node_history_vector[n]->value_reference[idx];
   return value;
}


void CSourceTerm::CreateHistoryNodeMemory(NODE_HISTORY* nh)
{
   size_t s_col = _no_an_sol * 2;
   size_t s_row = number_of_terms;

   nh->value_store = new double*[s_col];
   for (size_t i = 0; i < s_col; i++)
      nh->value_store[i] = new double[s_row];

   for (size_t k = 0; k < s_col; k++)
   {
      for (size_t j = 0; j < s_row; j++)
         nh->value_store[k][j] = 0.0;
   }
}


void CSourceTerm::DeleteHistoryNodeMemory()
{
   size_t size (node_history_vector.size());
   size_t s_row = _no_an_sol * 2;

   if (size > 0)
   {
      for (size_t j = 0; j < size; j++)
      {
         for (size_t i = 0; i < s_row; i++)
            delete node_history_vector[j]->value_store[i];
         delete node_history_vector[j]->value_store;
      }
      node_history_vector.clear();
   }
}


///**************************************************************************
//FEMLib-Method:
//07/2007 OK Implementation
//modifications:
//05/2010 TF restructured a little bit
//**************************************************************************/
//CSourceTerm* STGet(const std::string& pcs_name, const std::string& geo_type_name, const std::string& geo_name)
//{
//  for(size_t i=0;i<st_vector.size();i++) {
//    if((st_vector[i]->pcs_type_name.compare(pcs_name)==0) &&
//       (st_vector[i]->geo_type_name.compare(geo_type_name)==0) &&
//       (st_vector[i]->getGeoName().compare(geo_name)==0))
//      return st_vector[i];
//
////    if((st_vector[i]->pcs_type_name.compare(pcs_name)==0) &&
////    		(st_vector[i]->geo_type_name.compare(geo_type_name)==0) &&
////    		geo_type_name.compare ("POINT") == 0 &&
////    		(st_vector[i]->getGeoObjIdx() == compare(geo_name)==0))
////    	return st_vector[i];
//  }
//  return NULL;
//}

/**************************************************************************
 FEMLib-Method: 4.7.10 shift and average field variables
 10/2008 JOD Implementation
 **************************************************************************/
void GetCouplingFieldVariables(double* h_this, double* h_cond,
double* h_shifted, double* h_averaged, double* z_this, double* z_cond,
CSourceTerm* m_st, CNodeValue* cnodev)
{

   int nidx, nidx_cond;
   double gamma;
   CRFProcess* m_pcs_this = NULL;
   CRFProcess* m_pcs_cond = NULL;

   m_pcs_this = PCSGet(convertProcessTypeToString(m_st->getProcessType()));
   m_pcs_cond = PCSGet(m_st->pcs_type_name_cond);
                                                  // only one phase
   gamma = mfp_vector[0]->Density() * GRAVITY_CONSTANT;
   *z_this = m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->Z();
   *z_cond
      = m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->Z();
   nidx = m_pcs_this->GetNodeValueIndex (convertPrimaryVariableToString (m_st->getProcessPrimaryVariable())) + 1;
   nidx_cond = m_pcs_cond->GetNodeValueIndex(m_st->pcs_pv_name_cond) + 1;
   *h_this = m_pcs_this->GetNodeValue(cnodev->msh_node_number, nidx);
   *h_cond = m_pcs_cond->GetNodeValue(cnodev->msh_node_number_conditional,
      nidx_cond);

   if (m_st->getProcessType() == OVERLAND_FLOW)
   {
      if (m_st->node_averaging)                   // shift overland node on soil column top, averaging over nodes
      {
         *h_shifted = *h_this - *z_this + *z_cond;
         *h_averaged = 0;
         for (long i = 0; i
            < (long) cnodev->msh_node_numbers_averaging.size(); i++)
            *h_averaged
               += cnodev->msh_node_weights_averaging[i]
               * (m_pcs_this->GetNodeValue(
               cnodev->msh_node_numbers_averaging[i],
               nidx)
               - m_pcs_this->m_msh->nod_vector[cnodev->msh_node_numbers_averaging[i]]->Z());

         *h_averaged += *z_cond;
         *z_this = *z_cond;
      }                                           // end averaging
      else                                        // no averaging
      {
         *h_shifted = *h_this;
         *h_averaged = *h_this;
      }                                           // end no averaging

      if (m_st->pcs_pv_name_cond == "PRESSURE1")
      {
         *h_cond /= gamma;
         *h_cond += *z_cond;
      }
      if (m_st->pcs_type_name_cond == "GROUNDWATER_FLOW")
         h_cond = std::max(h_cond, z_this);       //groundwater level might not reach overland flow bottom, than only hydrostatic overland pressure
   }                                              // end overland flow
   else                                           // richards & groundwater flow
   {
      if (m_st->pcs_pv_name_cond == "PRESSURE1")  // JOD 4.10.01
      {
         *h_cond /= gamma;
         *h_cond += *z_cond;
      }
      if (m_st->node_averaging)                   // shift overland/groundwater node on soil column top, averaging over nodes
      {
         *h_shifted = *h_cond - *z_cond;
         *h_shifted += *z_this;
         *z_cond = *z_this;
      }                                           // end averaging
      else
         *h_shifted = *h_cond;

      if (m_st->getProcessPrimaryVariable() == PRESSURE)
      {
         *h_this /= gamma;
         *h_this += *z_this;
      }
   }                                              // end richards & groundwater flow

}


/**************************************************************************
 FEMLib-Method: 4.7.10
 10/2008 JOD Implementation
 **************************************************************************/
double CalcCouplingValue(double factor, double h_this, double h_cond,
double z_cond, CSourceTerm* m_st)
{

   if (m_st->getProcessType() == OVERLAND_FLOW)
   {
      if (m_st->no_surface_water_pressure)        // 4.10.01
         return factor * (h_cond - z_cond);
      else
         return factor * (h_cond - h_this);
   }                                              // richards & groundwater flow
   else
   {
      if (m_st->getProcessType() == GROUNDWATER_FLOW)
      {
         if (h_this < z_cond && m_st->pcs_type_name_cond == "OVERLAND_FLOW")
            return factor * (h_cond - z_cond);
         else
            return factor * h_cond;
      } else
      // richards flow
      return factor * (h_cond - z_cond);
   }
}
