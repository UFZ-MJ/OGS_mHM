/**************************************************************************
FEMLib-Object: MAT-MP
Task: MediumProperties
Programing:
01/2004 OK Implementation
**************************************************************************/
#include "Configure.h"

//#include "makros.h"
// C++ STL
//#include <iostream>
#include <cfloat>

// FEMLib
#include "tools.h"
//#include "rf_pcs.h"
//#include "femlib.h"
extern double* GEOGetELEJacobianMatrix(long number,double *detjac);
#include "mathlib.h"
//#include "rf_mfp_new.h"
#include "rf_msp_new.h"
//#include "material.h"
#include "rf_tim_new.h"
#include "rfmat_cp.h"
extern double gravity_constant;
//using SolidProp::CSolidProperties;
// LIB
#include "files0.h"
// this
#include "rf_mmp_new.h"
//#include "rf_react.h"
#ifdef RFW_FRACTURE
#include "fem_ele.h"
#endif
#include "gs_project.h"
// Gauss point veclocity
#include "fem_ele_std.h"
#include "fem_ele_vec.h"
// MSHLib
//#include "msh_lib.h"

using namespace std;

// MAT-MP data base lists
list<string>keywd_list;
list<string>mat_name_list;
list<char*>mat_name_list_char;
list<CMediumProperties*>db_mat_mp_list;
// MAT-MP list
vector<CMediumProperties*>mmp_vector;
list<CMediumPropertiesGroup*>mmp_group_list;

using FiniteElement::CElement;
using FiniteElement::ElementValue;
using FiniteElement::CFiniteElementStd;
using FiniteElement::ElementValue_DM;
////////////////////////////////////////////////////////////////////////////
// CMediumProperties
////////////////////////////////////////////////////////////////////////////
// WW
#define GAS_CONSTANT      8314.41
#define COMP_MOL_MASS_AIR    28.96
#define COMP_MOL_MASS_WATER  18.016
#define T_KILVIN_ZERO  273.15                     //AKS
/**************************************************************************
FEMLib-Method: CMediumProperties
Task: constructor
Programing:
02/2004 OK Implementation
last modification:
**************************************************************************/
CMediumProperties::CMediumProperties(void)
{
   name = "DEFAULT";
   mode = 0;
   selected = false;
   m_pcs = NULL;                                  //OK
   // GEO
   geo_dimension = 0;
   geo_area = 1.0;                                //OK
   porosity_model = -1;
   porosity_model_values[0] = 0.1;
   tortuosity_model = -1;
   tortuosity_model_values[0] = 1.;
   // MSH
   m_msh = NULL;                                  //OK
   // flow
   storage_model = -1;
   storage_model_values[0]= 0.;
   permeability_model = -1;
   permeability_tensor_type = 0;
   tortuosity_tensor_type = 0;
   permeability_tensor[0] = 1.e-13;
   conductivity_model = -1;
   flowlinearity_model = 0;
   capillary_pressure_model = -1;
   capillary_pressure_model_values[0] = 0.0;
   permeability_saturation_model[0] = -1;
   permeability_saturation_model_values[0] = 1.0;
   unconfined_flow_group = -1;
   permeability_stress_mode = -1;                 //WW
   c_coefficient = NULL;                          //WW
   // surface flow
   friction_coefficient = -1;
   //  friction_model = -1;
   // mass transport
   // heat transport
   heat_dispersion_model = -1;                    //WW
   heat_dispersion_longitudinal = 0.;
   heat_dispersion_transverse = 0.;
   lgpn = 0.0;
   mass_dispersion_transverse = 0.0;
   mass_dispersion_longitudinal = 0.0;
   heat_diffusion_model = -1;                     //WW
   geo_area = 1.0;
   geo_type_name = "DOMAIN";                      //OK
                                                  //WW
   saturation_max[0] = saturation_max[1] = saturation_max[2] = 1.0;
                                                  //WW
   saturation_res[0] = saturation_res[1] = saturation_res[2] = 0.0;
   permeability_tensor[9] = 1.0e-9;               // Minimum paermeability. 17.06.2008. WW
   vol_mat = 0.0;
   vol_bio = 0.0;
   vol_mat_model = 0;
   vol_bio_model = 0;
   foc = 0.0;
#ifdef RFW_FRACTURE
   frac_num = 0;
   fracs_set = 0;
#endif
}


/**************************************************************************
FEMLib-Method: CMediumProperties
Task: destructor
Programing:
02/2004 OK Implementation
last modification:
**************************************************************************/
CMediumProperties::~CMediumProperties(void)
{
   if(c_coefficient) delete[] c_coefficient;      //WW
   geo_name_vector.clear();
}


////////////////////////////////////////////////////////////////////////////
// IO functions
////////////////////////////////////////////////////////////////////////////

/**************************************************************************
FEMLib-Method: MMPRead
Task: master read functionn
Programing:
02/2004 OK Implementation
last modification:
**************************************************************************/
bool MMPRead(std::string base_file_name)
{
   //----------------------------------------------------------------------
   //OK  MMPDelete();
   //----------------------------------------------------------------------
   std::cout << "MMPRead ... " << std::flush;
   CMediumProperties *m_mat_mp = NULL;
   char line[MAX_ZEILE];
   std::string sub_line;
   std::string line_string;
   std::ios::pos_type position;
   //========================================================================
   // file handling
   std::string mp_file_name;
   mp_file_name = base_file_name + MMP_FILE_EXTENSION;
   std::ifstream mp_file (mp_file_name.data(),std::ios::in);
   if (!mp_file.good())
   {
      std::cout << "! Error in MMPRead: No material data !" << std::endl;
      return false;
   }
   mp_file.seekg(0L,std::ios::beg);
   //========================================================================
   // keyword loop
   while (!mp_file.eof())
   {
      mp_file.getline(line,MAX_ZEILE);
      line_string = line;
      if(line_string.find("#STOP")!=string::npos)
      {
         std::cout << "done, read " << mmp_vector.size() << " medium properties" << std::endl;
         return true;
      }
      //----------------------------------------------------------------------
                                                  // keyword found
      if(line_string.find("#MEDIUM_PROPERTIES")!=string::npos)
      {
         m_mat_mp = new CMediumProperties();
         position = m_mat_mp->Read(&mp_file);
                                                  //OK41
         m_mat_mp->number = (int)mmp_vector.size();
         mmp_vector.push_back(m_mat_mp);
         mp_file.seekg(position,std::ios::beg);
      }                                           // keyword found
   }                                              // eof
   return true;
   // Tests
}


/**************************************************************************
FEMLib-Method: CMediumProperties::Read
Task: read functionn
Programing:
02/2004 OK Template
08/2004 CMCD Implementation
10/2004 MX/OK Porosity model 3, swelling
11/2004 CMCD String streaming
07/2005 MB porosity_file, permeability_file, GEO_TYPE layer
10/2005 OK GEO_TYPE geo_name_vector
01/2006 YD PCS_TYPE
05/2007 PCH Tortuosity tensor
05/2007 WW Stress permeability coorector. Two models.
last modification:
**************************************************************************/
// Order of Key Words
/*
         0. $NAME
            (i)    _BORDEN
         1. $GEOTYPE
            (i)		_CLAY
            (ii)	_SILT
            (iii)	_SAND
            (iv)	_GRAVEL
            (v)		_CRYSTALINE
         2. $GEOMETRY
            (i)		_DIMENSION
(ii)	_AREA
3. $POROSITY
4. $TORTUOSITY
5. $MOBILE_IMOBILE_MODEL
6. $LITHOLOGY_GRAIN_CLASS
7. $FLOWLINEARITY
8. $SORPTION_MODEL
9. $STORAGE
11.$PERMEABILITY
12.$PERMEABILITY_FUNCTION_
(1)		DEFORMATION
(2) 	PRESSURE
(3)	    SATURATION
(4)	    STRESS
(5)		VELOCITY
(6)		POROSITY
13.$CAPILLARY_PRESSURE
14.$MASSDISPERSION
(i)		_LONGITUDINAL
(ii)	_TRANSVERSE
15.$HEATDISPERSION
(i)		_LONGITUDINAL
(ii)	_TRANSVERSE
19.$ELECTRIC_CONDUCTIVITY
20.$UNCONFINED_FLOW_GROUP
21.$FLUID_EXCHANGE_WITH_OTHER_CONTINUA
*/
std::ios::pos_type CMediumProperties::Read(std::ifstream *mmp_file)
{
   int i, j, k=0;
   std::string line_string;
   std::stringstream in;
   std::ios::pos_type position;
   std::string dollar("$");
   std::string hash("#");
   bool new_subkeyword = false;
   bool new_keyword = false;
   std::string m_string;
   // WW
   std::stringstream buff;
   std::vector<string> tokens;
   char *pch;
   char seps[] = "+\n";
   char seps1[] = "*";
   double f_buff;
   size_t indexChWin, indexChLinux;               //JT, DEC 2009
   std::string funfname;                          //JT, DEC 2009

   while (!new_keyword)
   {
      new_subkeyword = false;
      position = mmp_file->tellg();
      line_string = GetLineFromFile1(mmp_file);
      if(line_string.size() < 1) break;
      if(line_string.find(hash)!=std::string::npos)
      {
         new_keyword = true;
         break;
      }
      //--------------------------------------------------------------------
      //PCS                         //YD
      if(line_string.find("$PCS_TYPE")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> pcs_type_name;
         in.clear();
         continue;
      }

      //NAME
                                                  //subkeyword found
      if(line_string.find("$NAME")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> name;                              //sub_line
         in.clear();
         continue;
      }
      //--------------------------------------------------------------------
      //GEO_TYPE

                                                  //subkeyword found
      if(line_string.find("$GEO_TYPE")!=std::string::npos)
      {
         /*
               in.str(GetLineFromFile1(mmp_file));
               in >> geo_type_name;
               in.clear();
               if(geo_type_name.compare("DOMAIN")==0)
                 continue;
               while(!(geo_name.find("$")!=string::npos)&&(!(geo_name.find("#")!=string::npos)))
               {
                 position = mmp_file->tellg();
                 in.str(GetLineFromFile1(mmp_file));
                 in >> geo_name;
         in.clear();
         if(!(geo_name.find("$")!=string::npos))
         geo_name_vector.push_back(geo_name);
         }
         mmp_file->seekg(position,ios::beg);
         continue;
         */
         while(!(m_string.find("$")!=std::string::npos)&&(!(m_string.find("#")!=std::string::npos)))
         {
            position = mmp_file->tellg();
            in.str(GetLineFromFile1(mmp_file));
            in >> m_string >> geo_name;
            in.clear();
            if(!(m_string.find("$")!=std::string::npos)&&(!(m_string.find("#")!=std::string::npos)))
            {
               geo_type_name = m_string;
               geo_name_vector.push_back(geo_name);
            }
         }
         mmp_file->seekg(position,std::ios::beg);
         continue;
      }
      //....................................................................
      // ToDo to GeoLib
      //2i..GEOMETRY_DIMENSION
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$GEOMETRY_DIMENSION")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> geo_dimension;
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      // ToDo to GeoLib
      //2ii..GEOMETRY_AREA
      //------------------------------------------------------------------------
      indexChWin = indexChLinux = 0;              //JT, DEC 2009
                                                  //subkeyword found
      if(line_string.find("$GEOMETRY_AREA")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> line_string;
         if(line_string.find("FILE")!=string::npos)
         {
            in >> geo_area_file;
            // JT, Dec. 16, 2009, added lines below to correct and globalize the read of geometry area file
            std::string file_name = geo_area_file;
            CGSProject* m_gsp = NULL;
            m_gsp = GSPGetMember("mmp");
            if(m_gsp)
            {
               file_name = m_gsp->path + geo_area_file;
            }
            indexChWin = FileName.find_last_of('\\');
            indexChLinux = FileName.find_last_of('/');
            if(indexChWin==string::npos&&indexChLinux==std::string::npos)
               funfname = file_name;
            else if(indexChWin!=string::npos)
            {
               funfname = FileName.substr(0,indexChWin);
               funfname = funfname+"\\"+file_name;
            }
            else if(indexChLinux!=string::npos)
            {
               funfname = FileName.substr(0,indexChLinux);
               funfname = funfname+"/"+file_name;
            }
            geo_area_file = funfname;
            // End of new lines
         }
         else
         {
            geo_area = strtod(line_string.data(),NULL);
         }
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //3..POROSITY
      //------------------------------------------------------------------------
      //CB
                                                  //subkeyword found
      if((line_string.find("$POROSITY")!=string::npos) && (!(line_string.find("_DISTRIBUTION")!=std::string::npos)))
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> porosity_model;
         switch(porosity_model)
         {
            case 0:                               // n=f(x)
               in >> porosity_curve;
               break;
            case 1:                               // n=const
               in >> porosity_model_values[0];
               break;
            case 2:                               // f(normal effective stress for fracture systems)
               in >> porosity_model_values[0];
               in >> porosity_model_values[1];
               in >> porosity_model_values[2];
               in >> porosity_model_values[3];
               porosity_pcs_name_vector.push_back("PRESSURE1");
               break;
            case 3:                               // Chemical swelling model
               in >> porosity_model_values[0];    // Initial porosity
               in >> porosity_model_values[1];    // Specific surface[m^2/g]
               in >> porosity_model_values[2];    // Expansive min. fragtion
               in >> porosity_model_values[3];    // m
               in >> porosity_model_values[4];    // I
               in >> porosity_model_values[5];    // S^l_0
               in >> porosity_model_values[6];    // beta
               porosity_pcs_name_vector.push_back("SATURATION2");
               porosity_pcs_name_vector.push_back("TEMPERATURE1");
               break;
            case 4:                               // Chemical swelling model (constrained swelling, constant I)
               in >> porosity_model_values[0];    // Initial porosity
               in >> porosity_model_values[1];    // Specific surface[m^2/g]
               in >> porosity_model_values[2];    // Expansive min. fragtion
               in >> porosity_model_values[3];    // m
               in >> porosity_model_values[4];    // I
               in >> porosity_model_values[5];    // S^l_0
               in >> porosity_model_values[6];    // beta
               in >> porosity_model_values[7];    // n_min
                                                  //for richard flow only
               porosity_pcs_name_vector.push_back("SATURATION1");
               porosity_pcs_name_vector.push_back("TEMPERATURE1");
               break;
            case 5:                               // Chemical swelling model (free swelling, constant I)
               in >> porosity_model_values[0];    // Initial porosity
               in >> porosity_model_values[1];    // Specific surface[m^2/g]
               in >> porosity_model_values[2];    // Expansive min. fragtion
               in >> porosity_model_values[3];    // m
               in >> porosity_model_values[4];    // I
               in >> porosity_model_values[5];    // S^l_0
               in >> porosity_model_values[6];    // beta
               porosity_pcs_name_vector.push_back("SATURATION2");
               porosity_pcs_name_vector.push_back("TEMPERATURE1");
               break;
            case 6:                               // Chemical swelling model (constrained swelling)
               in >> porosity_model_values[0];    // Initial porosity
               in >> porosity_model_values[1];    // Specific surface[m^2/g]
               in >> porosity_model_values[2];    // Expansive min. fragtion
               in >> porosity_model_values[3];    // m
               in >> porosity_model_values[4];    // I
               in >> porosity_model_values[5];    // S^l_0
               in >> porosity_model_values[6];    // beta
               in >> porosity_model_values[7];    // n_min
               porosity_pcs_name_vector.push_back("SATURATION2");
               porosity_pcs_name_vector.push_back("TEMPERATURE1");
               break;
            case 7:                               // n=f(stress_mean) WW
               in >> porosity_curve;
               break;
            case 10:                              // Chemical swelling model (constrained swelling, constant I)
            {
               int m;
               in >> porosity_model_values[0];    // Initial porosity
               in >> m;                           // m
               if (m > 15)
                  std::cout
                     << "Maximal number of solid phases is now limited to be 15!!!"
                     << std::endl;
               for (int i = 0; i < m + 1; i++)
                                                  // molar volume [l/mol]
                     in >> porosity_model_values[i + 1];
               break;
            }
            case 11:                              //MB: read from file ToDo
               // in >> porosity_file; // CB
               in >> porosity_model_values[0];    // CB some dummy default value is read
               // CB $POROSITY_DISTRIBUTION should be given as keyword in *.mmp file,
               //     porosities then are to be read in from file by fct.
               //     CMediumProperties::SetDistributedELEProperties
               break;
#ifdef GEM_REACT
            case 15:
               in >> porosity_model_values[0];    // set a default value for GEMS calculation
                                                  // save this seperately;
               KC_porosity_initial = porosity_model_values[0];

               // KG44: TODO!!!!!!!!!!!!! check the above  ***************

               break;
#endif
#ifdef BRNS
            case 16:
               in >> porosity_model_values[0];    // set a default value for BRNS calculation
               break;
#endif
            default:
               std::cerr << "Error in MMPRead: no valid porosity model" << std::endl;
               break;

         }
         in.clear();
         continue;
      }
                                                  //subkeyword found
      if(line_string.find("$VOL_MAT")!=string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         // in >> idummy >> this->vol_mat; CB
         in >> vol_mat_model >> this->vol_mat;
         switch(vol_mat_model)
         {
            case 1:                               // do nothing
               break;
            case 2:                               // do nothing
               break;
            default:
               std::cout << "Error in MMPRead: no valid vol_mat_model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
                                                  //subkeyword found
      if(line_string.find("$VOL_BIO")!=string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         // in >> idummy >> this->vol_bio; CB
         in >> vol_bio_model >> this->vol_bio;
         switch(vol_bio_model)
         {
            case 1:                               // do nothing
               break;
            case 2:                               // do nothing
               break;
            default:
               std::cout << "Error in MMPRead: no valid vol_bio_model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //4..TORTUOSITY
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$TORTUOSITY")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> tortuosity_tensor_type_name;
         switch(tortuosity_tensor_type_name[0])
         {
            case '0':                             // n=f(x) <- Case zero
               break;
            case '1':                             // n=const
               tortuosity_model = 1;
               in >> tortuosity_model_values[0];
               break;
            case 'I':                             // isotropic
               tortuosity_tensor_type = 0;
               tortuosity_model = 1;
               in >> tortuosity_model_values[0];
                                                  //CMCD to pick up 2D and 3D Isotropic case.
               tortuosity_model_values[1]=tortuosity_model_values[2]=tortuosity_model_values[0];
               break;
            case 'O':                             // orthotropic  <- Case Alphabet O
               tortuosity_model = 1;
               tortuosity_tensor_type = 1;
               if(geo_dimension==0)
                  std::cout << "Error in CMediumProperties::Read: no geometric dimension" << std::endl;
               if(geo_dimension==2)
               {
                  in >> tortuosity_model_values[0];
                  in >> tortuosity_model_values[1];
               }
               if(geo_dimension==3)
               {
                  in >> tortuosity_model_values[0];
                  in >> tortuosity_model_values[1];
                  in >> tortuosity_model_values[2];
               }
               break;
            case 'A':                             // anisotropic
               tortuosity_model = 1;
               tortuosity_tensor_type = 2;
               if(geo_dimension==0)
                  std::cout << "Error in CMediumProperties::Read: no geometric dimension" << std::endl;
               if(geo_dimension==2)
               {
                  in >> tortuosity_model_values[0];
                  in >> tortuosity_model_values[1];
                  in >> tortuosity_model_values[2];
                  in >> tortuosity_model_values[3];
               }
               if(geo_dimension==3)
               {
                  in >> tortuosity_model_values[0];
                  in >> tortuosity_model_values[1];
                  in >> tortuosity_model_values[2];
                  in >> tortuosity_model_values[3];
                  in >> tortuosity_model_values[4];
                  in >> tortuosity_model_values[5];
                  in >> tortuosity_model_values[6];
                  in >> tortuosity_model_values[7];
                  in >> tortuosity_model_values[8];
               }
               break;
            case 'F':                             //SB: read from file
               tortuosity_model = 2;              //OK
               in >> permeability_file;
               break;
            default:
               std::cout << "Error in MMPRead: no valid tortuosity tensor type" << std::endl;
               break;
         }
         in.clear();
         continue;
      }

      //------------------------------------------------------------------------
      //5..MOBILE_IMOBILE_MODEL
      //------------------------------------------------------------------------
      //To do as necessary

      //------------------------------------------------------------------------
      //6..LITHOLOGY_GRAIN_CLASS
      //------------------------------------------------------------------------
      //To do as necessary

      //------------------------------------------------------------------------
      //7..FLOWLINEARITY
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$FLOWLINEARITY")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> flowlinearity_model;
         switch(flowlinearity_model)
         {
            case 0:                               // k=f(x)
               break;
            case 1:                               // Read in alpha
                                                  //Alpha
               in >> flowlinearity_model_values[0];
               break;
            case 2:                               // For equivalent flow in trianglular elements
                                                  //Alpha
               in >> flowlinearity_model_values[0];
                                                  //Number of Fractures in Equivalent Medium
               in >> flowlinearity_model_values[1];
                                                  //Reynolds Number above which non linear flow occurs.
               in >> flowlinearity_model_values[2];
               pcs_name_vector.push_back("PRESSURE1");
               break;
            default:
               std::cout << "Error in MMPRead: no valid flow linearity model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //8..SORPTION_MODEL
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$ORGANIC_CARBON")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> foc;
         in.clear();
         continue;
      }

      //------------------------------------------------------------------------
      //9..STORAGE
      //------------------------------------------------------------------------

                                                  //subkeyword found
      if(line_string.find("$STORAGE")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> storage_model;
         switch(storage_model)
         {
            case 0:                               // S=f(x)
               in >> storage_model_values[0];     //Function of pressure defined by curve
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 1:                               // S=const
               in >> storage_model_values[0];     //Constant value in Pa
               break;
            case 2:
               in >> storage_model_values[0];     //S0
               in >> storage_model_values[1];     //increase
               in >> storage_model_values[2];     //sigma(z0)
               in >> storage_model_values[3];     //d_sigma/d_z
               break;
            case 3:
               in >> storage_model_values[0];     //curve number (as real number)
               in >> storage_model_values[1];     //sigma(z0)
               in >> storage_model_values[2];     //_sigma/d_z
               break;
            case 4:
               in >> storage_model_values[0];     //curve number (as real number)
               in >> storage_model_values[1];     //time collation
               in >> storage_model_values[2];     //solid density
               in >> storage_model_values[3];     //curve fitting factor, default 1
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 5:                               //Storativity is a function of normal effective stress defined by curve, set up for KTB.
               in >> storage_model_values[0];     //curve number
               in >> storage_model_values[1];     //time collation
               in >> storage_model_values[2];     //Default storage value for material groups > 0
               in >> storage_model_values[3];     //Angular difference between Y direction and Sigma 1
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 6:                               //Storativity is a function of normal effective stress defined by curve and distance from borehole, set up for KTB.
               in >> storage_model_values[0];     //curve number
               in >> storage_model_values[1];     //time collation
               in >> storage_model_values[2];     //Default storage value for material groups > 0
               in >> storage_model_values[3];     //Angular difference between Y direction and Sigma 1
               in >> storage_model_values[4];     //Borehole (x) coordinate
               in >> storage_model_values[5];     //Borehole (y) coordinate
               in >> storage_model_values[6];     //Borehole (z) coordinate
               in >> storage_model_values[7];     //Maximum thickness of shear zone
               in >> storage_model_values[8];     //Fracture density
               pcs_name_vector.push_back("PRESSURE1");
               break;
            default:
               cout << "Error in MMPRead: no valid storativity model" << endl;
               break;
            case 7:                               //RW/WW
               in >> storage_model_values[0];     //Biot's alpha
               in >> storage_model_values[1];     //Skempton's B coefficient
               in >> storage_model_values[2];     //macroscopic drained bulk modulus
               double val_l = storage_model_values[0]* (1.-storage_model_values[0]*storage_model_values[1])
                  /storage_model_values[1]/storage_model_values[2];
               storage_model_values[1] = val_l;
               break;
         }
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //10..CONDUCTIVITY_MODEL
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$CONDUCTIVITY_MODEL")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> conductivity_model;
         switch(conductivity_model)
         {
            case 0:                               // K=f(x)
               break;
            case 1:                               // K=const
               in >> conductivity;
               break;
            case 2:                               // Manning
               break;
            case 3:                               // Chezy
               break;
            default:
               std::cout << "Error in MMPRead: no valid conductivity model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }

                                                  //subkeyword found
      if(line_string.find("$UNCONFINED")!=string::npos)
      {
         unconfined_flow_group = 1;
         continue;
      }

      //------------------------------------------------------------------------
      //11..PERMEABILITY_TENSOR
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$PERMEABILITY_TENSOR")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> permeability_tensor_type_name;
         switch(permeability_tensor_type_name[0])
         {
            case 'I':                             // isotropic
               permeability_tensor_type = 0;
               permeability_model = 1;
               in >> permeability_tensor[0];
                                                  //CMCD to pick up 2D and 3D Isotropic case.
               permeability_tensor[1]=permeability_tensor[2]=permeability_tensor[0];
               break;
            case 'O':                             // orthotropic
               permeability_tensor_type = 1;
               if(geo_dimension==0)
                  std::cout << "Error in CMediumProperties::Read: no geometric dimension" << std::endl;
               if(geo_dimension==2)
               {
                  in >> permeability_tensor[0];
                  in >> permeability_tensor[1];
               }
               if(geo_dimension==3)
               {
                  in >> permeability_tensor[0];
                  in >> permeability_tensor[1];
                  in >> permeability_tensor[2];
               }
               break;
            case 'A':                             // anisotropic
               permeability_tensor_type = 2;
               if(geo_dimension==0)
                  std::cout << "Error in CMediumProperties::Read: no geometric dimension" << std::endl;
               if(geo_dimension==2)
               {
                  in >> permeability_tensor[0];
                  in >> permeability_tensor[1];
                  in >> permeability_tensor[2];
                  in >> permeability_tensor[3];
               }
               if(geo_dimension==3)
               {
                  in >> permeability_tensor[0];
                  in >> permeability_tensor[1];
                  in >> permeability_tensor[2];
                  in >> permeability_tensor[3];
                  in >> permeability_tensor[4];
                  in >> permeability_tensor[5];
                  in >> permeability_tensor[6];
                  in >> permeability_tensor[7];
                  in >> permeability_tensor[8];
               }
               break;
            case 'F':                             //SB: read from file
               permeability_model = 2;            //OK
               in >> permeability_file;
               break;
            default:
               std::cout << "Error in MMPRead: no valid permeability tensor type" << std::endl;
               break;
         }
         in.clear();
         continue;
      }

      //------------------------------------------------------------------------
      //12. $PERMEABILITY_FUNCTION
      //				(i)		_DEFORMATION
      //				(ii)	_PRESSURE
      //				(iii)	_SATURATION
      //				(iv)	_STRESS
      //				(v)		_VELOCITY
      //				(vi)    _POROSITY
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      //12.1 PERMEABILITY_FUNCTION_DEFORMATION
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$PERMEABILITY_FUNCTION_DEFORMATION")!=std::string::npos)
      {
#ifdef RFW_FRACTURE
         relative_permeability_function.push_back("PERMEABILITY_FUNCTION_DEFORMATION");
#endif
         in.str(GetLineFromFile1(mmp_file));
         in >> permeability_model;
         switch(permeability_model)
         {
            case 0:                               // k=f(x)
               break;
            case 1:                               // k=const
               in >> permeability;
               break;
            default:
               std::cout << "Error in MMPRead: no valid permeability model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //12.2 PERMEABILITY_FUNCTION_PRESSURE
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$PERMEABILITY_FUNCTION_PRESSURE")!=std::string::npos)
      {
#ifdef RFW_FRACTURE
         relative_permeability_function.push_back("PERMEABILITY_FUNCTION_PRESSURE");
#endif
         in.str(GetLineFromFile1(mmp_file));
         in >> permeability_pressure_model;
         switch(permeability_pressure_model)
         {
            case 0:                               // k=f(x)
               break;
            case 1:                               // k=const
               in >> permeability_pressure_model_values[0];
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 2:                               // Permebility is a function of effective stress
               in >> permeability_pressure_model_values[0];
               in >> permeability_pressure_model_values[1];
               in >> permeability_pressure_model_values[2];
               in >> permeability_pressure_model_values[3];
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 3:                               // Permeability is a function of non linear flow
               in >> permeability_pressure_model_values[0];
               in >> permeability_pressure_model_values[1];
               in >> permeability_pressure_model_values[2];
               in >> permeability_pressure_model_values[3];
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 4:                               // Function of effective stress from a curve
               in >> permeability_pressure_model_values[0];
               in >> permeability_pressure_model_values[1];
               in >> permeability_pressure_model_values[2];
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 5:                               // Function of overburden converted to effective stress and related to a curve.
               in >> permeability_pressure_model_values[0];
               in >> permeability_pressure_model_values[1];
               in >> permeability_pressure_model_values[2];
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 6:                               // Normal effective stress calculated according to fracture orientation, related over a curve : KTB site
               in >> permeability_pressure_model_values[0];
               in >> permeability_pressure_model_values[1];
               in >> permeability_pressure_model_values[2];
               in >> permeability_pressure_model_values[3];
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 7:                               // Normal effective stress calculated according to fracture orientation, related over a curve : KTB site, and distance to the borehole.
               in >> permeability_pressure_model_values[0];
               in >> permeability_pressure_model_values[1];
               in >> permeability_pressure_model_values[2];
               in >> permeability_pressure_model_values[3];
               in >> permeability_pressure_model_values[4];
               in >> permeability_pressure_model_values[5];
               in >> permeability_pressure_model_values[6];
               in >> permeability_pressure_model_values[7];
               in >> permeability_pressure_model_values[8];
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 8:                               // Effective stress related to depth and curve, Urach
               in >> permeability_pressure_model_values[0];
               in >> permeability_pressure_model_values[1];
               pcs_name_vector.push_back("PRESSURE1");
               break;
            case 9:                               // Effective stress related to depth and curve, and to temperature, Urach
               in >> permeability_pressure_model_values[0];
               in >> permeability_pressure_model_values[1];
               in >> permeability_pressure_model_values[2];
               in >> permeability_pressure_model_values[3];
               in >> permeability_pressure_model_values[4];
               pcs_name_vector.push_back("PRESSURE1");
               pcs_name_vector.push_back("TEMPERATURE1");
               break;
            default:
               std::cout << "Error in MMPRead: no valid permeability model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //12.3 PERMEABILITY_FUNCTION_SATURATION
      //------------------------------------------------------------------------
      //....................................................................
                                                  //subkeyword found
      if(line_string.find("$PERMEABILITY_SATURATION")!=std::string::npos)
      {
#ifdef RFW_FRACTURE
         relative_permeability_function.push_back("PERMEABILITY_SATURATION");
#endif
         int no_fluid_phases =(int)mfp_vector.size();
         for(i=0;i<no_fluid_phases;i++)
         {
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_saturation_model[i];
            switch(permeability_saturation_model[i])
            {
               case 0:                            // k=f(x)
                  in >> permeability_saturation_model_values[i];
                  break;
               case 1:                            // Constant of 1. 26.05.2008. WW
                  in >> saturation_max[i];
                  break;
               case 2:                            // k=S_eff    // Constant of 1. 26.05.2008. WW
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  break;
               case 22:                           // k=1-S_eff // changed to 22 from 21. WW
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  break;
               case 3:                            // power function
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  in >> saturation_exp[i];
                  break;
               case 33:                           // power function
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  in >> saturation_exp[i];
                  break;
               case 4:                            // van Genuchten
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  in >> saturation_exp[i];
                  break;
               case 44:                           // Van Genuchten for no wetting fluid, e.g., gas. WW
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  in >> saturation_exp[i];
                  in >> permeability_tensor[9];   // Minimum permeability. WW
                  break;
               case 6:                            //Brooks-Corey WW
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  in >> saturation_exp[i];
                  break;
               case 66:                           //Brooks-Corey for no wetting fluid, e.g., gas. WW
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  in >> saturation_exp[i];
                  in >> permeability_tensor[9];   // Minimum permeability. WW
                  break;
               case 14:                           // van Genuchten for liquid MX 03.2005 paper swelling pressure
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  in >> saturation_exp[i];
                  break;
               case 15:                           // van Genuchten for gas MX 03.2005 paper swelling pressure
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  in >> saturation_exp[i];
                  break;
               case 16:                           // van Genuchten/Mualem for light oil (B.R. Thoms Masters Thesis 2003) JOD
                  in >> saturation_res[i];
                  in >> saturation_max[i];
                  in >> saturation_exp[i];
                  in >> saturation_alpha[i];
                  break;
               default:
                  std::cout << "Error in MMPRead: no valid permeability saturation model" << std::endl;
                  break;
            }
            in.clear();
         }
         continue;
      }
      /* See case 3 above   CMCD 11/2004...who put this in?, please come and see me so we can put it in together.
          if(line_string.find("$PERMEABILITY_FUNCTION_SATURATION")!=string::npos) { //subkeyword found
           int no_fluid_phases =(int)mfp_vector.size();
           for(i=0;i<no_fluid_phases;i++){
             *mmp_file >> permeability_saturation_model[i];
             *mmp_file >> delimiter;
              switch(permeability_saturation_model[i]) {
                case 0: // k=f(x)
                  break;
                case 3: // power function
                 *mmp_file >> saturation_res[i];
      *mmp_file >> saturation_max[i];
      *mmp_file >> saturation_exp[i];
      break;
      default:
      cout << "Error in MMPRead: no valid permeability saturation model" << endl;
      break;
      }
      mmp_file->ignore(MAX_ZEILE,'\n');
      }
      continue;
      }*/
      //------------------------------------------------------------------------
      //12.4 PERMEABILITY_FUNCTION_STRESS
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$PERMEABILITY_FUNCTION_STRESS")!=std::string::npos)
      {
#ifdef RFW_FRACTURE
         relative_permeability_function.push_back("PERMEABILITY_FUNCTION_STRESS");
#endif
         in.str(GetLineFromFile1(mmp_file));
         in >> permeability_stress_mode;
         switch(permeability_stress_mode)
         {
            case 0:                               // k=f(x)
               break;
            case 1:                               // k=const
               in >> permeability;
               break;
            case 2:                               // Modified LBNL model. WW
               c_coefficient = new double [18];
               in >> c_coefficient[0]             // b_0
                  >> c_coefficient[1]             // alpha
                  >> c_coefficient[2]             // br_1
                  >> c_coefficient[3]             // br_2
                  >> c_coefficient[4]             // br_3
                  >> c_coefficient[5]>>ws;        // fracture frequency
               for(i=6; i<18; i++) c_coefficient[i] = 0.0;
               for(i=0; i<3; i++)
               {
                  in.clear();
                  in.str(GetLineFromFile1(mmp_file));
                  in>>m_string;
                  if(m_string.find("_XX")!=string::npos)
                     k=0;
                  else if(m_string.find("_YY")!=string::npos)
                     k=1;
                  else if(m_string.find("_ZZ")!=string::npos)
                     k=2;
                  m_string.clear();
                  in>>m_string;
                  pch = strtok (const_cast<char*> (m_string.c_str()),seps);
                  buff<<pch;
                  buff>>c_coefficient[6+k*4];
                  buff.clear();
                  while (pch != NULL)
                  {
                     pch = strtok (NULL, seps);
                     if(pch==NULL) break;
                     string token = pch;
                     tokens.push_back(token);
                  }
                  for(j=0; j<(int)tokens.size(); j++)
                  {
                     pch = strtok (const_cast<char*> (tokens[j].c_str()),seps1);
                     buff<<pch;
                     buff>>f_buff;
                     buff.clear();
                     pch = strtok (NULL,seps1);
                     switch(pch[0])
                     {
                        case 'x':  c_coefficient[k*4+7]=f_buff; break;
                        case 'y':  c_coefficient[k*4+8]=f_buff; break;
                        case 'z':  c_coefficient[k*4+9]=f_buff; break;
                     }
                  }
                  tokens.clear();
               }
               break;
            case 3:                               // Barton-Bandis  WW
               c_coefficient = new double [24];
               in >> c_coefficient[0]             // JRC
                  >> c_coefficient[1]             // JCS        //an0
                  >> c_coefficient[2]             // UCS        //Kn
                  >> c_coefficient[7]             // sig_h
                  >> c_coefficient[8]>>ws;        // sig_H
               for(i=3; i<7; i++) c_coefficient[i] = 0.0;
               for(i=9; i<24; i++) c_coefficient[i] = 0.0;
               for(i=0; i<3; i++)
               {
                  in.clear();
                  in.str(GetLineFromFile1(mmp_file));
                  in>>m_string;
                  if(m_string.find("_XX")!=string::npos)
                     k=0;
                  else if(m_string.find("_YY")!=string::npos)
                     k=1;
                  else if(m_string.find("_ZZ")!=string::npos)
                     k=2;
                  m_string.clear();
                  in>>m_string;
                  pch = strtok (const_cast<char*> (m_string.c_str()),seps);
                  buff<<pch;
                  buff>>c_coefficient[9+k*4];
                  buff.clear();
                  while (pch != NULL)
                  {
                     pch = strtok (NULL, seps);
                     if(pch==NULL) break;
                     string token = pch;
                     tokens.push_back(token);
                  }
                  for(j=0; j<(int)tokens.size(); j++)
                  {
                     pch = strtok (const_cast<char*> (tokens[j].c_str()),seps1);
                     buff<<pch;
                     buff>>f_buff;
                     buff.clear();
                     pch = strtok (NULL,seps1);
                     switch(pch[0])
                     {
                        case 'x':  c_coefficient[k*4+10]=f_buff; break;
                        case 'y':  c_coefficient[k*4+11]=f_buff; break;
                        case 'z':  c_coefficient[k*4+12]=f_buff; break;
                     }
                  }
                  tokens.clear();
               }
               //
               CalStressPermeabilityFactor3_Coef();
               //
               break;

            default:
               std::cout << "Error in MMPRead: no valid permeability model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //12.5 PERMEABILITY_FUNCTION_VELOCITY
      //------------------------------------------------------------------------
                                                  //WW
      if(line_string.find("$PERMEABILITY_FUNCTION_VELOCITY")!=std::string::npos)
      {
         //WW   if(line_string.find("$PERMEABILITY_FUNCTION_STRESS")!=string::npos) { //subkeyword found
#ifdef RFW_FRACTURE
         relative_permeability_function.push_back("PERMEABILITY_FUNCTION_STRESS");
#endif
         in.str(GetLineFromFile1(mmp_file));
         in >> permeability_model;
         switch(permeability_model)
         {
            case 0:                               // k=f(x)
               break;
            case 1:                               // k=const
               in >> permeability;
               break;
            default:
               std::cout << "Error in MMPRead: no valid permeability model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //12.6 PERMEABILITY_FUNCTION_POROSITY
      //------------------------------------------------------------------------

                                                  //subkeyword found
      if(line_string.find("$PERMEABILITY_FUNCTION_POROSITY")!=std::string::npos)
      {
#ifdef RFW_FRACTURE
         relative_permeability_function.push_back("PERMEABILITY_FUNCTION_POROSITY");
#endif
         in.str(GetLineFromFile1(mmp_file));
         in >> permeability_porosity_model;
         switch(permeability_porosity_model)
         {
            case 0:                               // k=f(x)
               break;
            case 1:                               // k=const
               in >> permeability_porosity_model_values[0];
               mmp_file->ignore(MAX_ZEILE,'\n');
               break;
            case 2:                               // Model from Ming Lian
               in >> permeability_porosity_model_values[0];
               in >> permeability_porosity_model_values[1];
               break;
            case 3:                               // HS: 11.2008, Kozeny-Carman relationship
               // we set the tensor type first to isotropic and constant as initial value
               permeability_tensor_type = 0;
               permeability_model = 3;            // this means permeability depends on K-C relationship
                                                  // initial values
               in >> permeability_porosity_model_values[0];
               break;
            case 4:                               // HS: 11.2008, Kozeny-Carman_normalized relationship
               // we set the tensor type first to isotropic and constant as initial value
               permeability_tensor_type = 0;
               permeability_model = 4;            // this means permeability depends on K-C_normalized relationship
                                                  // initial values
               in >> permeability_porosity_model_values[0];
               break;
            case 5:                               // HS: 01.2010, Clement 1996 model
               // M. Thullner et al. 2004, J Contaminant Hydrology 70: 37-62, pp42
               permeability_tensor_type = 0;
               permeability_model = 5;            // Clement original model
                                                  // this is initial porosity
               in >> permeability_porosity_model_values[0];
                                                  // this is initial permeability
               in >> permeability_porosity_model_values[1];
               break;
            case 6:                               // HS: 01.2010, ,modified Clement, biomass colonies clogging
               // M. Thullner et al. 2004, J Contaminant Hydrology 70: 37-62, pp42
               permeability_tensor_type = 0;
               permeability_model = 6;            // modified Clement, biomass growing in colonies
                                                  // this is initial porosity
               in >> permeability_porosity_model_values[0];
                                                  // this is initial permeability
               in >> permeability_porosity_model_values[1];
                                                  // this is parameter a
               in >> permeability_porosity_model_values[2];
               break;
            case 7:                               // HS: 01.2010, ,modified Clement, biofilm clogging
               // M. Thullner et al. 2004, J Contaminant Hydrology 70: 37-62, pp42
               permeability_tensor_type = 0;
               permeability_model = 7;            // modified Clement, biomass growing in biofilm
                                                  // this is initial porosity
               in >> permeability_porosity_model_values[0];
                                                  // this is initial permeability
               in >> permeability_porosity_model_values[1];
                                                  // this is parameter b
               in >> permeability_porosity_model_values[2];
                                                  // this is prarameter k_fmin
               in >> permeability_porosity_model_values[3];
               break;
            default:
               std::cout << "Error in MMPRead: no valid permeability model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }

#ifdef RFW_FRACTURE
      //------------------------------------------------------------------------
      //12.7 PERMEABILITY_FRAC_APERTURE
      //------------------------------------------------------------------------

                                                  //subkeyword found
      if(line_string.find("$PERMEABILITY_FRAC_APERTURE")!=string::npos)
      {
         relative_permeability_function.push_back("PERMEABILITY_FRAC_APERTURE");
         in.str(GetLineFromFile1(mmp_file));
         in >> frac_perm_average_type >>roughness;
         in.clear();
         continue;
      }
#endif

      //....................................................................
                                                  //subkeyword found
      if(line_string.find("$CAPILLARY_PRESSURE")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> capillary_pressure_model;
         switch(capillary_pressure_model)
         {
            case 0:                               // k=f(x)
               in >> capillary_pressure_model_values[0];
               break;
            case 1:                               // const
               in >> capillary_pressure;
               break;
            case 4:                               // van Genuchten
               in >> capillary_pressure_model_values[0];
               break;                             //WW

            case 6:                               //  Brook & Corey. 2.05.2008. WW
               in >> capillary_pressure_model_values[0];
               break;                             //WW
            case 16:                              // van Genuchten, separate fit (thoms) JOD
               in >> permeability_exp[0];
               in >> permeability_alpha[0];
               break;
            default:
                                                  //WW
               std::cout << "Error in MMPRead: no valid capillary model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
      //....................................................................
                                                  //Dual Richards
      if(line_string.find("$TRANSFER_COEFFICIENT")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> transfer_coefficient;              //(-)
         //   in >> unsaturated_hydraulic_conductivity;      //(L/T)
         in.clear();
         continue;
      }
      //....................................................................
                                                  //Dual Richards
      if(line_string.find("$SPECIFIC_STORAGE")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> specific_storage;                  //(Pa-1)
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //14 MASSDISPERSION_
      //			(1) LONGITUDINAL
      //			(2) TRANSVERSE
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$MASS_DISPERSION")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> mass_dispersion_model;
         switch(mass_dispersion_model)
         {
            case 0:                               // f(x)
               break;
            case 1:                               // Constant value
               in >> mass_dispersion_longitudinal;
               in >> mass_dispersion_transverse;
               in >> lgpn;
               if(lgpn > 0 ) cout << "      Limiting Grid Peclet Numbers to " << lgpn << endl;
               break;
            default:
               std::cout << "Error in CMediumProperties::Read: no valid mass dispersion model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //15 HEATDISPERSION
      //			(1) LONGTIDUINAL
      //			(2) TRANSVERSE
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$HEAT_DISPERSION")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> heat_dispersion_model;
         switch(heat_dispersion_model)
         {
            case 0:                               // f(x)
               break;
            case 1:                               // Constant value
               in >> heat_dispersion_longitudinal;
               in >> heat_dispersion_transverse;
               break;
            default:
               std::cout << "Error in CMediumProperties::Read: no valid heat dispersion model" << std::endl;
               break;
         }
         in.clear();
         continue;
      }
                                                  //subkeyword found
      if(line_string.find("$DIFFUSION")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> heat_diffusion_model;
         in.clear();
         continue;
      }
                                                  //subkeyword found
      if(line_string.find("$EVAPORATION")!=string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> evaporation;
         in >> heatflux;
         in >> vaporfraction;
         in.clear();
         continue;
      }
      //------------------------------------------------------------------------
      //16. Surface water
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$SURFACE_FRICTION")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> friction_coefficient >> friction_exp_slope >> friction_exp_depth;
         in.clear();
         continue;
      }

                                                  //subkeyword found
      if(line_string.find("$WIDTH")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> overland_width;
         in.clear();
         continue;
      }

                                                  //subkeyword found
      if(line_string.find("$RILL")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> rill_height >> rill_epsilon;
         in.clear();
         continue;
      }

                                                  //subkeyword found
      if(line_string.find("$CHANNEL")!=std::string::npos)
      {
         channel = 1;
         continue;
      }

      //------------------------------------------------------------------------
      //19 ELECTRIC_CONDUCTIVITY
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      //20 UNCONFINED_FLOW_GROUP
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      //21 FLUID_EXCHANGE_WITH_OTHER_CONTINUA
      //------------------------------------------------------------------------

      //------------------------------------------------------------------------
      //11..PERMEABILITY_DISTRIBUTION
      //------------------------------------------------------------------------
      size_t indexChWin, indexChLinux;            //WW
      indexChWin = indexChLinux = 0;
      std::string funfname;
                                                  //subkeyword found
      if(line_string.find("$PERMEABILITY_DISTRIBUTION")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> permeability_file;
         string file_name = permeability_file;
         CGSProject* m_gsp = NULL;
         m_gsp = GSPGetMember("mmp");
         if(m_gsp)
         {
            file_name = m_gsp->path + permeability_file;
         }
         //-------WW
         indexChWin = FileName.find_last_of('\\');
         indexChLinux = FileName.find_last_of('/');
         if(indexChWin==string::npos&&indexChLinux==std::string::npos)
            funfname = file_name;
         else if(indexChWin!=string::npos)
         {
            funfname = FileName.substr(0,indexChWin);
            funfname = funfname+"\\"+file_name;
         }
         else if(indexChLinux!=string::npos)
         {
            funfname = FileName.substr(0,indexChLinux);
            funfname = funfname+"/"+file_name;
         }
         permeability_file = funfname;
         //--------------------------------------
                                                  //WW
         std::ifstream mmp_file(funfname.data(),std::ios::in);
         if (!mmp_file.good())
         {
            std::cout << "Fatal error in MMPRead: no PERMEABILITY_DISTRIBUTION file" << std::endl;
         }
         mmp_file.close();
         permeability_model = 2;
         in.clear();
         continue;
      }

      //------------------------------------------------------------------------
      //11..POROSITY_DISTRIBUTION
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$POROSITY_DISTRIBUTION")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> porosity_file;
         string file_name = porosity_file;
         CGSProject* m_gsp = NULL;
         m_gsp = GSPGetMember("mmp");
         if(m_gsp)
         {
            file_name = m_gsp->path + porosity_file;
         }
         //else{ //CB this is to get the correct path in case the exe is not run from within the project folder
         //  pos = (int)FileName.find_last_of('\\', -1) + 1;
         //  file_name = FileName.substr(0,pos) + porosity_file;
         //}
         //-------CB as above by WW
         indexChWin = FileName.find_last_of('\\');
         indexChLinux = FileName.find_last_of('/');
         if(indexChWin==string::npos&&indexChLinux==std::string::npos)
            funfname = file_name;
         else if(indexChWin!=string::npos)
         {
            funfname = FileName.substr(0,indexChWin);
            funfname = funfname+"\\"+file_name;
         }
         else if(indexChLinux!=string::npos)
         {
            funfname = FileName.substr(0,indexChLinux);
            funfname = funfname+"/"+file_name;
         }
         porosity_file = funfname;
                                                  //WW
         ifstream mmp_file(funfname.data(),ios::in);
         if (!mmp_file.good())
         {
            std::cout << "Fatal error in MMPRead: no POROSITY_DISTRIBUTION file" << std::endl;
         }
         mmp_file.close();
         porosity_model = 11;
         in.clear();
         continue;
      }

      //------------------------------------------------------------------------
#ifdef RFW_FRACTURE
      //------------------------------------------------------------------------
      //FRACTURE DATA, namely how many fractures, what are their names? RFW 11/2005
      //------------------------------------------------------------------------
                                                  //subkeyword found
      if(line_string.find("$FRACTURE_DATA")!=std::string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> frac_num;
         frac_names.clear();
         for(long i=0; i!=frac_num; ++i)
         {
            string read_frac;
            in >> read_frac;

            frac_names.push_back(read_frac);
         }
         frac_perm.resize(frac_names.size());
         avg_aperture.resize(frac_names.size());
         closed_fraction.resize(frac_names.size());
         in.clear();

         continue;
      }
      //-------------------------------------------------------------------------
#endif
   }
   return position;
}


/**************************************************************************
FEMLib-Method: MMPWrite
Task: master write functionn
Programing:
02/2004 OK Implementation
last modification:
**************************************************************************/
void MMPWrite(std::string file_base_name)
//void MATWriteMediumProperties(fstream *mp_file)
{
   CMediumProperties *m_mat_mp = NULL;
   std::string sub_line;
   std::string line_string;
   //========================================================================
   // File handling
   std::string mp_file_name = file_base_name + MMP_FILE_EXTENSION;
   std::fstream mp_file (mp_file_name.c_str(),std::ios::trunc|ios::out);
   mp_file.setf(std::ios::scientific,std::ios::floatfield);
   mp_file.precision(12);
   if (!mp_file.good()) return;
   mp_file << "GeoSys-MMP: Material Medium Properties ------------------------------------------------\n";
   //========================================================================
   // MAT-MP list
   int no_mat =(int)mmp_vector.size();
   int i;
   for(i=0;i<no_mat;i++)
   {
      m_mat_mp = mmp_vector[i];
      m_mat_mp->Write(&mp_file);
   }
   mp_file << "#STOP";
   mp_file.close();
}


/**************************************************************************
FEMLib-Method: CMediumProperties::Write
Task: read functionn
Programing:
02/2004 OK Template
01/2005 OK new formats
05/2005 OK CAPILLARY_PRESSURE,
05/2005 OK Surface flow parameter
06/2005 MB file names
10/2005 OK/YD bugfix model 4
10/2005 OK geo_name_vector
ToDo: MB CONDUCTIVITY_MODEL, PERMEABILITY_SATURATION
**************************************************************************/
void CMediumProperties::Write(std::fstream* mmp_file)
{
   //KEYWORD
   *mmp_file  << "#MEDIUM_PROPERTIES\n";
   //--------------------------------------------------------------------
   //NAME
   if(name.length()>0)
   {
      *mmp_file << " $NAME" << endl;
      *mmp_file << "  ";
      *mmp_file << name << endl;
   }
   //--------------------------------------------------------------------
   //GEO_TYPE
   if(geo_type_name.compare("DOMAIN")==0)         //OK
   {
      *mmp_file << " $GEO_TYPE" << endl;
      *mmp_file << "  DOMAIN" << endl;
   }
   else if((int)geo_name_vector.size()>0)         //OK
   {
      *mmp_file << " $GEO_TYPE" << endl;
      for(int i=0;i<(int)geo_name_vector.size();i++)
      {
         *mmp_file << "  ";
         *mmp_file << geo_type_name;
         *mmp_file << " ";
         *mmp_file << geo_name_vector[i] << endl;
      }
   }
   //--------------------------------------------------------------------
   //DIMENSION
   *mmp_file << " $GEOMETRY_DIMENSION" << endl;
   *mmp_file << "  ";
   *mmp_file << geo_dimension << endl;
   *mmp_file << " $GEOMETRY_AREA" << endl;
   *mmp_file << "  ";
   *mmp_file << geo_area << endl;
   //--------------------------------------------------------------------
   //PROPERTIES
   //....................................................................
   //POROSITY
   if(porosity_model>-1)
   {
      *mmp_file << " $POROSITY" << endl;
      *mmp_file << "  ";
      *mmp_file << porosity_model << " ";
      switch (porosity_model)
      {
         case 1:
            *mmp_file << porosity_model_values[0] << endl;
            break;
         case 11:                                 //MB ToDo
            *mmp_file << porosity_file << endl;   //MB
            break;
      }
   }
   //....................................................................
   //TORTUOSITY //OK4104
   if(tortuosity_model>-1)
   {
      *mmp_file << " $TORTUOSITY" << endl;
      *mmp_file << "  ";
      *mmp_file << tortuosity_model << " ";
      switch (tortuosity_model)
      {
         case 1:
            *mmp_file << tortuosity_model_values[0] << endl;
            break;
      }
   }
   //....................................................................
   //CONDUCTIVITY //MB ToDo
   if(conductivity_model>-1)
   {
      *mmp_file << " $CONDUCTIVITY_MODEL" << endl;
      *mmp_file << "  ";
      *mmp_file << conductivity_model << " ";
      switch (conductivity_model)
      {
         case 1:
            *mmp_file << conductivity;
            break;
      }
      *mmp_file << endl;
   }
   //....................................................................
   //STORAGE
   if(storage_model>-1)
   {
      *mmp_file << " $STORAGE" << endl;
      switch (storage_model)
      {
         case 1:
            *mmp_file << "  " << storage_model << " " << storage_model_values[0] << endl;
            break;
      }
   }
   //....................................................................
   //PERMEABILITY_TENSOR
   if(permeability_tensor_type>-1)
   {
      *mmp_file << " $PERMEABILITY_TENSOR" << endl;
      switch (permeability_tensor_type)
      {
         case 0:
            *mmp_file << "  " << "ISOTROPIC" << " " << permeability_tensor[0] << endl;
            break;
         case 3:                                  //MB for heterogeneous fields //OK
            *mmp_file << "  " << "FILE" << " " << permeability_file << endl;
            break;
      }
   }
   //....................................................................
   //PERMEABILITY_DISTRIBUTION
   if(permeability_model==2)
   {
      *mmp_file << " $PERMEABILITY_DISTRIBUTION" << endl;
      *mmp_file << "  " << permeability_file << endl;
   }
   //....................................................................
   //PERMEABILITY
   if(permeability_saturation_model[0]>-1)
   {
      *mmp_file << " $PERMEABILITY_SATURATION" << endl;
      for(int i=0;i<(int)mfp_vector.size();i++)
      {
         *mmp_file << "  ";
         *mmp_file << permeability_saturation_model[0] << " ";
         switch (permeability_saturation_model[0])
         {
            case 0:
               *mmp_file << (int)permeability_saturation_model_values[0] << endl;
               break;
            case 4:                               //VG
               *mmp_file << saturation_res[0] << " ";
               *mmp_file << saturation_max[0] << " ";
               *mmp_file << saturation_exp[0] << endl;
               break;
         }
      }
   }
   //....................................................................
   //CAPILLARY PRESSURE
   if(capillary_pressure_model>-1)
   {
      *mmp_file << " $CAPILLARY_PRESSURE" << endl;
      *mmp_file << "  ";
      *mmp_file << capillary_pressure_model << " ";
      switch (capillary_pressure_model)
      {
         case 0:
            *mmp_file << (int)capillary_pressure_model_values[0] << endl;
            break;
         case 4:                                  //VG
            *mmp_file << capillary_pressure_model_values[0] << endl;
            break;
         case 9:                                  // HydroSphere
            *mmp_file << saturation_res[0] << " ";
            *mmp_file << capillary_pressure_model_values[0] << " ";
            *mmp_file << capillary_pressure_model_values[1] << endl;
            break;
      }
   }
   //....................................................................
   //HEAT DISPERSION
   if(heat_dispersion_model>-1)
   {
      *mmp_file << " $HEAT_DISPERSION" << endl;
      switch (heat_dispersion_model)
      {
         case 1:
                                                  //CMCD permeability
            *mmp_file << "  " << heat_dispersion_model << " " << heat_dispersion_longitudinal <<" "<<heat_dispersion_transverse<<endl;
            break;
      }
   }
   //....................................................................
   //MASS DISPERSION
   if(mass_dispersion_model>-1)
   {
      *mmp_file << " $MASS_DISPERSION" << endl;
      switch(mass_dispersion_model)
      {
         case 1:
                                                  //CMCD permeability
            *mmp_file << "  " << mass_dispersion_model << " " << mass_dispersion_longitudinal <<" "<<mass_dispersion_transverse<<endl;
            break;
      }
   }

   //----------------------------------------------------------------------
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2004 OK Implementation
last modification:
**************************************************************************/
void MMPWriteTecplot(std::string msh_name)
{
   CMediumProperties *m_mmp = NULL;
   int i;
   int no_mat =(int)mmp_vector.size();
   for(i=0;i<no_mat;i++)
   {
      m_mmp = mmp_vector[i];
      m_mmp->WriteTecplot(msh_name);
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
09/2005 OK MSH
10/2005 OK OO-ELE
last modification:
**************************************************************************/
void CMediumProperties::WriteTecplot(std::string msh_name)
{
   long i,j;
   std::string element_type;
   //----------------------------------------------------------------------
   // GSP
   CGSProject* m_gsp = NULL;
   m_gsp = GSPGetMember("msh");
   if(!m_gsp)
   {
      return;
   }
   //--------------------------------------------------------------------
   // file handling
   std::string mat_file_name = m_gsp->path + "MAT_" + name + TEC_FILE_EXTENSION;
   std::fstream mat_file (mat_file_name.data(),ios::trunc|ios::out);
   mat_file.setf(ios::scientific,ios::floatfield);
   mat_file.precision(12);
   //--------------------------------------------------------------------
   // MSH
   CFEMesh* m_msh = NULL;
   CNode* m_nod = NULL;
   CElem* m_ele = NULL;
   m_msh = FEMGet(msh_name);
   if(!m_msh)
      return;
   //--------------------------------------------------------------------
   if (!mat_file.good()) return;
   mat_file.seekg(0L,std::ios::beg);
   //--------------------------------------------------------------------
   j=0;
   for(i=0;i<(long)m_msh->ele_vector.size();i++)
   {
      m_ele = m_msh->ele_vector[i];
      if(m_ele->GetPatchIndex()==number)
      {
         j++;
      }
   }
   long no_elements = j-1;
   //--------------------------------------------------------------------
   mat_file << "VARIABLES = X,Y,Z,MAT" << std::endl;
   long no_nodes = (long)m_msh->nod_vector.size();
   mat_file << "ZONE T = " << name << ", " \
      << "N = " << no_nodes << ", " \
      << "E = " << no_elements << ", " \
      << "F = FEPOINT" << ", " << "ET = BRICK" << std::endl;
   for(i=0;i<no_nodes;i++)
   {
      m_nod = m_msh->nod_vector[i];
      mat_file \
         << m_nod->X() << " " << m_nod->Y() << " " << m_nod->Z() << " " \
         << number << std::endl;
   }
   j=0;
   for(i=0;i<(long)m_msh->ele_vector.size();i++)
   {
      m_ele = m_msh->ele_vector[i];
      if(m_ele->GetPatchIndex()==number)
      {
         switch(m_ele->GetElementType())
         {
            case MshElemType::LINE:
               mat_file \
                  << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[0]+1 << std::endl;
               j++;
               element_type = "ET = QUADRILATERAL";
               break;
            case MshElemType::QUAD:
               mat_file \
                  << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << " " << m_ele->nodes_index[3]+1 << std::endl;
               j++;
               element_type = "ET = QUADRILATERAL";
               break;
            case MshElemType::HEXAHEDRON:
               mat_file \
                  << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << " " << m_ele->nodes_index[3]+1 << " " \
                  << m_ele->nodes_index[4]+1 << " " << m_ele->nodes_index[5]+1 << " " << m_ele->nodes_index[6]+1 << " " << m_ele->nodes_index[7]+1 << std::endl;
               j++;
               element_type = "ET = BRICK";
               break;
            case MshElemType::TRIANGLE:
               mat_file \
                  << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << endl;
               j++;
               element_type = "ET = TRIANGLE";
               break;
            case MshElemType::TETRAHEDRON:
               mat_file \
                  << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << " " << m_ele->nodes_index[3]+1 << std::endl;
               j++;
               element_type = "ET = TETRAHEDRON";
               break;
            case MshElemType::PRISM:
               mat_file \
                  << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[0]+1 << " " << m_ele->nodes_index[1]+1 << " " << m_ele->nodes_index[2]+1 << " " \
                  << m_ele->nodes_index[3]+1 << " " << m_ele->nodes_index[3]+1 << " " << m_ele->nodes_index[4]+1 << " " << m_ele->nodes_index[5]+1 << std::endl;
               j++;
               element_type = "ET = BRICK";
               break;
            default:
               std::cerr << "CMediumProperties::WriteTecplot MshElemType not handled" << std::endl;
         }
      }
   }
}


////////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////////

/**************************************************************************
FEMLib-Method: CMediumProperties
Task: get instance by name
Programing:
02/2004 OK Implementation
last modification:
**************************************************************************/
CMediumProperties* MMPGet(const std::string &mat_name)
{
   CMediumProperties *m_mat = NULL;
   int no_mat =(int)mmp_vector.size();
   int i;
   for(i=0;i<no_mat;i++)
   {
      m_mat = mmp_vector[i];
      if(mat_name.compare(m_mat->name)==0)
         return m_mat;
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method: CMediumProperties
Task: get instance by name
Programing:
02/2004 OK Implementation
last modification:
**************************************************************************/
CMediumProperties* CMediumProperties::GetByGroupNumber(int group_number)
{
   CMediumProperties *m_mat = NULL;
   int no_mat =(int)mmp_vector.size();
   int i;
   for(i=0;i<no_mat;i++)
   {
      m_mat = mmp_vector[i];
      if(m_mat->number==group_number)
         return m_mat;
   }
   return NULL;
}


////////////////////////////////////////////////////////////////////////////
// Properties functions
////////////////////////////////////////////////////////////////////////////

/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/1997 CT Erste Version
09/1998 CT Kurven ueber Schluesselwort #CURVES
08/1999 CT Erweiterung auf n-Phasen begonnen
09/2002 OK case 13
05/2003 MX case 14
08/2004 OK Template for MMP implementation
05/2005 MX Case 14, 15
05/2005 OK MSH
last modification:
ToDo: GetSoilRelPermSatu
10/2010 TF changed access to process type
**************************************************************************/
double CMediumProperties::PermeabilitySaturationFunction(long number,double*gp,double theta,\
int phase)
{
   //OK411
   theta = theta;
   gp = gp;
   number = number;

   int no_fluid_phases =(int)mfp_vector.size();
   if(!m_pcs)
   {
      cout << "CMediumProperties::PermeabilitySaturationFunction: no PCS data" << endl;
      return 1.0;
   }
   //  if(!(m_pcs->pcs_type_name.compare("RICHARDS_FLOW")==0)&&no_fluid_phases!=2) TF
   if(!(m_pcs->getProcessType () == RICHARDS_FLOW) && no_fluid_phases!=2)
      return 1.0;
   // static int nidx0,nidx1;
   double saturation=0.0,saturation_eff;          //OK411
   int gueltig;
   //---------------------------------------------------------------------
   if(mode==2)                                    // Given value
   {
      saturation = argument;
   }
   else
   {
      /*OK411
          string pcs_name_this = "SATURATION";
          char phase_char[1];
          sprintf(phase_char,"%i",phase+1);
          pcs_name_this.append(phase_char);
          nidx0 = PCSGetNODValueIndex(pcs_name_this,0);
          nidx1 = PCSGetNODValueIndex(pcs_name_this,1);
          if(mode==0){ // Gauss point values
            saturation = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
                     + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
          }
      else{ // Node values
      saturation = (1.-theta)*GetNodeVal(number,nidx0) \
      + theta*GetNodeVal(number,nidx1);
      }
      */
   }
   //----------------------------------------------------------------------
   switch(permeability_saturation_model[phase])
   {
      case 0:                                     // k = f(x) user-defined function
                                                  //CMCD permeability_saturation_model_values[phase]
         permeability_saturation = GetCurveValue((int)permeability_saturation_model_values[phase],0,saturation,&gueltig);
         break;
      case 1:                                     // linear function
         break;
      case 2:                                     // linear function from ! liquid saturation
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
            saturation = saturation_res[phase] + MKleinsteZahl;
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         permeability_saturation = saturation_eff;
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      case 22:                                    // Non-wetting. 1-Se WW
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
            saturation = saturation_res[phase] + MKleinsteZahl;
         saturation_eff = 1.0- (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         permeability_saturation = saturation_eff;
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      case 3:                                     // Nichtlineare Permeabilitaets-Saettigungs-Beziehung
         saturation_eff = (saturation-saturation_res[phase])/(saturation_max[phase]-saturation_res[phase]);
         saturation_eff = MRange(0.,saturation_eff,1.);
         permeability_saturation = pow(saturation_eff,saturation_exp[phase]);
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      case 4:                                     // Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 Page 894 Equation 19 (Burdine's Equation)
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         //SB relperm[1] = pow(s_eff, 2.) * pow((1 - pow((1 - pow(s_eff, 1. / m)), m)),2.0); //SB:todo added ^2N1
         //permeability_saturation = pow(saturation_eff,2.)
         //                        * (1.-pow((1.- pow(saturation_eff,1./saturation_exp[phase])),saturation_exp[phase]));
         // YD Advances in Water Resource 19 (1995) 25-38
         permeability_saturation = pow(saturation_eff,0.5) \
            * pow(1.-pow(1-pow(saturation_eff,1./saturation_exp[phase]),saturation_exp[phase]),2.0);
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      case 5:                                     // Haverkamp Problem: aus SOIL SCI. SOC. AM. J., VOL. 41, 1977, Pages 285ff
         break;
      case 6:                                     // Brooks/Corey: Alte Variante nach Helmig
         break;
      case 7:                                     // Sprungfunktion
         break;
      case 8:                                     // Brooks/Corey: Journal of the Irrigation and Drainage Division June/1966 Pages 61ff
         break;
      case 11:                                    // Drei Phasen ueber Kurven
         break;
      case 12:                                    // Drei Phasen nichtlineare Permeabilitaets-Saettigungs-Beziehung
         break;
      case 13:                                    //OK Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 (Mualems Equation)
         break;
      case 14:                                    //MX Van Genuchten: Wasser aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 (Mualems Equation)
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         //MX relperm[1] = pow(s_eff, 0.5) * pow((1 - pow((1 - pow(s_eff, 1. / m)), m)),2.0);
         permeability_saturation = pow(saturation_eff,0.5) \
            * pow((1.-pow((1.- pow(saturation_eff,1./saturation_exp[phase])),saturation_exp[phase])),2.0);
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      case 15:                                    //MX Van Genuchten: Gas aus Water Resour. Res. VOL. 23, pp2197-2206 1987
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         //MX relperm[1] = pow((1-s_eff), 0.5) * pow((1 - pow(s_eff, 1. / m)),2m);
         permeability_saturation = pow((1-saturation_eff),0.5) \
            * pow((1.-pow(saturation_eff,1./saturation_exp[phase])),2.0*saturation_exp[phase]);
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      case 16:                                    // like 4  JOD vanGebuchten: Thoms MastersThesis 2003

         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         //SB relperm[1] = pow(s_eff, 2.) * pow((1 - pow((1 - pow(s_eff, 1. / m)), m)),2.0); //SB:todo added ^2N1
         //permeability_saturation = pow(saturation_eff,2.)
         //                        * (1.-pow((1.- pow(saturation_eff,1./saturation_exp[phase])),saturation_exp[phase]));

         permeability_saturation = pow(saturation_eff,0.5) \
            * pow(1.-pow(1-pow(saturation_eff,1./permeability_exp[phase]),permeability_exp[phase]),2.0);
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      default:
         std::cout << "Error in CFluidProperties::PermeabilitySaturationFunction: no valid material model" << std::endl;
         break;
   }
   return permeability_saturation;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/1997 CT Erste Version
09/1998 CT Kurven ueber Schluesselwort #CURVES
08/1999 CT Erweiterung auf n-Phasen begonnen
09/2002 OK case 13
05/2003 MX case 14
08/2004 OK Template for MMP implementation
05/2005 MX Case 14, 15
05/2005 OK MSH
08/2005 WW Remove interpolation
03/2007 WW Brooks/Corey:
last modification:
ToDo: GetSoilRelPermSatu
**************************************************************************/
double CMediumProperties::PermeabilitySaturationFunction(const double Saturation,\
int phase)
{
   int no_fluid_phases =(int)mfp_vector.size();
   if(!m_pcs)
   {
      std::cout << "CMediumProperties::PermeabilitySaturationFunction: no PCS data" << std::endl;
      return 1.0;
   }
   //  if(!(m_pcs->pcs_type_name.find("RICHARDS")!=string::npos) && no_fluid_phases!=2) TF
   if(!(m_pcs->getProcessType () == RICHARDS_FLOW) && no_fluid_phases!=2)
      return 1.0;
   double saturation, saturation_eff;
   int gueltig;
   saturation = Saturation;
   //----------------------------------------------------------------------
   switch(permeability_saturation_model[phase])
   {
      case 0:                                     // k = f(x) user-defined function
                                                  //CMCD permeability_saturation_model_values[phase]
         permeability_saturation = GetCurveValue((int)permeability_saturation_model_values[phase],0,saturation,&gueltig);
         break;
      case 1:                                     // constant WW. //linear function
         return 1.0;
         break;
      case 2:                                     // linear function from ! liquid saturation
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
            saturation = saturation_res[phase] + MKleinsteZahl;
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         permeability_saturation = saturation_eff;
         permeability_saturation = MRange(permeability_tensor[9],permeability_saturation,saturation_max[phase]);
         break;
      case 22:                                    // Non-wetting. WW 1-Se
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
            saturation = saturation_res[phase] + MKleinsteZahl;
         saturation_eff = 1.- (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         permeability_saturation = saturation_eff;
                                                  //WW
         permeability_saturation = MRange(permeability_tensor[9],permeability_saturation,saturation_max[phase]);
         break;
      case 3:                                     // Nichtlineare Permeabilitaets-Saettigungs-Beziehung
         saturation_eff = (saturation-saturation_res[phase])/(saturation_max[phase]-saturation_res[phase]);
         saturation_eff = MRange(0.,saturation_eff,1.);
         permeability_saturation = pow(saturation_eff,saturation_exp[phase]);
                                                  //WW
         permeability_saturation = MRange(permeability_tensor[9],permeability_saturation,1.);
         break;
      case 33:                                    // Nichtlineare Permeabilitaets-Saettigungs-Beziehung
         saturation_eff = (saturation-saturation_res[phase])/(saturation_max[phase]-saturation_res[phase]);
         saturation_eff = MRange(0.,saturation_eff,1.);
         permeability_saturation = pow(1.-saturation_eff,saturation_exp[phase]);
         permeability_saturation = MRange(permeability_tensor[9],permeability_saturation,1.);
         break;
      case 4:                                     // Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 Page 894 Equation 19 (Burdine's Equation)
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         //SB relperm[1] = pow(s_eff, 2.) * pow((1 - pow((1 - pow(s_eff, 1. / m)), m)),2.0); //SB:todo added ^2N1
         //permeability_saturation = pow(saturation_eff,2.)
         //                        * (1.-pow((1.- pow(saturation_eff,1./saturation_exp[phase])),saturation_exp[phase]));
         // YD Advances in Water Resource 19 (1995) 25-38
         permeability_saturation = pow(saturation_eff,0.5) \
            * pow(1.-pow(1-pow(saturation_eff,1./saturation_exp[phase]),
            saturation_exp[phase]),2.0);
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      case 44:                                    // Van Genuchtenfor non wetting fluid (e.g. gas):  WW
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
            saturation = saturation_res[phase] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         permeability_saturation = pow(1.0-saturation_eff,1.0/3.0) \
            * pow(1-pow(saturation_eff,1./saturation_exp[phase]),2.0*saturation_exp[phase]);
         permeability_saturation = MRange(permeability_tensor[9],permeability_saturation,1.);
         break;
      case 5:                                     // Haverkamp Problem: aus SOIL SCI. SOC. AM. J., VOL. 41, 1977, Pages 285ff
         break;
      case 6:                                     // Brooks Corey:  WW; changed 27.10.2010. BG + SB, version 5.0.08
         if (saturation > (saturation_max[0] - MKleinsteZahl))
            saturation = saturation_max[0] - MKleinsteZahl;
         if (saturation < (saturation_res[0] + MKleinsteZahl))
            saturation = saturation_res[0] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[0]) / (1.-saturation_res[0]-saturation_res[1]);
         permeability_saturation = pow(saturation_eff,(2.0+3.0*saturation_exp[0])/saturation_exp[0]);
         permeability_saturation = MRange(DBL_EPSILON,permeability_saturation,1.);
         break;
      case 66:                                    // Brooks Corey for non wetting fluid (e.g. gas):  WW; changed 27.10.2010. BG + SB, version 5.0.08
         if (saturation > (saturation_max[0] - MKleinsteZahl))
            saturation = saturation_max[0] - MKleinsteZahl;
         if (saturation < (saturation_res[0] + MKleinsteZahl))
            saturation = saturation_res[0] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[0]) / (1.-saturation_res[0]-saturation_res[1]);
         permeability_saturation = pow(1.0-saturation_eff, 2.0)*
            (1.0-pow(saturation_eff,(2.0+saturation_exp[phase])/saturation_exp[phase]));
         permeability_saturation = MRange(permeability_tensor[9],permeability_saturation,1.);
         break;
      case 7:                                     // Sprungfunktion
         break;
      case 11:                                    // Drei Phasen ueber Kurven
         break;
      case 12:                                    // Drei Phasen nichtlineare Permeabilitaets-Saettigungs-Beziehung
         break;
      case 13:                                    //OK Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 (Mualems Equation)
         break;
      case 14:                                    //MX Van Genuchten: Wasser aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 (Mualems Equation)
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         //MX relperm[1] = pow(s_eff, 0.5) * pow((1 - pow((1 - pow(s_eff, 1. / m)), m)),2.0);
         permeability_saturation = pow(saturation_eff,0.5) \
            * pow((1.-pow((1.- pow(saturation_eff,1./saturation_exp[phase])),saturation_exp[phase])),2.0);
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      case 15:                                    //MX Van Genuchten: Gas aus Water Resour. Res. VOL. 23, pp2197-2206 1987
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         //MX relperm[1] = pow((1-s_eff), 0.5) * pow((1 - pow(s_eff, 1. / m)),2m);
         permeability_saturation = pow((1-saturation_eff),0.5) \
            * pow((1.-pow(saturation_eff,1./saturation_exp[phase])),2.0*saturation_exp[phase]);
         permeability_saturation = MRange(0.,permeability_saturation,1.);
         break;
      case 16:                                    // JOD

         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         //
         saturation_eff = (saturation - saturation_res[phase]) / (saturation_max[phase] - saturation_res[phase]);
         //SB relperm[1] = pow(s_eff, 2.) * pow((1 - pow((1 - pow(s_eff, 1. / m)), m)),2.0); //SB:todo added ^2N1
         //permeability_saturation = pow(saturation_eff,2.)
         //                        * (1.-pow((1.- pow(saturation_eff,1./saturation_exp[phase])),saturation_exp[phase]));
         // YD Advances in Water Resource 19 (1995) 25-38
         permeability_saturation = pow(saturation_eff,0.5) \
            * pow(1.-pow(1-pow(saturation_eff,1./permeability_exp[phase]),permeability_exp[phase]),2.0);
         permeability_saturation = MRange(0.,permeability_saturation,1.);

         break;
      default:
         std::cout << "Error in CFluidProperties::PermeabilitySaturationFunction: no valid material model" << std::endl;
         break;
   }
   return permeability_saturation;
}


/**************************************************************************
FEMLib-Method:
Task: Calculate heat capacity of porous medium
c rho = (c rho)^f + (c rho)^s
Programing:
04/2002 OK Implementation
08/2004 OK MMP implementation
03/2005 WW Case of no fluids
10/2005 YD/OK: general concept for heat capacity
10/2010 TF changed access to process type
**************************************************************************/
double CMediumProperties::HeatCapacity(long number, double theta,
CFiniteElementStd* assem)
{
   CSolidProperties *m_msp = NULL;
   double heat_capacity_fluids, specific_heat_capacity_solid;
   double density_solid;
   double porosity, Sat, PG;
   int group;
   double T0, T1 = 0.0;
   //  double H0,H1;
   // ???
   bool FLOW = false;                             //WW
   for (size_t ii = 0; ii < pcs_vector.size(); ii++)
   {
      // if(pcs_vector[ii]->pcs_type_name.find("FLOW")!=string::npos) {
      if (isFlowProcess (pcs_vector[ii]->getProcessType ()))
      {
         FLOW = true;
         break;
      }

      //WW     if(pcs_vector[ii]->pcs_type_name.find("RICHARDS")!=string::npos)
      //WW       nphase = 3;
   }

   //----------------------------------------------------------------------
   switch (assem->SolidProp->GetCapacityModel())
   {
      //....................................................................
      case 0:                                     // f(x) user-defined function
         break;
         //....................................................................
      case 1:                                     // const
                                                  //OK411
         group = m_pcs->m_msh->ele_vector[number]->GetPatchIndex();
         m_msp = msp_vector[group];
         specific_heat_capacity_solid = m_msp->Heat_Capacity();
         density_solid = fabs(m_msp->Density());
         if (FLOW)
         {
            porosity = assem->MediaProp->Porosity(number, theta);
            heat_capacity_fluids = MFPCalcFluidsHeatCapacity(assem);
         }
         else
         {
            heat_capacity_fluids = 0.0;
            porosity = 0.0;
         }
         heat_capacity = porosity * heat_capacity_fluids + (1.
            - porosity) * specific_heat_capacity_solid
            * density_solid;
         break;
      case 2:                                     //boiling model for YD
         //YD/OK: n c rho = n S^g c^g rho^g + n S^l c^l rho^l + (1-n) c^s rho^s
                                                  //assem->GetNodalVal(1); WW
         T0 = assem->interpolate(assem->NodalVal0);
         //This following lines moved from fem_ele_std but wrong.. WW
         /*
          if(assem->FluidProp->heat_phase_change_curve>0){ //
          if(assem->FluidProp->heat_phase_change_curve>0||assem->heat_phase_change)
          { //
          if(fabs(assem->TG-T0)<1.0e-8) T1 +=1.0e-8;
          H0 = assem->interpolate(assem->NodalVal2);
          H1 = assem->interpolate(assem->NodalVal3);
          heat_capacity = (H1-H0)/(assem->TG-T0);
          }
          else //WW
          {
         */
         if (FLOW)
         {
            PG = assem->interpolate(assem->NodalValC1);
            if (assem->cpl_pcs->type == 1212)     // Multi-phase WW
               PG *= -1.0;
            Sat = SaturationCapillaryPressureFunction(-PG, 0);
         } else
         Sat = 1.0;
         T1 = assem->TG;
         if ((T1 - T0) < DBL_MIN)
            T1 *= -1;
         heat_capacity = assem->SolidProp->Heat_Capacity(T1, Porosity(
            assem), Sat);
         //  }
         break;
      case 3:                                     // D_THM1 - Richards model //WW
         T1 = assem->TG;
         heat_capacity = assem->SolidProp->Heat_Capacity(T1) * fabs(
            assem->SolidProp->Density()) + Porosity(assem)
            * MFPCalcFluidsHeatCapacity(assem);
         break;
         //....................................................................
      default:
         std::cout
            << "Error in CMediumProperties::HeatCapacity: no valid material model"
            << std::endl;
         break;
         //....................................................................
   }
   return heat_capacity;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
04/2002 OK Implementation
08/2004 MB case 5 tet, case 6 prism
09/2004 OK MMP implementation
03/2005 WW Case of no fluids and symmetry
11/2005 CMCD
03/2007 WW Conductivity for multi-phase flow
last modification:
ToDo:
10/2010 TF changed access to process type
**************************************************************************/
double* CMediumProperties::HeatConductivityTensor(int number)
{
   int i, dimen;
   CSolidProperties *m_msp = NULL;
   double heat_conductivity_fluids,Kx[3];
   // TF unused variable - comment fix compile warning
//   double *tensor = NULL;
   //double a, b, Pc, T, Mw, rhow, rho_gw,rho_ga,rho_g, p_gw, mat_fac_w, mat_fac_g, A, B,H_vap, dp_gw, dPc, dA, dB, dT, q,Tc=647.3,expfactor;
   double a, b, rhow, rho_gw,rho_ga,rho_g, p_gw, mat_fac_w, mat_fac_g, A, B,H_vap, dp_gw, dPc, dA, dB, dT, q;
   // TF unused variable - comment fix compile warning
//   double Tc=647.3;
   double expfactor;
   double dens_arg[3]; //AKS
   // TF unused variable - comment fix compile warning
//   ElementValue* gp_ele = ele_gp_value[Fem_Ele_Std->Index];
   //  double porosity =  this->porosity;  //MX
   double Sw, porosity = this->porosity_model_values[0];
   bool FLOW = false;                             //WW
   //  int heat_capacity_model = 0;
   CFluidProperties *m_mfp;                       //WW
   // long group = Fem_Ele_Std->GetMeshElement()->GetPatchIndex();
   m_mfp = Fem_Ele_Std->FluidProp;                //WW

   for (size_t ii = 0; ii < pcs_vector.size(); ii++)
   {
      //		if (pcs_vector[ii]->pcs_type_name.find("FLOW") != string::npos) TF
      if (isFlowProcess (pcs_vector[ii]->getProcessType ()))
         FLOW = true;
   }
   if (FLOW)                                      //WW
   {
      if (Fem_Ele_Std->cpl_pcs->type == 1212)     // Multi-phase WW
      {
         double PG = Fem_Ele_Std->interpolate(
            Fem_Ele_Std->NodalValC1);             // Capillary pressure
         double
            Sw =
            Fem_Ele_Std->MediaProp->SaturationCapillaryPressureFunction(
            PG, 0);
         //
         m_mfp = mfp_vector[0];
         heat_conductivity_fluids = Sw * m_mfp->HeatConductivity();
         m_mfp = mfp_vector[1];
         heat_conductivity_fluids += (1.0 - Sw)
            * m_mfp->HeatConductivity();
      }
      else
      {
         if (Fem_Ele_Std->FluidProp->density_model == 14
            && Fem_Ele_Std->MediaProp->heat_diffusion_model
            == 273 && Fem_Ele_Std->cpl_pcs)
         {
            dens_arg[0] = Fem_Ele_Std->interpolate(
               Fem_Ele_Std->NodalValC1);          //Pressure
            dens_arg[1] = Fem_Ele_Std->interpolate(
               Fem_Ele_Std->NodalVal1) + 273.15;  //Temperature
            dens_arg[2] = Fem_Ele_Std->Index;     //ELE index
            heat_conductivity_fluids
               = Fem_Ele_Std->FluidProp->HeatConductivity(
               dens_arg);
         }
         else
         {
            heat_conductivity_fluids
               = Fem_Ele_Std->FluidProp->HeatConductivity();
         }
         Sw = 1;

         if (Fem_Ele_Std->cpl_pcs->type != 1)
         {
            double PG = Fem_Ele_Std->interpolate(
               Fem_Ele_Std->NodalValC1);          // Capillary pressure

            if (PG < 0.0)
            {
               Sw
                  = Fem_Ele_Std->MediaProp->SaturationCapillaryPressureFunction(
                  -PG, 0);
               heat_conductivity_fluids *= Sw;
               if (Fem_Ele_Std->GasProp != 0)
                  heat_conductivity_fluids
                     += (1. - Sw)
                     * Fem_Ele_Std->GasProp->HeatConductivity();
            }

         }

      }
   }
   else
   {
      heat_conductivity_fluids = 0.0;
      porosity = 0.0;
   }
   dimen = m_pcs->m_msh->GetCoordinateFlag() / 10;
   int group = m_pcs->m_msh->ele_vector[number]->GetPatchIndex();

   for (i = 0; i < dimen * dimen; i++)            //MX
      heat_conductivity_tensor[i] = 0.0;

   m_msp = msp_vector[group];
   m_msp->HeatConductivityTensor(dimen, heat_conductivity_tensor,
      group);                                     //MX

   for (i = 0; i < dimen * dimen; i++)
      heat_conductivity_tensor[i] *= (1.0 - porosity);
   for (i = 0; i < dimen; i++)
      heat_conductivity_tensor[i*dimen+i] += porosity*heat_conductivity_fluids;

   if(evaporation==647)
   {
      int GravityOn = 1;
      if((Fem_Ele_Std->coordinate_system)%10!=2&&(!Fem_Ele_Std->axisymmetry))
         GravityOn = 0;
      double  PG2 = Fem_Ele_Std->interpolate(Fem_Ele_Std->NodalVal_p2);
                                                  // Capillary pressure
      double PG = Fem_Ele_Std->interpolate(Fem_Ele_Std->NodalValC1);
                                                  //Temperature
      double TG=Fem_Ele_Std->interpolate(Fem_Ele_Std->NodalVal1)+273.15;
      double Sw = Fem_Ele_Std->MediaProp->SaturationCapillaryPressureFunction(PG,0);
      heat_conductivity_fluids=0.0;
      H_vap = 2257000;                            //pow((Tc - TG),0.38)*2.65E+5;
      a=19.81;
      b=4975.9;
      m_mfp = mfp_vector[0];
      rhow=m_mfp->Density();
      expfactor = COMP_MOL_MASS_WATER/(rhow*GAS_CONSTANT*TG);
      rho_gw = m_mfp->vaporDensity(TG)*exp(-PG*expfactor);
      p_gw = rho_gw*GAS_CONSTANT*TG/COMP_MOL_MASS_WATER;
      dens_arg[0] = PG2-p_gw;
      dens_arg[1] = TG;
      m_mfp = mfp_vector[1];
      rho_ga = m_mfp->Density(dens_arg);
      rho_g = rho_ga+rho_gw;
      m_mfp = mfp_vector[0];
      mat_fac_w = PermeabilitySaturationFunction(Sw,0)/m_mfp->Viscosity();
      m_mfp = mfp_vector[1];
      mat_fac_g = PermeabilitySaturationFunction(Sw,1)/m_mfp->Viscosity();
      A=b+((PG*COMP_MOL_MASS_WATER)/(rhow*GAS_CONSTANT));
      B=a-log(p_gw/30024.895431831395);
      q=heatflux;

      for(i=0;i<dimen;i++) Kx[i]=0.0;
      dPc=(q/(H_vap*1.0e-13))*((1/(rhow*mat_fac_w)) + (1/(rho_g*mat_fac_g)));
      dA=COMP_MOL_MASS_WATER*dPc/(rhow*GAS_CONSTANT);
      dp_gw=q/(H_vap*rho_gw*mat_fac_g*1.0e-13) ;
      dB=-dp_gw/p_gw;
      dT=(B*dA - A*dB)/(pow(B,2)+(0/TG));
      heat_conductivity_fluids = 2*q/dT;
      Kx[0]=heat_conductivity_fluids;
      if(GravityOn)
      {
         dPc -= (rhow-rho_g)*gravity_constant;
         dA=COMP_MOL_MASS_WATER*dPc/(rhow*GAS_CONSTANT);
         dp_gw -= rho_g*gravity_constant;
         dB=-dp_gw/p_gw;
         dT=(B*dA - A*dB)/(pow(B,2)+(0/TG));
         heat_conductivity_fluids = 2*q/dT;
         Kx[dimen-1]=heat_conductivity_fluids;
      }
      for(i=0;i<dimen*dimen;i++) heat_conductivity_tensor[i] = 0.0;
      for(i=0;i<dimen;i++)
         heat_conductivity_tensor[i*dimen+i] = Kx[i];

   }
   return heat_conductivity_tensor;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
04/2002 OK Implementation
09/2004 OK MMP implementation
12/2005 CMCD
last modification:
ToDo:
**************************************************************************/
double* CMediumProperties::HeatDispersionTensorNew(int ip)
{
   static double heat_dispersion_tensor[9];
   double *heat_conductivity_porous_medium;
   double vg, D[9];
   double dens_arg[3];                            //AKS
   double heat_capacity_fluids=0.0;
   double fluid_density;
   double alpha_t, alpha_l;
   long index = Fem_Ele_Std->GetMeshElement()->GetIndex();
   CFluidProperties *m_mfp;
   int Dim = m_pcs->m_msh->GetCoordinateFlag()/10;
   int i;
   ElementValue* gp_ele = ele_gp_value[index];

   // Materials
                                                  //MX, add index
   heat_conductivity_porous_medium = HeatConductivityTensor(index);
   m_mfp = Fem_Ele_Std->FluidProp;
   if(Fem_Ele_Std->FluidProp->density_model==14 ) //used density changing with p, T
   {
      dens_arg[0]=Fem_Ele_Std->interpolate(Fem_Ele_Std->NodalValC1);
      dens_arg[1]=Fem_Ele_Std->interpolate(Fem_Ele_Std->NodalVal1)+T_KILVIN_ZERO;
      dens_arg[2]=Fem_Ele_Std->Index;
      fluid_density = Fem_Ele_Std->FluidProp->Density(dens_arg);
   }
   else
   {
      fluid_density = m_mfp->Density();
   }
   heat_capacity_fluids = m_mfp->specific_heat_capacity;

   //Global Velocity
   double velocity[3]={0.,0.,0.};
   gp_ele->getIPvalue_vec(ip, velocity);          //gp velocities
   vg = MBtrgVec(velocity,3);

   //Dl in local coordinates
   alpha_l = heat_dispersion_longitudinal;
   alpha_t = heat_dispersion_transverse;

   if (abs(vg) > MKleinsteZahl                    //For the case of diffusive transport only
                                                  //WW
      &&(alpha_l>MKleinsteZahl||alpha_t>MKleinsteZahl) )
   {
      switch (Dim)
      {
         case 1:                                  // line elements
            heat_dispersion_tensor[0] =   heat_conductivity_porous_medium[0] + alpha_l*heat_capacity_fluids*fluid_density*vg;
            break;
         case 2:
            D[0] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[0]*velocity[0])/vg;
            D[1] = ((alpha_l - alpha_t)*(velocity[0]*velocity[1]))/vg;
            D[2] = ((alpha_l - alpha_t)*(velocity[1]*velocity[0]))/vg;
            D[3] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[1]*velocity[1])/vg;
            for (i = 0; i<4; i++)
               heat_dispersion_tensor[i] =   heat_conductivity_porous_medium[i] + (D[i]*heat_capacity_fluids*fluid_density);
            break;
         case 3:
            D[0] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[0]*velocity[0])/vg;
            D[1] = ((alpha_l - alpha_t)*(velocity[0]*velocity[1]))/vg;
            D[2] = ((alpha_l - alpha_t)*(velocity[0]*velocity[2]))/vg;
            D[3] = ((alpha_l - alpha_t)*(velocity[1]*velocity[0]))/vg;
            D[4] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[1]*velocity[1])/vg;
            D[5] = ((alpha_l - alpha_t)*(velocity[1]*velocity[2]))/vg;
            D[6] = ((alpha_l - alpha_t)*(velocity[2]*velocity[0]))/vg;
            D[7] = ((alpha_l - alpha_t)*(velocity[2]*velocity[1]))/vg;
            D[8] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[2]*velocity[2])/vg;
            for (i = 0; i<9; i++)
               heat_dispersion_tensor[i] =   heat_conductivity_porous_medium[i] + (D[i]*heat_capacity_fluids*fluid_density);
            break;
      }
   }
   else
   {
      for (i = 0; i<Dim*Dim; i++)
         heat_dispersion_tensor[i] =   heat_conductivity_porous_medium[i];
   }
   return heat_dispersion_tensor;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
04/2002 OK Implementation
09/2004 OK MMP implementation
10/2004 SB Adapted to mass dispersion
11/2005 CMCD
05/2007 PCH Diffusion Coefficient corrected and Anisotropic diffusion
         coeefient computed with tortuosity added
ToDo:
**************************************************************************/
double* CMediumProperties::MassDispersionTensorNew(int ip)
{
   static double advection_dispersion_tensor[9];  //Name change due to static conflict
   int component = Fem_Ele_Std->pcs->pcs_component_number;
   int i;
   long index = Fem_Ele_Std->GetMeshElement()->GetIndex();
   double molecular_diffusion[9], molecular_diffusion_value;
   double vg;
   double D[9];
   double alpha_l,alpha_t;
   double theta = Fem_Ele_Std->pcs->m_num->ls_theta;
   double g[3]={0.,0.,0.};
   double l_char=0.0;                             //OK411 volume=0.0;
   int set=0;
   ElementValue* gp_ele = ele_gp_value[index];
   CompProperties *m_cp = cp_vec[component];
   MshElemType::type eleType = m_pcs->m_msh->ele_vector[number]->GetElementType();
   int Dim = m_pcs->m_msh->GetCoordinateFlag()/10;
   //----------------------------------------------------------------------
   // Materials
   molecular_diffusion_value = m_cp->CalcDiffusionCoefficientCP(index,theta, m_pcs) * TortuosityFunction(index,g,theta);
   molecular_diffusion_value *= Porosity(index,theta);
   //CB
   molecular_diffusion_value *= PCSGetEleMeanNodeSecondary_2(index, Fem_Ele_Std->pcs->flow_pcs_type, "SATURATION1", 1);
   //if(PCSGet("RICHARDS_FLOW")) molecular_diffusion_value *= PCSGetEleMeanNodeSecondary(index, "RICHARDS_FLOW", "SATURATION1", 1);
   for (i = 0; i<Dim*Dim; i++)
      molecular_diffusion[i] = 0.0;
   for (i = 0; i<Dim; i++)
   {
      molecular_diffusion[i*Dim+i] = molecular_diffusion_value;
   }
   //----------------------------------------------------------------------

   //Anisotropic diffusion coefficient
   if(tortuosity_tensor_type_name[0] == 'A')
   {
      switch (Dim)
      {
         //--------------------------------------------------------------------
         case 1:                                  // line elements
            ;                                     // Do nothing
            break;
            //--------------------------------------------------------------------
         case 2:
            molecular_diffusion[0] = molecular_diffusion_value*tortuosity_model_values[0];
            molecular_diffusion[1] = molecular_diffusion_value*tortuosity_model_values[1];
            molecular_diffusion[2] = molecular_diffusion_value*tortuosity_model_values[2];
            molecular_diffusion[3] = molecular_diffusion_value*tortuosity_model_values[3];
            break;
            //--------------------------------------------------------------------
         case 3:
            molecular_diffusion[0] = molecular_diffusion_value*tortuosity_model_values[0];
            molecular_diffusion[1] = molecular_diffusion_value*tortuosity_model_values[1];
            molecular_diffusion[2] = molecular_diffusion_value*tortuosity_model_values[2];
            molecular_diffusion[3] = molecular_diffusion_value*tortuosity_model_values[3];
            molecular_diffusion[4] = molecular_diffusion_value*tortuosity_model_values[4];
            molecular_diffusion[5] = molecular_diffusion_value*tortuosity_model_values[5];
            molecular_diffusion[6] = molecular_diffusion_value*tortuosity_model_values[6];
            molecular_diffusion[7] = molecular_diffusion_value*tortuosity_model_values[7];
            molecular_diffusion[8] = molecular_diffusion_value*tortuosity_model_values[8];
      }
   }

   //Global Velocity
   double velocity[3]={0.,0.,0.};
   gp_ele->getIPvalue_vec(ip, velocity);          //gp velocities
   vg = MBtrgVec(velocity,3);
   //  if(index < 10) cout <<" Velocity in MassDispersionTensorNew(): "<<velocity[0]<<", "<<velocity[1]<<", "<<velocity[2]<<", "<<vg<<endl;
   // test bubble velocity
   if(m_cp->bubble_velocity_model == 1)
   {
      if(index < 0) cout <<" Velocity in MassDispersionTensorNew(): "<<velocity[0]<<", "<<velocity[1]<<", "<<velocity[2]<<", "<<vg<<endl;
      velocity[0] = m_cp->bubble_velocity[0];
      velocity[1] = m_cp->bubble_velocity[1];
      velocity[2] = m_cp->bubble_velocity[2];
      vg = MBtrgVec(velocity,3);
      if(index == 100) std::cout <<" Bubble velocity in MassDispersionTensorNew(): "<<velocity[0]<<", "<<velocity[1]<<", "<<velocity[2]<<", "<<vg<<std::endl;
   }
   // end bubble velocity
   //Dl in local coordinates
   alpha_l = mass_dispersion_longitudinal;
   alpha_t = mass_dispersion_transverse;

   // hard stabilization
   if(this->lgpn > 0.0)
   {
      CElem* m_ele = NULL;
      m_ele = m_pcs->m_msh->ele_vector[index];
      if(eleType == 2)   l_char = sqrt(m_ele->GetVolume());
      if(eleType == 4)   l_char = sqrt(m_ele->GetVolume());
      // cout << " Element number: " << index << ", Volume: " << m_ele->GetVolume() << ", l_char: " << l_char << endl;
      set=0;
      if(alpha_l < l_char/lgpn)
      {
         set = 1;                                 //flag for output
         alpha_l = l_char/lgpn;
      }
      if(alpha_t < l_char/lgpn)
      {
         set = 1;
         alpha_t = l_char/lgpn;
      }

                                                  //cout << " alpha_L = " << alpha_l << " < l_char/Pe; setting alpha_L = " << l_char/lgpn << " for element " << index << endl;
      if((set > 0)&(aktueller_zeitschritt==1)&(component<1)&(ip<1)) std::cout << "element " << index << " " << l_char << " " << alpha_l << " " << alpha_t <<  std::endl;
   }
   //----------------------------------------------------------------------

   if (abs(vg) > MKleinsteZahl)                   //For the case of diffusive transport only.
   {
      switch (Dim)
      {
         //--------------------------------------------------------------------
         case 1:                                  // line elements
            advection_dispersion_tensor[0] = molecular_diffusion[0] + alpha_l*vg;
            break;
            //--------------------------------------------------------------------
         case 2:
            D[0] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[0]*velocity[0])/vg;
            D[1] = ((alpha_l - alpha_t)*(velocity[0]*velocity[1]))/vg;
            D[2] = ((alpha_l - alpha_t)*(velocity[1]*velocity[0]))/vg;
            D[3] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[1]*velocity[1])/vg;

            if(tortuosity_tensor_type_name[0] == 'A')
            {
               for (i = 0; i<Dim*Dim; i++)
                  advection_dispersion_tensor[i] = molecular_diffusion[i] + D[i];
            }
            else
            {
               advection_dispersion_tensor[0] = molecular_diffusion[0] + D[0];
                                                  //SB added - CHP, all parts of tensor required
               advection_dispersion_tensor[1] = D[1];
                                                  //SB added - CHP, all parts of tensor required
               advection_dispersion_tensor[2] = D[2];
               advection_dispersion_tensor[3] = molecular_diffusion[3] + D[3];
            }
            break;
            //--------------------------------------------------------------------
         case 3:
            D[0] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[0]*velocity[0])/vg;
            D[1] = ((alpha_l - alpha_t)*(velocity[0]*velocity[1]))/vg;
            D[2] = ((alpha_l - alpha_t)*(velocity[0]*velocity[2]))/vg;
            D[3] = ((alpha_l - alpha_t)*(velocity[1]*velocity[0]))/vg;
            D[4] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[1]*velocity[1])/vg;
            D[5] = ((alpha_l - alpha_t)*(velocity[1]*velocity[2]))/vg;
            D[6] = ((alpha_l - alpha_t)*(velocity[2]*velocity[0]))/vg;
            D[7] = ((alpha_l - alpha_t)*(velocity[2]*velocity[1]))/vg;
            D[8] = (alpha_t * vg) + (alpha_l - alpha_t)*(velocity[2]*velocity[2])/vg;
            if(tortuosity_tensor_type_name[0] == 'A')
            {
               for (i = 0; i<Dim*Dim; i++)
                  advection_dispersion_tensor[i] = molecular_diffusion[i] + D[i];
            }
            else
            {
               /* SB changed - CHP, check - aktually the main coordinate method should work
               advection_dispersion_tensor[0] = molecular_diffusion[0] + D[0];
               advection_dispersion_tensor[4] = molecular_diffusion[4] + D[4];
               advection_dispersion_tensor[8] = molecular_diffusion[8] + D[8];
               */
               for (i = 0; i<Dim*Dim; i++)
                  advection_dispersion_tensor[i] =  D[i];
               advection_dispersion_tensor[0] += molecular_diffusion[0];
               advection_dispersion_tensor[4] += molecular_diffusion[4];
               advection_dispersion_tensor[8] += molecular_diffusion[8];

            }
      }
   }
   else
   {
      for (i = 0; i<Dim*Dim; i++)
         advection_dispersion_tensor[i] = molecular_diffusion[i];
   }
   return advection_dispersion_tensor;
}


////////////////////////////////////////////////////////////////////////////
// DB functions
////////////////////////////////////////////////////////////////////////////

/**************************************************************************
FEMLib-Method: CMediumProperties
Task: get instance by name
Programing:
02/2004 OK Implementation
last modification:
**************************************************************************/
CMediumProperties* CMediumProperties::GetDB(std::string mat_name)
{
   CMediumProperties *m_mat = NULL;
   std::list<CMediumProperties*>::const_iterator p_mat = db_mat_mp_list.begin();
   while(p_mat!=db_mat_mp_list.end())
   {
      m_mat = *p_mat;
      if(mat_name.compare(m_mat->name)==0)
         return m_mat;
      ++p_mat;
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method: CMediumProperties
Task: set properties
Programing:
02/2004 OK Implementation
last modification:
**************************************************************************/
void CMediumProperties::SetDB(std::string mat_name,std::string prop_name,double value)
{
   CMediumProperties *m_mat = NULL;
   m_mat = GetDB(mat_name);
   switch (GetPropertyType(prop_name))
   {
      case 0:
         m_mat->conductivity = value;
         break;
      case 1:
         m_mat->permeability = value;
         break;
      case 2:
         m_mat->porosity = value;
         break;
   }
}


/**************************************************************************
FEMLib-Method: CMediumProperties::GetPropertyType
Task: get property type
Programing:
02/2004 OK Implementation
last modification:
**************************************************************************/
int CMediumProperties::GetPropertyType(std::string prop_name)
{
   int counter=0;
   string keyword_name;
   list<string>::const_iterator p_keywords = keywd_list.begin();
   while(p_keywords!=keywd_list.end())
   {
      keyword_name = *p_keywords;
      if(prop_name.compare(keyword_name)==0)
         return counter;
      ++p_keywords,
         counter++;
   }
   return -1;
}


////////////////////////////////////////////////////////////////////////////
// MAT-MP data base
////////////////////////////////////////////////////////////////////////////
std::string read_MAT_name(std::string in, std::string *z_rest_out)
{
   std::string mat_name;
   std::string z_rest;
   std::string delimiter(";");
                                                  //wenn eine mg gefunden wird nach dem Schlüsselwort
   if(in.find_first_not_of(delimiter)!=std::string::npos)
   {
      z_rest = in.substr(in.find_first_not_of(delimiter));
                                                  //string für die erste (1) material group
      mat_name = z_rest.substr(0,z_rest.find_first_of(delimiter));
      *z_rest_out = z_rest.substr(mat_name.length());
      return mat_name;
   }
   else
      return "";
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2004 OK/JG Implementation
last modification:
**************************************************************************/
void MATLoadDB(std::string csv_file_name)
{
   std::string keyword("MATERIALS");
   std::string in;
   std::string line;
   std::string z_rest;
   std::string mat_name;
   std::string mat_name_tmp("MAT_NAME");
   //========================================================================
   // Read available MAT-MP keywords
   read_keywd_list();
   //========================================================================
   // File handling
   ifstream csv_file(csv_file_name.data(),ios::in);
   if (!csv_file.good()) return;
   csv_file.seekg(0L,ios::beg);                   // rewind
   //========================================================================
   // Read MATERIALS group names
   comp_keywd_list(csv_file_name);
   /*
      //iostream
      // 2.1 Create MAT-MP instances
      // search "MATERIALS"
      SOIL_PROPERTIES *sp = NULL;
      sp = MATAddSoilPropertiesObject(); //OK
      strcpy(sp->mat_type_name,"Brooklawn");
      // 2.2 Insert to db_mat_mp_list
      db_material_mp_list.push_back(sp);
      // 2.3 Read material properties
      string line;
   mat_read_hydraulic_conductivity_mean(line,sp->mat_type_name);
   // select corresponding MAT-MP
   }
   */
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2004 JG Implementation
last modification:
**************************************************************************/
void read_keywd_list(void)
{
   //File handling=============================
   std::ifstream eingabe("mat_mp_keywords.dat",ios::in);
   if (eingabe.good())
   {
      eingabe.seekg(0L,ios::beg);
      //==========================================
      std::string keyword("MATERIALS");
      std::string in;
      std::string line;
      std::string z_rest;
      std::string mat_name;
      std::string delimiter(";");
      std::string keywd_tmp("KEYWD");
      char line_char[MAX_ZEILE];
      // Read MATERIALS group names
      //1 get string with keywords
      while (!eingabe.eof())
      {
         eingabe.getline(line_char, MAX_ZEILE);
         line = line_char;
         if(line.find_first_not_of(delimiter)!=string::npos)
         {
                                                  //schneidet delimiter ab
            in = line.substr(line.find_first_not_of(delimiter));
            keywd_tmp = in.substr(0,in.find_first_of(delimiter));
            keywd_list.push_back(keywd_tmp);
         }
         //keywd_list.remove(keyword);
      }                                           // eof
   }                                              // if eingabe.good
   else
   {
      printf("No keyword file: mat_mp_keywords.dat");
   }
   return;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2004 JG/OK Implementation
last modification:
**************************************************************************/
void comp_keywd_list(std::string csv_file_name)
{
   std::string mat_names("MATERIALS");
   std::string in;
   std::string line;
   std::string z_rest;
   std::string mat_name;
   std::string mat_name_tmp("MAT_NAME");
   char line_char[MAX_ZEILE];
   std::string in1;                               //zwischenstring zum abschneiden der einheit
   std::string delimiter(";");
   std::string keyword("MATERIALS");
   double kwvalue;

   //File handling------------------------------------
   std::string sp;
   std::ifstream eingabe(csv_file_name.data(),ios::in);
   if (eingabe.good())
   {
      eingabe.seekg(0L,ios::beg);                 //rewind um materialgruppen auszulesen
      eingabe.getline(line_char, MAX_ZEILE);
      line = line_char;
      //-----------------------------------------------
      // 1 - read MAT names
      // ReadMATNames(csv_file_name);
      if(line.find(mat_names)!=string::npos)
      {
         in = line.substr(mat_names.length()+1);
         // 2 Read material group names
         while(!mat_name_tmp.empty())
         {
            mat_name_tmp = read_MAT_name(in,&z_rest);
            if(mat_name_tmp.empty())
               break;
            else
            {
               mat_name = mat_name_tmp;
               mat_name_list.push_back(mat_name);
               in = z_rest;
            }
         }                                        // while mat_name
      }                                           // keyword found
      //-----------------------------------------------
      // 2 - create MAT-SP instances
      CMediumProperties *m_mat_mp = NULL;
      list<string>::const_iterator p_mat_names=mat_name_list.begin();
      while(p_mat_names!=mat_name_list.end())
      {
         mat_name = *p_mat_names;
         m_mat_mp = new CMediumProperties;
         m_mat_mp->name = mat_name;
         db_mat_mp_list.push_back(m_mat_mp);
         ++p_mat_names;
      }
      //-----------------------------------------------
      // 3 - read MAT properties

      //1 get string where keyword hits
      CMediumProperties *m_mat_mp1 = NULL;
      while (!eingabe.eof())
      {
         eingabe.getline(line_char, MAX_ZEILE);
         line = line_char;
         string sp;
         list<string>::const_iterator pm =keywd_list.begin();
         while(pm!=keywd_list.end())
         {
            sp = *pm;
            if(line.find(sp)!=string::npos)
            {
                                                  //schneidet keyword ab
               in1 = line.substr(line.find_first_not_of(delimiter)+sp.length()+1);
                                                  //schneidet keyword ab
               in = in1.substr(in1.find_first_of(delimiter));
               //Schleife über alle MAT-Gruppen
               list<string>::const_iterator p_mat_names=mat_name_list.begin();
               while(p_mat_names!=mat_name_list.end())
               {
                  mat_name = *p_mat_names;
                  m_mat_mp1 = m_mat_mp1->GetDB(mat_name);
                  mat_name_tmp = read_MAT_name(in,&z_rest);
                  kwvalue = atof(mat_name_tmp.data());
                  //Val = strtod(pDispInfo->item.strText, NULL);
                  m_mat_mp1->SetDB(mat_name,sp,kwvalue);
                  ++p_mat_names;
               }
            }                                     //
            ++pm;
         }                                        // kwlist
      }                                           // eof
   }                                              // if eingabe.good
   return;
}


///*************************************************************************************************
////////////////////////////////////////////////////////////////////////////
// Properties functions
////////////////////////////////////////////////////////////////////////////

/*************************************************************************
 ROCKFLOW - Funktion: SetSoilPropertiesDefaultsClay

 Aufgabe:

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

 Ergebnis:
   - void -

Programmaenderungen:
10/2003   CMCD   Erste Version

*************************************************************************/

void CMediumProperties::SetMediumPropertiesDefaultsClay(void)
{
   int i;
   geo_dimension = 3;
   geo_area = 1.0;
   porosity = 0.4;
   tortuosity = 1.0;
   storage = 1.5e-10;                             //m/pa
   for (i = 0; i < 9; i++) permeability_tensor[i]=1.e-17;
   heat_dispersion_longitudinal = 0.01;
   heat_dispersion_transverse = 0.01;
   mass_dispersion_longitudinal = 0.01;
   mass_dispersion_transverse = 0.01;
   for (i = 0; i < 9; i++) heat_conductivity_tensor[i] = 3.2;

   return;
}


/*************************************************************************
 ROCKFLOW - Funktion: SetSoilPropertiesDefaultsSilt

 Aufgabe:

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

 Ergebnis:
   - void -

Programmaenderungen:
10/2003   CMCD   Erste Version

*************************************************************************/

void CMediumProperties::SetMediumPropertiesDefaultsSilt(void)
{
   int i;
   geo_dimension = 3;
   geo_area = 1.0;
   porosity = 0.4;
   tortuosity = 1.0;
   storage = 5.e-7;                               //m/pa
                                                  //m²
   for (i = 0; i < 9; i++) permeability_tensor[i]=1.e-14;
   heat_dispersion_longitudinal = 0.01;
   heat_dispersion_transverse = 0.01;
   mass_dispersion_longitudinal = 0.01;
   mass_dispersion_transverse = 0.01;
   for (i = 0; i < 9; i++) heat_conductivity_tensor[i] = 3.2;

   return;
}


/*************************************************************************
 ROCKFLOW - Funktion: SetSoilPropertiesDefaultsSand

 Aufgabe:

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

 Ergebnis:
   - void -

Programmaenderungen:
10/2003   CMCD   Erste Version

*************************************************************************/
void CMediumProperties::SetMediumPropertiesDefaultsSand(void)
{
   int i;
   geo_dimension = 3;
   geo_area = 1.0;
   porosity = 0.3;
   tortuosity = 1.0;
   storage = 5.e-8;                               //m/pa
                                                  //m²
   for (i = 0; i < 9; i++) permeability_tensor[i]=1.e-11;
   heat_dispersion_longitudinal = 0.01;
   heat_dispersion_transverse = 0.01;
   mass_dispersion_longitudinal = 0.01;
   mass_dispersion_transverse = 0.01;
   for (i = 0; i < 9; i++) heat_conductivity_tensor[i] = 3.2;

   return;
}


/*************************************************************************
 ROCKFLOW - Funktion: SetSoilPropertiesDefaultsGravel

 Aufgabe:

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

 Ergebnis:
   - void -

Programmaenderungen:
10/2003   CMCD   Erste Version

*************************************************************************/

void CMediumProperties::SetMediumPropertiesDefaultsGravel(void)
{
   int i;
   geo_dimension = 3;
   geo_area = 1.0;
   porosity = 0.32;
   tortuosity = 1.0;
   storage = 5.e-9;                               //m/pa
                                                  //m²
   for (i = 0; i < 9; i++) permeability_tensor[i]=1.e-9;
   heat_dispersion_longitudinal = 0.01;
   heat_dispersion_transverse = 0.01;
   mass_dispersion_longitudinal = 0.01;
   mass_dispersion_transverse = 0.01;
   for (i = 0; i < 9; i++) heat_conductivity_tensor[i] = 3.2;

   return;
}


/*************************************************************************
 ROCKFLOW - Funktion: SetMediumPropertiesDefaultsCrystaline

 Aufgabe:

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

 Ergebnis:
   - void -

Programmaenderungen:
10/2003   CMCD   Erste Version

*************************************************************************/
void CMediumProperties::SetMediumPropertiesDefaultsCrystalline(void)
{
   int i;
   geo_dimension = 3;
   geo_area = 1.0;
   porosity = 0.05;
   tortuosity = 1.0;
   storage = 5.e-11;                              //m/pa
                                                  //m²
   for (i = 0; i < 9; i++) permeability_tensor[i]=1.e-14;
   heat_dispersion_longitudinal = 0.01;
   heat_dispersion_transverse = 0.01;
   mass_dispersion_longitudinal = 0.01;
   mass_dispersion_transverse = 0.01;
   for (i = 0; i < 9; i++) heat_conductivity_tensor[i] = 3.2;

   return;
}


/*************************************************************************
 ROCKFLOW - Funktion: SetMediumPropertiesDefaultsCrystaline

 Aufgabe:

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

 Ergebnis:
   - void -

Programmaenderungen:
10/2003   CMCD   Erste Version

*************************************************************************/
void CMediumProperties::SetMediumPropertiesDefaultsBordenAquifer(void)
{
   int i;
   geo_dimension = 3;
   geo_area = 1.0;
   porosity = 0.10;
   tortuosity = 1.0;
   storage = 5.e-11;                              //m/pa
                                                  //m²
   for (i = 0; i < 9; i++) permeability_tensor[i]=1.e-12;
   heat_dispersion_longitudinal = 0.01;
   heat_dispersion_transverse = 0.01;
   mass_dispersion_longitudinal = 0.01;
   mass_dispersion_transverse = 0.01;
   for (i = 0; i < 9; i++) heat_conductivity_tensor[i] = 3.2;

   return;
}


/*------------------------------------------------------------------------*/
/* MAT-MP geometric properties */
/* 3 porosity */
/*------------------------------------------------------------------------*/

/**************************************************************************
FEMLib-Method:
Task: Porosity calc function
  Case overview
  0 Curve function
  1 Constant Value
  2 Function of normal stress from Geomechnical model
  3 Free swelling
  4 Constraint swelling
Programing:
07/2004 OK C++ Implementation
08/2004	CMCD Re-written based on MATCalcPorosity
10/2004 MX Swelling processes
04/2007 WW Porosity by gauss stress value
last modification:
*************************************************************************/
double CMediumProperties::Porosity(long number,double theta)
{
   static int nidx0, nidx1;
   double primary_variable[PCS_NUMBER_MAX];
   int gueltig;
#ifdef GEM_REACT
   int idx;
#endif
   double porosity_sw;
   CFiniteElementStd* assem = m_pcs->GetAssember();
   string str;
   ///
   ElementValue_DM* gval = NULL;

   // CB Get idx of porosity in elements mat vector for het porosity
   size_t por_index (0);
   if (porosity_model == 11)
   {
      for (por_index = 0; por_index
         < m_pcs->m_msh->mat_names_vector.size(); por_index++)
      {
         if (m_pcs->m_msh->mat_names_vector[por_index].compare("POROSITY")
            == 0)
         {
            break;
         }
      }
   }

   // Functional dependencies
   CRFProcess *pcs_temp;

   const size_t no_pcs_names(this->porosity_pcs_name_vector.size());
   for (size_t i = 0; i < no_pcs_names; i++)
   {
      str = porosity_pcs_name_vector[i];
      pcs_temp = PCSGet(str, true);               //MX
      nidx0 = pcs_temp->GetNodeValueIndex(porosity_pcs_name_vector[i]);
      nidx1 = nidx0 + 1;
      if (mode == 0)                              // Gauss point values
      {
         assem->ComputeShapefct(1);
         primary_variable[i] = (1. - theta) * assem->interpolate(nidx0,
            pcs_temp) + theta * assem->interpolate(nidx1, pcs_temp);
      }                                           // Node values
      else if (mode == 1)
      {
         primary_variable[i] = (1. - theta) * pcs_temp->GetNodeValue(number,
            nidx0) + theta * pcs_temp->GetNodeValue(number, nidx1);
      }                                           // Element average value
      else if (mode == 2)
      {
         primary_variable[i] = (1. - theta) * assem->elemnt_average(nidx0,
            pcs_temp) + theta * assem->elemnt_average(nidx1, pcs_temp);
      }
      else if(mode==1)                            // Node values
      {
         primary_variable[i] = (1.-theta)*pcs_temp->GetNodeValue(number,nidx0) \
            + theta*pcs_temp->GetNodeValue(number,nidx1);
   }
      else if(mode==2)                            // Element average value
      {
         primary_variable[i] = (1.-theta)*assem->elemnt_average(nidx0,pcs_temp)
            + theta*assem->elemnt_average(nidx1,pcs_temp);
      }
   }

   //----------------------------------------------------------------------
   // Material models
   /*
   if(permeability_stress_mode==3) //Barton-Bandis WW
   {
      int i;
      double w[3], TG = 0.0;
      if(assem->cpl_pcs)
         TG = assem->interpolate(assem->NodalValC1)+273.15;
      else
         TG = 296.0;
      CalStressPermeabilityFactor(w, TG);
      porosity = 0.0;
   for(i=0; i<geo_dimension; i++)
   {
   if(i==0)
   porosity = pow(18.0e+6*permeability_tensor[i]*w[i]/(c_coefficient[21+i]*c_coefficient[21+i]),
   1.0/(float)geo_dimension);
   else
   porosity *= pow(18.0e+6*permeability_tensor[i]*w[i]/(c_coefficient[21+i]*c_coefficient[21+i]),
   1.0/(float)geo_dimension);
   }
   //
   return porosity;
   }
   */
   switch (porosity_model)
   {
      case 0:                                     // n = f(x)
         porosity = GetCurveValue(fct_number, 0, primary_variable[0], &gueltig);
         break;
      case 1:                                     // n = const
         porosity = porosity_model_values[0];
         break;
      case 2:                                     // n = f(sigma_eff), Stress dependance
         porosity = PorosityEffectiveStress(number, primary_variable[0]);
         break;
      case 3:                                     // n = f(S), Free chemical swelling
         porosity = PorosityVolumetricFreeSwellingConstantIonicstrength(number,
            primary_variable[0], primary_variable[1]);
         break;
      case 4:                                     // n = f(S), Constrained chemical swelling
         porosity = PorosityEffectiveConstrainedSwellingConstantIonicStrength(
            number, primary_variable[0], primary_variable[1], &porosity_sw);
         break;
      case 5:                                     // n = f(S), Free chemical swelling, I const
         porosity = PorosityVolumetricFreeSwelling(number, primary_variable[0],
            primary_variable[1]);
         break;
      case 6:                                     // n = f(S), Constrained chemical swelling, I const
         porosity = PorosityEffectiveConstrainedSwelling(number,
            primary_variable[0], primary_variable[1], &porosity_sw);
         break;
      case 7:                                     // n = f(mean stress) WW
         gval = ele_value_dm[number];
         primary_variable[0] = -gval->MeanStress(assem->gp) / 3.0;
         porosity = GetCurveValue(porosity_curve, 0, primary_variable[0],
            &gueltig);
         break;
      case 10:
                                                  /* porosity change through dissolution/precipitation */
         porosity = PorosityVolumetricChemicalReaction(number);
         break;
      case 11:                                    // n = temp const, but spatially distributed CB
         //porosity = porosity_model_values[0];
         porosity = m_msh->ele_vector[number]->mat_vector(por_index);
         break;
#ifdef GEM_REACT
      case 15:
         porosity = porosity_model_values[0];     // default value as backup

         for (size_t i=0; i < pcs_vector.size(); i++)
         {
            //		if ((pcs_vector[i]->pcs_type_name.find("FLOW") != string::npos)) {
            if (isFlowProcess(pcs_vector[i]->getProcessType()))
            {
               idx=pcs_vector[i]->GetElementValueIndex ( "POROSITY" );
               porosity = pcs_vector[i]->GetElementValue(number, idx);
               if (porosity <1.e-6 || porosity > 1.0) {cout <<"Porosity: error getting porosity for model 15. porosity: " <<porosity << " at node "<< number << endl; porosity = porosity_model_values[0];}
            }
         }

         break;
#endif
#ifdef BRNS
      case 16:
         porosity = porosity_model_values[0];     // default value as backup
         if ( aktueller_zeitschritt > 1 )
         {
            for (size_t i=0; i < pcs_vector.size(); i++)
            {
               pcs_temp = pcs_vector[i];
               //	            if ( pcs_temp->pcs_type_name.compare("GROUNDWATER_FLOW") == 0 ||
               //                     pcs_temp->pcs_type_name.compare("LIQUID_FLOW") == 0         ) {
               if ( pcs_temp->getProcessType() == GROUNDWATER_FLOW
                  || pcs_temp->getProcessType() == LIQUID_FLOW)
               {
                  int idx;
                  idx=pcs_temp->GetElementValueIndex ( "POROSITY" );

                  porosity = pcs_temp->GetElementValue(number, idx);
                  if (porosity <1.e-6)
                     cout << "error for porosity1 " << porosity << " node "<< number << endl;
               }
            }
         }
         break;
#endif
      default:
         cout << "Unknown porosity model!" << endl;
         break;
   }
   return (porosity);
}


/*------------------------------------------------------------------------*/
/* MAT-MP geometric properties */
/* 3 porosity */
/*------------------------------------------------------------------------*/

/**************************************************************************
FEMLib-Method:
Task: Porosity calc function
  Case overview
  0 Curve function
  1 Constant Value
  2 Function of normal stress from Geomechnical model
  3 Free swelling
  4 Constraint swelling
Programing:
07/2004 OK C++ Implementation
08/2004	CMCD Re-written based on MATCalcPorosity
10/2004 MX Swelling processes
04/2007 WW Porosity by gauss stress value
last modification:
*************************************************************************/
                                                  //WW
double CMediumProperties::Porosity(CElement* assem)
{
   static int nidx0,nidx1;
   double primary_variable[PCS_NUMBER_MAX];
   int gueltig;
   double porosity_sw, theta;
   std::string str;
   ///
   ElementValue_DM* gval = NULL;

   //----------------------------------------------------------------------
   // Functional dependencies
   number = assem->GetElementIndex();
   CRFProcess *pcs_temp;

   const size_t no_pcs_names (porosity_pcs_name_vector.size());
   for(size_t i=0;i<no_pcs_names;i++)
   {
      str = porosity_pcs_name_vector[i];
      pcs_temp = PCSGet(str, true);               //MX
      theta=pcs_temp->m_num->ls_theta;            //WW
      nidx0 = pcs_temp->GetNodeValueIndex(porosity_pcs_name_vector[i]);
      nidx1 = nidx0+1;
      if(mode==0)                                 // Gauss point values
      {
         assem->ComputeShapefct(1);
         primary_variable[i] = (1.-theta)*assem->interpolate(nidx0,pcs_temp)
            + theta*assem->interpolate(nidx1,pcs_temp);
      }
      else if(mode==1)                            // Node values
      {
         primary_variable[i] = (1.-theta)*pcs_temp->GetNodeValue(number,nidx0) \
            + theta*pcs_temp->GetNodeValue(number,nidx1);
      }
      else if(mode==2)                            // Element average value
      {
         primary_variable[i] = (1.-theta)*assem->elemnt_average(nidx0,pcs_temp)
            + theta*assem->elemnt_average(nidx1,pcs_temp);
      }
   }
   //----------------------------------------------------------------------
   // Material models
   switch (porosity_model)
   {
      case 0:                                     // n = f(x)
         porosity = GetCurveValue(fct_number,0,primary_variable[0],&gueltig);
         break;
      case 1:                                     // n = const
         porosity = porosity_model_values[0];
         break;
      case 2:                                     // n = f(sigma_eff), Stress dependance
         porosity = PorosityEffectiveStress(number,primary_variable[0]);
         break;
      case 3:                                     // n = f(S), Free chemical swelling
         porosity = PorosityVolumetricFreeSwellingConstantIonicstrength(number,primary_variable[0],primary_variable[1]);
         break;
      case 4:                                     // n = f(S), Constrained chemical swelling
         porosity = PorosityEffectiveConstrainedSwellingConstantIonicStrength(number,primary_variable[0],primary_variable[1],&porosity_sw);
         break;
      case 5:                                     // n = f(S), Free chemical swelling, I const
         porosity = PorosityVolumetricFreeSwelling(number,primary_variable[0],primary_variable[1]);
         break;
      case 6:                                     // n = f(S), Constrained chemical swelling, I const
         porosity = PorosityEffectiveConstrainedSwelling(number,primary_variable[0],primary_variable[1],&porosity_sw);
         break;
      case 7:                                     // n = f(mean stress) WW
         gval = ele_value_dm[number];
         primary_variable[0] = -gval->MeanStress(assem->GetGPindex())/3.0;
         porosity = GetCurveValue(porosity_curve,0,primary_variable[0],&gueltig);
         break;
      case 10:
                                                  /* porosity change through dissolution/precipitation */
         porosity = PorosityVolumetricChemicalReaction(number);
         break;
      case 11:                                    // n = const, but spatially distributed CB
         porosity = porosity_model_values[0];
         break;
#ifdef GEM_REACT
      case 15:

         porosity = porosity_model_values[0];     // default value as backup

         //                CRFProcess* mf_pcs = NULL;
         for (size_t i=0; i < pcs_vector.size(); i++)
         {
            pcs_temp = pcs_vector[i];
            //			if ((pcs_temp->pcs_type_name.compare("GROUNDWATER_FLOW") == 0) || (pcs_temp->pcs_type_name.compare("RICHARDS_FLOW") == 0)||(pcs_temp->pcs_type_name.compare("MULTI_PHASE_FLOW") == 0))
            if ( (pcs_temp->getProcessType() == GROUNDWATER_FLOW)
               || (pcs_temp->getProcessType() == RICHARDS_FLOW)
                                                  // TF
               || (pcs_temp->getProcessType() == MULTI_PHASE_FLOW) )
            {
               int idx=pcs_temp->GetElementValueIndex ( "POROSITY" );
               porosity = pcs_temp->GetElementValue(number, idx);
               if (porosity <1.e-6 || porosity > 1.0)
               {
                  std:: cout <<" error getting porosity model 15 " <<porosity << " node "<< number << endl;
                  porosity = porosity_model_values[0];
               }
            }
         }
         // KG44: TODO!!!!!!!!!!!!! check the above  ***************
         break;
#endif

      default:
         DisplayMsgLn("Unknown porosity model!");
         break;
   }
   //----------------------------------------------------------------------
   return (porosity);
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 MX Implementation
last modification:
**************************************************************************/
double CMediumProperties::PorosityEffectiveConstrainedSwellingConstantIonicStrength(long index,
double saturation,double temperature, double *porosity_sw)
{
   porosity_sw = porosity_sw;                     // WW Remove this argument
   /*  Soil Properties */

   double mat_mp_m, beta;
   double porosity = 0.0;
   static double theta;
   double porosity_n, porosity_IL, d_porosity, \
      density_rock, fmon, porosity_max;
   double n_sw, det_n_sw;
   double old_value, new_value;
   static double porosity_IP0;
   static double S_0;
   /* Component Properties */
   static double  satu, epsilon;
   static double satu_0=0.20;
   static double porosity_min=0.05;
   static double ion_strength;
   static double F_const=96484.6, epsilon_0=8.854e-12;
   static double R=8.314510, psi=1.0;
   /* Fluid properies */
   int phase=1;

   //--------------------------------------------------------------------
   // MMP medium properties
   S_0 =      porosity_model_values[1];           // Specific surface area [m^2/g]
   fmon =     porosity_model_values[2];           // Anteil quelfaehige mineral [-]
   mat_mp_m = porosity_model_values[3];           // Schichtenanzahl eines quellfähigen Partikels [-]
   beta =     porosity_model_values[6];           // modifications coefficient (z.B. free swelling Weimar beta=3.0)
   //--------------------------------------------------------------------
   // MSP solid properties
   CSolidProperties *m_msp = NULL;
   // long group = ElGetElementGroupNumber(index);
   long group = m_pcs->m_msh->ele_vector[index]->GetPatchIndex();
   m_msp = msp_vector[group];
   density_rock  = fabs(m_msp->Density(1));
   //--------------------------------------------------------------------

   /* Component Properties */
   ion_strength = porosity_model_values[4];       /*Ionic strength [M]*/
   satu_0 = porosity_model_values[5];             /*Initial saturation, if the sample is homogenous [-]*/
   porosity_min = porosity_model_values[7];       /*minimal porosity after swelling compaction*/

   /* Field v0ariables */
   theta = m_pcs->m_num->ls_theta;
   phase=1;
   satu = saturation;                             /*only for fluid, phase=1*/

   /*-----------------------------------------------------------------------*/
   /* Interlayer Porositaet berechnen */
   if (temperature < 0.0 || temperature > 600.0)
      temperature = 298.0;                        //ToDo, MX
   epsilon =87.0+exp(-0.00456*(temperature-273.0));
   porosity_n = porosity_model_values[0];

   /* Maximal inter layer porosity */
   porosity_max=fmon * psi * mat_mp_m * S_0 * (density_rock * 1.0e3) \
      * sqrt(epsilon * epsilon_0 * R * temperature / (2000.0 * F_const * F_const * ion_strength ));
   d_porosity=porosity_max*(pow(satu, beta)-pow(satu_0, beta));

   /*-----------Interlayer porosity calculation------------------*/

   /*  porosity_IL = porosity_IL*satu; */
   porosity_IL  =porosity_max*pow(satu,beta);
   // ElSetElementVal(index,PCSGetELEValueIndex("POROSITY_IL"),porosity_IL);
   m_pcs->SetElementValue(index, m_pcs->GetElementValueIndex("POROSITY_IL")+1,porosity_IL);

   /*-----------Effective porosity calculation------------------*/

   porosity = porosity_n-d_porosity;

   if (porosity<porosity_min)
      porosity =porosity_min;

   // ElSetElementVal(index,PCSGetELEValueIndex("POROSITY"),porosity);
   m_pcs->SetElementValue(index, m_pcs->GetElementValueIndex("POROSITY")+1,porosity);

   /*-----------Swelling strain rate, xie, wang and Kolditz, (23)------------------*/

   porosity_IP0  =porosity_n - porosity_max * pow(satu_0,beta);
   n_sw = d_porosity-(porosity_IP0 - porosity_min);

   if (n_sw<=0.0)
   {
      n_sw = 0.0;
   }
   int indx_sw0 =  m_pcs->GetElementValueIndex("n_sw");
   int indx_sw1 = indx_sw0 + 1;

   if ( aktuelle_zeit== dt)
   {
      m_pcs->SetElementValue(index,indx_sw0,n_sw);
      m_pcs->SetElementValue(index,indx_sw1,n_sw);
   }
   if ( aktuelle_zeit >= 2*dt)
   {
      m_pcs->SetElementValue(index,indx_sw1,n_sw);
   }

   old_value =  m_pcs->GetElementValue(index,indx_sw0);
   new_value =  m_pcs->GetElementValue(index,indx_sw1);

   if ( index == 399 && aktuelle_zeit== 10*dt)
   {

      index = index;                              //only for debug
   }

   // change rate of swelling strain
   det_n_sw = new_value - old_value;
   if (det_n_sw<-1.0e-8)
      det_n_sw = det_n_sw;
   if ( aktuelle_zeit== dt)
      m_pcs->SetElementValue(index, m_pcs->GetElementValueIndex("n_sw_rate"),det_n_sw);
   m_pcs->SetElementValue(index, m_pcs->GetElementValueIndex("n_sw_rate")+1,det_n_sw);

   return porosity;
}


//------------------------------------------------------------------------
//4..TORTUOSITY
//------------------------------------------------------------------------

/* 4B.4.3 non-linear flow */

/**************************************************************************
 9 Storage

**************************************************************************
 11 Permeability
**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
08/2004 OK MMP implementation
           based on GetPermeabilityTensor by OK
10/2004 SB adapted to het_file
last modification:
10/2010 TF changed access to process type
**************************************************************************/
double* CMediumProperties::PermeabilityTensor(long index)
{
   static double tensor[9];
   int perm_index=0;

   int idx_k, idx_n;
   double /*k_old, n_old,*/ k_new, n_new, k_rel, n_rel;

   // HS: move the following loop into the "if ( permeability_tensor_type == 0 )" scope.----
   // this is not necessary for in-isotropic case;
   // if(permeability_model==2)
   //    for(perm_index=0;perm_index<(int)m_pcs->m_msh->mat_names_vector.size();perm_index++)
   //        if(m_pcs->m_msh->mat_names_vector[perm_index].compare("PERMEABILITY")==0)
   //              break;
   // end of comment out codes--------------------------------------------------------------

   // -------------------------------------------------------------------------------------------------------
   // Start of K-C relationship. This will write a value into tensor[0] first, the values depends on the
   // value in permability_model. for 3 and 4, it gets the values from K-C relationship.
   // this will only influence then case when permeability_tensor_type = 0
   // -------------------------------------------------------------------------------------------------------
   if ( permeability_tensor_type == 0 )
   {                                              // only when permeability is isotropic
      tensor[0] = permeability_tensor[0];
      if( permeability_model == 2 )
      {                                           // here get the initial permeability values from material perperty class;
         // get the index:-------------------------------------------------------------------
         for(perm_index=0;perm_index<(int)m_pcs->m_msh->mat_names_vector.size();perm_index++)
            if(m_pcs->m_msh->mat_names_vector[perm_index].compare("PERMEABILITY")==0)
               break;
         // end of getting the index---------------------------------------------------------

         tensor[0] = m_msh->ele_vector[index]->mat_vector(perm_index);
                                                  //CMCD
         int edx = m_pcs->GetElementValueIndex("PERMEABILITY");
                                                  //CMCD
         m_pcs->SetElementValue(index,edx,tensor[0]);
      }
      else if ( permeability_model == 3 )
      {                                           // HS: 11.2008, for K-C relationship
         // get indexes
         idx_k = m_pcs->GetElementValueIndex("PERMEABILITY");
         idx_n = m_pcs->GetElementValueIndex("POROSITY");

         // get values of k0, n0, and n.
         k_new = m_pcs->GetElementValue( index, idx_k + 1 );
         n_new = m_pcs->GetElementValue( index, idx_n + 1 );

         // if first time step, get the k_new from material class
         if ( aktueller_zeitschritt == 1)
         {                                        // for the first time step
            // get the permeability.
            KC_permeability_initial = k_new = tensor[0];
         }

         // save old permeability
         m_pcs->SetElementValue( index, idx_k, k_new  );

         // calculate new permeability
         k_new = CMediumProperties::KozenyCarman( KC_permeability_initial,
            KC_porosity_initial,
            n_new );

         // save new permeability
         m_pcs->SetElementValue( index, idx_k+1, k_new   );

         // now gives the newly calculated value to tensor[]
         tensor[0] = k_new ;
      }
      else if ( permeability_model == 4 )
      {                                           // HS: 11.2008, for K-C_normalized relationship
         // get indexes
         idx_k = m_pcs->GetElementValueIndex("PERMEABILITY");
         idx_n = m_pcs->GetElementValueIndex("POROSITY");

         // get values of k0, n0, and n.
         k_new = m_pcs->GetElementValue( index, idx_k + 1 );
         n_new = m_pcs->GetElementValue( index, idx_n + 1 );

         // if first time step, get the k_new from material class
         if ( aktueller_zeitschritt == 0)
         {                                        // for the first time step
            // get the permeability.
            KC_permeability_initial = k_new = tensor[0];
         }

         // save old permeability
         m_pcs->SetElementValue( index, idx_k, k_new  );

         // calculate new permeability
         k_new = CMediumProperties::KozenyCarman_normalized( KC_permeability_initial,
            KC_porosity_initial,
            n_new );

         // save new permeability
         m_pcs->SetElementValue( index, idx_k+1, k_new   );

         // now gives the newly calculated value to tensor[]
         tensor[0] = k_new ;
      }
      else if ( permeability_model == 5 )
      {                                           // HS: 11.2008, for Clement clogging model

         // if first time step, do nothing. otherwise,
         if ( aktueller_zeitschritt > 1 )
         {
            for (size_t i=0; i < pcs_vector.size() ; i++)
            {
               m_pcs_tmp = pcs_vector[i];
               //	                if ( m_pcs_tmp->pcs_type_name.compare("GROUNDWATER_FLOW") == 0 ||
               //                                   m_pcs_tmp->pcs_type_name.compare("LIQUID_FLOW") == 0)
               if ( m_pcs_tmp->getProcessType () == GROUNDWATER_FLOW ||
                  m_pcs_tmp->getProcessType () == LIQUID_FLOW)
                  break;
            }
            // get index
            idx_k = m_pcs_tmp->GetElementValueIndex("PERMEABILITY");
            idx_n = m_pcs_tmp->GetElementValueIndex("POROSITY");

            // get values of n.
            n_new = m_pcs_tmp->GetElementValue( index, idx_n + 1 );

            // calculate new permeability
            // k_rel(n) = n_rel^{19/6}
            // first relative porosity change
            n_rel = n_new / permeability_porosity_model_values[0];
            // then relative permeability change
            k_rel = pow(n_rel, 19.0/6.0);
            // finially permeability
            k_new = k_rel * permeability_porosity_model_values[1];
            // save new permeability
            m_pcs->SetElementValue( index, idx_k+ 1, k_new );

            // now gives the newly calculated value to tensor[]
            tensor[0] = k_new ;
         }
      }
      else if ( permeability_model == 6 )
      {                                           // HS: 11.2008, for Clement biomass colony clogging

         // if first time step, do nothing. otherwise,
         if ( aktueller_zeitschritt > 1 )
         {
            for (size_t i=0; i < pcs_vector.size() ; i++)
            {
               m_pcs_tmp = pcs_vector[i];
               //	                if ( m_pcs_tmp->pcs_type_name.compare("GROUNDWATER_FLOW") == 0 ||
               //                                   m_pcs_tmp->pcs_type_name.compare("LIQUID_FLOW") == 0) TF
               if ( m_pcs_tmp->getProcessType () == GROUNDWATER_FLOW ||
                  m_pcs_tmp->getProcessType () == LIQUID_FLOW)
                  break;
            }
            // get index
            idx_k = m_pcs_tmp->GetElementValueIndex("PERMEABILITY");
            idx_n = m_pcs_tmp->GetElementValueIndex("POROSITY");

            // get values of n.
            n_new = m_pcs_tmp->GetElementValue( index, idx_n + 1 );

            // calculate new permeability
            // k_rel(n) = n_rel^{19/6}
            // first relative porosity change
            n_rel = n_new / permeability_porosity_model_values[0];

            // then relative permeability change
            k_rel = permeability_porosity_model_values[2] * \
               pow( ( n_rel - permeability_porosity_model_values[0] )/(1 - permeability_porosity_model_values[0]) , 3.0) \
               + ( 1 - permeability_porosity_model_values[2] ) * \
               pow( ( n_rel - permeability_porosity_model_values[0] )/(1 - permeability_porosity_model_values[0]) , 2.0);
            // finially permeability
            k_new = k_rel * permeability_porosity_model_values[1];
            // save new permeability
            m_pcs->SetElementValue( index, idx_k+1, k_new  );

            // now gives the newly calculated value to tensor[]
            tensor[0] = k_new ;
         }

      }
      else if ( permeability_model == 7 )
      {                                           // HS: 11.2008, for Clement biofilm clogging
         // if first time step, do nothing. otherwise,
         if ( aktueller_zeitschritt > 1 )
         {
            for (size_t i=0; i < pcs_vector.size() ; i++)
            {
               m_pcs_tmp = pcs_vector[i];
               //	                if ( m_pcs_tmp->pcs_type_name.compare("GROUNDWATER_FLOW") == 0 ||
               //	                		m_pcs_tmp->pcs_type_name.compare("LIQUID_FLOW") == 0) TF
               if ( m_pcs_tmp->getProcessType() == GROUNDWATER_FLOW ||
                  m_pcs_tmp->getProcessType () == LIQUID_FLOW)
                  break;
            }
            // get index
            idx_k = m_pcs_tmp->GetElementValueIndex("PERMEABILITY");
            idx_n = m_pcs_tmp->GetElementValueIndex("POROSITY");

            // get values of n.
            n_new = m_pcs_tmp->GetElementValue( index, idx_n + 1 );

            // calculate new permeability
            // k_rel(n) = n_rel^{19/6}
            // first relative porosity change
            n_rel = n_new / permeability_porosity_model_values[0];

            // then relative permeability change
            k_rel = (pow( ( n_rel - permeability_porosity_model_values[0] )/(1 - permeability_porosity_model_values[0]) , permeability_porosity_model_values[2] ) \
               + permeability_porosity_model_values[3] ) / \
               ( 1 + permeability_porosity_model_values[3]) ;
            // finially permeability
            k_new = k_rel * permeability_porosity_model_values[1];
            // save new permeability
            m_pcs->SetElementValue( index, idx_k+1, k_new  );

            // now gives the newly calculated value to tensor[]
            tensor[0] = k_new ;
         }
      }
   }
   // end of K-C relationship-----------------------------------------------------------------------------------

   switch(geo_dimension)
   {
      case 1:                                     // 1-D
         // HS: tensor[0] already set, doing nothing here;
         break;
      case 2:                                     // 2-D
         if(permeability_tensor_type==0)
         {
            // tensor[0] = permeability_tensor[0]; // HS: done already;
            tensor[1] = 0.0;
            tensor[2] = 0.0;
            // tensor[3] = permeability_tensor[0];
            tensor[3] = tensor[0];                // HS: use the existing value;

            // HS: this is not needed any more--------------------------------
            // if(permeability_model==2) {
            //SB 4218	tensor[0] = GetHetValue(index,"permeability");
            // 	tensor[0] = m_msh->ele_vector[index]->mat_vector(perm_index);
            //	tensor[3] = tensor[0];
            // }
            // end of comment out section-------------------------------------
         }
         else if(permeability_tensor_type==1)
         {
            tensor[0] = permeability_tensor[0];
            tensor[1] = 0.0;
            tensor[2] = 0.0;
            tensor[3] = permeability_tensor[1];
         }
         else if(permeability_tensor_type==2)
         {
            tensor[0] = permeability_tensor[0];
            tensor[1] = permeability_tensor[1];
            tensor[2] = permeability_tensor[2];
            tensor[3] = permeability_tensor[3];
         }
         break;
      case 3:                                     // 3-D
         if(permeability_tensor_type==0)
         {
            // tensor[0] = permeability_tensor[0]; // HS: not needed. already done before
            tensor[1] = 0.0;
            tensor[2] = 0.0;
            tensor[3] = 0.0;
            // tensor[4] = permeability_tensor[0]; // HS: using the line below instead;
            tensor[4] = tensor[0];                // using the existing value
            tensor[5] = 0.0;
            tensor[6] = 0.0;
            tensor[7] = 0.0;
            // tensor[8] = permeability_tensor[0]; // HS: using the line below instead;
            tensor[8] = tensor[0];                // using the existing value
            // HS: this is not needed any more--------------------------------
            // if(permeability_model==2) {
            // SB 4218	tensor[0] = GetHetValue(index,"permeability");
            // 	tensor[0] = m_pcs->m_msh->ele_vector[index]->mat_vector(perm_index);
            // 	tensor[4] = tensor[0];
            // 	tensor[8] = tensor[0];
            // }
            // end of comment out section-------------------------------------
         }
         else if(permeability_tensor_type==1)
         {
            tensor[0] = permeability_tensor[0];
            tensor[1] = 0.0;
            tensor[2] = 0.0;
            tensor[3] = 0.0;
            tensor[4] = permeability_tensor[1];
            tensor[5] = 0.0;
            tensor[6] = 0.0;
            tensor[7] = 0.0;
            tensor[8] = permeability_tensor[2];
         }
         else if(permeability_tensor_type==2)
         {
            tensor[0] = permeability_tensor[0];
            tensor[1] = permeability_tensor[1];
            tensor[2] = permeability_tensor[2];
            tensor[3] = permeability_tensor[3];
            tensor[4] = permeability_tensor[4];
            tensor[5] = permeability_tensor[5];
            tensor[6] = permeability_tensor[6];
            tensor[7] = permeability_tensor[7];
            tensor[8] = permeability_tensor[8];
         }
         break;
   }
#ifdef RFW_FRACTURE
   double Krel = RelativePermeability(index);     // rfw cmcd

   for (int i = 0; i< geo_dimension*geo_dimension;i++)
      tensor[i]*=Krel;
#endif
   return tensor;
}


#ifdef RFW_FRACTURE
                                                  //rfw cmcd
double CMediumProperties::RelativePermeability (long index)
{
   double Krel = 1.0;
   string name;
   int size = (int)relative_permeability_function.size();
   for(int i=0; i<size; ++i)
   {
      string name = relative_permeability_function[i];
      //if (name == "PERMEABILITY_FUNCTION_DEFORMATION") Krel *= Call to function not yet in existence;
      //if (name == "PERMEABILITY_FUNCTION_PRESSURE") Krel *= PermeabilityPressureFunction(index,double *gp,double theta);
      //NOTE: the function PermeabilitySaturationFunction is overloaded, this must be dealt with somehow
      //if (name == "PERMEABILITY_SATURATION") Krel *= PermeabilitySaturationFunction(long number,double*gp,double theta, int phase);
      //if (name == "PERMEABILITY_FUNCTION_STRESS") Krel *= Call to function not yet in existence;
      //if (name == "PERMEABILITY_FUNCTION_VELOCITY") Krel *= Call to function not yet in existence;
      //if (name == "PERMEABILITY_FUNCTION_POROSITY") Krel *= PermeabilityPorosityFunction(index,double *gp,double theta);
#ifdef RFW_FRACTURE
      if (name == "PERMEABILITY_FRAC_APERTURE") Krel *= PermeabilityFracAperture(index);
#endif
   }
   return Krel;
}
#endif
//------------------------------------------------------------------------
//12.(i) PERMEABILITY_FUNCTION_DEFORMATION
//------------------------------------------------------------------------
//
//------------------------------------------------------------------------
//12.(ii) PERMEABILITY_FUNCTION_PRESSURE
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//12.(i) PERMEABILITY_FUNCTION_SATURATION
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//12.(vi) PERMEABILITY_FUNCTION_POROSITY
//------------------------------------------------------------------------

//---------------------------------------------------------------------------------
//12.(vii) PERMEABILITY_FUNCTION_FRAC_APERTURE
//RFW 07/2005
//---------------------------------------------------------------------------------
#ifdef RFW_FRACTURE
double CMediumProperties::PermeabilityFracAperture(long index)
{

   CElem *elem = NULL;
   CElem *elem2 = NULL;
   CFEMesh* m_msh = NULL;
   CSolidProperties* mat_pointer;
   CMediumProperties *m_mmp = NULL;
   double ApertureSum=0,ApertureAvg=0, normalised_perm, min_aperture,
      permeability, c =0./*c is the fraction of fracture that is closed*/, WeightSum=0;
   double Sum_Squared_Diffs=0, Std_Dev =0, closed=0, total_count=0;
   vector<double> aperture_list;
   int i;

   m_msh = fem_msh_vector[0];                     // Something must be done later on here.
   elem = m_msh->ele_vector[index];
   int group = elem->GetPatchIndex();
   m_mmp = mmp_vector[group];

   double mini=100, roughness_corr=1, tortuosity_corr=1;

   if(elem->InFrac())
   {
      if(!elem->PermeabilityIsSet())              //if 1
      {
         mat_pointer = msp_vector[group];
         min_aperture = mat_pointer->Get_Youngs_Min_Aperture(elem);

         if(m_mmp->fracs_set!=1)                  //if 2
         {
            m_mmp->fracs_set=1;
            switch(frac_perm_average_type[0])
            {
               //******************************************************************************************************
               case 'A':                          // -------------------arithmetic average----------------------
                  for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
                  {
                     elem2 = m_msh->ele_vector[i];
                                                  //are they part of the same fracture?
                     if(elem->GetFracNum() == elem2->GetFracNum() )
                     {
                                                  //change this to simply calling CalculateFracAperture, it will check if aperture has been set
                        if (elem2->ApertureIsSet())
                        {

                           if( (elem2->GetAperture()-min_aperture) >= 0.00000001)
                           {
                              ApertureSum += (elem2->GetAperture()-min_aperture)*elem2->GetWeight();
                              aperture_list.push_back( (elem2->GetAperture()-min_aperture) );
                              WeightSum += elem2->GetWeight();
                           }
                           else
                           {
                              closed = closed + elem2->GetWeight();
                           }
                           total_count = total_count + elem2->GetWeight();
                        }
                     }
                  }

                  ApertureAvg = (ApertureSum / WeightSum);

                  //roughness correction, after Zimmerman and Bovarsson 1996
                  if(roughness[0]=='C' || roughness[0]=='c')
                  {
                     for(i=0; i!=(int)aperture_list.size(); ++i)
                     {
                        Sum_Squared_Diffs += pow((aperture_list[i] - ApertureAvg),2);
                     }
                     Std_Dev = pow((Sum_Squared_Diffs/(double)aperture_list.size()),0.5);
                  }
                  else Std_Dev = 0;

                  //calculating the peremability, permeability is the average fracture perm, normalised_perm is the perm of the current frac element
                  c = closed/total_count;         //closed fractioncorrection afer Z and B 1996, and Walsh 1981
                  if(index==1)
                     cout<<"\nTotal weight = "<< total_count<<endl;

                  tortuosity_corr = (1-c)/(1+c);
                  if(Std_Dev < 0.73*ApertureAvg)
                  {
                     roughness_corr = 1 - 1.5*pow(Std_Dev,2.0)/pow(ApertureAvg, 2.0);
                  }
                  else
                  {
                     roughness_corr = 0.2;
                  }
                  break;
                  //******************************************************************************************************
               case 'G':                          // -------------------geometric average----------------------
                  for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
                  {
                     elem2 = m_msh->ele_vector[i];
                                                  //are they part of the same fracture?
                     if(elem->GetFracNum() == elem2->GetFracNum() )
                     {
                        if (elem2->ApertureIsSet())
                        {
                           if( (elem2->GetAperture()-min_aperture) > 0.0000001)
                           {
                              ApertureSum = ApertureSum + elem2->GetWeight()*log(elem2->GetAperture()-min_aperture);
                              aperture_list.push_back( (elem2->GetAperture()-min_aperture));
                              WeightSum += elem2->GetWeight();
                           }
                           else
                              closed = closed + elem2->GetWeight();

                           total_count = total_count + elem2->GetWeight();
                        }
                     }
                  }
                  ApertureAvg = exp(ApertureSum/WeightSum);

                  //roughness correction, after Zimmerman and Bovarsson 1996
                  //This roughness correction may not really be appropriate here, as this is and arithmetic standard deviation; I'm
                  //not very confident about the math here.
                  if(roughness[0]=='C' || roughness[0]=='c')
                  {
                     for(i=0; i!=(int)aperture_list.size(); ++i)
                     {
                        Sum_Squared_Diffs += pow((aperture_list[i] - ApertureAvg),2);
                     }
                     Std_Dev = pow((Sum_Squared_Diffs/(double)aperture_list.size()),0.5);
                  }

                  //calculating the peremability, permeability is the average fracture perm, normalised_perm is the perm of the current frac element
                  c = closed/total_count;         //closed fraction
                  if(Std_Dev < 0.73*ApertureAvg)
                  {
                     roughness_corr = 1 - 1.5*pow(Std_Dev,2.0)/pow(ApertureAvg, 2.0);
                  }
                  else
                  {
                     roughness_corr = 0.2;
                  }
                  break;
                  //******************************************************************************************************
               case 'H':                          // -------------------harmonic average----------------------
                  for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
                  {
                     elem2 = m_msh->ele_vector[i];
                     if(elem->GetFracNum() == elem2->GetFracNum() )
                     {
                        if (elem2->ApertureIsSet())
                        {
                           if( (elem2->GetAperture()-min_aperture) > 0.0000001)
                           {
                              ApertureSum = ApertureSum + ( elem2->GetWeight()/((elem2->GetAperture()-min_aperture)) );
                              WeightSum += elem2->GetWeight();
                           }
                           else
                              closed = closed + elem2->GetWeight();

                           total_count = total_count + elem2->GetWeight();
                        }
                     }
                  }

                  ApertureAvg = WeightSum / ApertureSum;

                  //no roughness correction is available here, the harmonic mean is already an extreme case
                                                  //error message should print once
                  if( (roughness[0]=='C' || roughness[0]=='c') && index==0)
                  {
                     cout<<"Error in CMediumProperties:: Roughness correction not implemented for harmonic mean.\n";
                  }

                  //calculating the peremability, permeability is the average fracture perm, normalised_perm is the perm of the current frac element
                  c = closed/total_count;         //closed fraction, also not used for perm calculation with harmonic mean

                  break;
                  //******************************************************************************************************
               default:                           // ---------------------------error------------------------------
                  cout << "Error in CMediumProperties::PermeabilityFracAperture.  No averaging type is set.\n";
                  abort();
                  break;
                  //******************************************************************************************************
            }                                     //end of switch case
            //the values are stored in the mmp

            permeability = (  pow(ApertureAvg, 2.0) / 12  ) * tortuosity_corr * roughness_corr;
            normalised_perm = permeability * ApertureAvg / (elem->GetAperture());

            m_mmp->frac_perm[elem->GetFracNum()] = permeability;
            m_mmp->avg_aperture[elem->GetFracNum()] = ApertureAvg;
            m_mmp->closed_fraction[elem->GetFracNum()] = c;
         }                                        //end if 2
         else
         {
            permeability = m_mmp->frac_perm[elem->GetFracNum()];
            ApertureAvg = m_mmp->avg_aperture[elem->GetFracNum()];

            normalised_perm = permeability * ApertureAvg / (elem->GetAperture());
         }

         //the values stored in elem are used for flow calculation
         elem->SetPermeability(normalised_perm);

         if(index==0)
         {
            cout <<"\n************************************************************\n";
            cout <<"Some numbers from PermeabilityFracAperture:\n";
            cout << " \nAvg Aperture =\t"<<ApertureAvg<<"\nStdDev Aperture =\t"<<Std_Dev<<"   \nAvg Permeability =\t" << permeability <<"\n";
            cout << " \nClosure ratio =\t"<<c<<" \nRoughness_corr =\t"<<roughness_corr<<"   \nTortuosity_corr =\t" << tortuosity_corr <<"\n";
            cout <<"\n************************************************************\n";
         }
      }                                           //end if 1
      else                                        // if the permeability has aleady been set for this timestep
      {
         normalised_perm = elem->GetPermeability();
      }
      return normalised_perm;
   }
   else                                           // RFW 18/11/2005
   {
      normalised_perm = 1e-9;                     //this is a kind of default value but will rarely be used, need to fix this in the future
      elem->SetPermeability(normalised_perm);
      return normalised_perm;
   }
}


/*************************************************************************
 ROCKFLOW - Funktion: CSolidProperties::CalculateFracAperture
 Aufgabe:
   Calculate the aperture of the fracture at a given element, currently
   aperture size is only in y-direction, and only for 2D triangles.
 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: const long index   :  element index
   E: const double dy    :  step size for search for frac boundary
 Ergebnis:
   Returns the aperture as a double
 Programmaenderungen:
04/2005 RFW Implementierung
*************************************************************************/
double CMediumProperties::CalculateFracAperture(CElem* elem, double delta_y)
{
   // TEST *************
   //cout << "ELEMENT: "<<index <<endl;
   // TEST *************

   vector<double> intercept;
   double aperture=0, dy, dx;
   vector<double> centroid;
   CElem *last_elem=NULL, *current_elem=NULL;

   vec<CElem*> neighbor_last, neighbor_current, neighbor_neighbor_last;
   //the next line is here because Wenqings vec class destructor does not like empty vec's, and there is no guarantee these vecs will be filled in any given call to this function
   neighbor_last.resize(1); neighbor_current.resize(1); neighbor_neighbor_last.resize(1);
   long num_face1, num_face2;
   bool at_frac_bound;
   long step, current_index, last_index=elem->GetIndex(), group, index = elem->GetIndex();
   vector<CElem*> neighbor_vec;
   bool neighbors, shared_neighbors, neighbors_neighbor;
   CFEMesh* msh_pointer = NULL;
   msh_pointer = fem_msh_vector[0];               // Something must be done later on here.

   if(elem->InFrac())                             // RFW 18/11/2005
   {
      if(!elem->ApertureIsSet())                  //if the aperture has not already been set for this timestep
      {
         elem->CalcDispGravityCenter(centroid);

         for(int k=-1; k!=3; k=k+2  )             //start big for
         {
            //assigning search directions
            dx = elem->GetFracDx()*delta_y*k;
            dy = elem->GetFracDy()*delta_y*k;

            at_frac_bound = false;
            current_index=index;
            step = 1;
            group=0;
            neighbors=false;  shared_neighbors=false;  neighbors_neighbor=false;
            at_frac_bound=false;
            neighbor_vec.clear();

            while(!at_frac_bound && step<10000)   //second condition to avoid infinite loop
            {

               last_index = current_index;

               current_index = MSHWhatElemIsPointIn(centroid[0]+(step*dx), centroid[1]+(step*dy), current_index);
               step = step + 1;
               //************************
               if(current_index >= 0)
               {
                  if(msh_pointer->ele_vector[current_index]->InFrac() && current_index >=0)
                     { }
                     else if(current_index >= 0)  //else 1, current element is outside fracture
                  {
                     at_frac_bound = true;

                     if(last_index == current_index)
                        cout << "PROBLEM IN: CalculateFracAperture, ~line 5280 of rf_mmp_new.\n";

                     last_elem = msh_pointer->ele_vector[last_index];
                     num_face1 = last_elem->GetFacesNumber();
                     neighbor_last.resize(num_face1);
                     last_elem->GetNeighbors(neighbor_last);

                     current_elem = msh_pointer->ele_vector[current_index];
                     num_face2 = current_elem->GetFacesNumber();
                     neighbor_current.resize(num_face2);
                     current_elem->GetNeighbors(neighbor_current);

                     for(long i =0; i!=num_face1; ++i)
                     {
                        if(neighbor_last[i] ==current_elem)
                        {
                           neighbors = true;
                           neighbor_vec.push_back(last_elem);
                           neighbor_vec.push_back(current_elem);
                        }
                     }

                     if(!neighbors)
                     {                            //else 2
                        for(long i =0; i!=num_face1; ++i)
                        {
                           for(long j =0; j!=num_face2; ++j)
                           {
                              //one of the neighbours of the elements borders one of the neighbors of the other
                              if( neighbor_last[i]==neighbor_current[j] && (neighbor_last[i]!=NULL && neighbor_current[j]!=NULL) )
                              {
                                 shared_neighbors = true;

                                 if(neighbor_current[j]->InFrac())
                                 {
                                    neighbor_vec.push_back(neighbor_current[j]);
                                    neighbor_vec.push_back(current_elem);
                                 }
                                 else
                                 {
                                    neighbor_vec.push_back(last_elem);
                                    neighbor_vec.push_back(neighbor_current[j]);
                                 }
                              }
                           }
                        }
                     }

                     if(!neighbors)
                     {
                        bool match = false;
                        vector<CNode*> matching_nodes;
                        for(long i =0; i!=num_face1; ++i)
                        {
                           for(long j =0; j!=num_face2; ++j)
                           {
                              if(neighbor_last[i]!=NULL && neighbor_current[j]!=NULL)
                              {
                                 match = MSHGetCommonNodes(neighbor_last[i], neighbor_current[j], matching_nodes);
                                 //one of the neighbours of the elements borders one of the neighbors of the other
                                 if(match && matching_nodes.size()>1 )
                                 {
                                    neighbors_neighbor = true;
                                    if( neighbor_current[j]->InFrac() && neighbor_last[i]->InFrac())
                                    {
                                       neighbor_vec.push_back(neighbor_current[j]);
                                       neighbor_vec.push_back(current_elem);
                                    }
                                    else if( !neighbor_current[j]->InFrac() && !neighbor_last[i]->InFrac())
                                    {
                                       neighbor_vec.push_back(last_elem);
                                       neighbor_vec.push_back(neighbor_last[i]);
                                    }
                                    else if( !neighbor_current[j]->InFrac() && neighbor_last[i]->InFrac() )
                                    {
                                       neighbor_vec.push_back(neighbor_last[i]);
                                       neighbor_vec.push_back(neighbor_current[j]);
                                    }
                                    else
                                       cout << "ERROR 4 in CalculateFracAperture, element: "<<index<<endl ;
                                 }
                              }
                           }
                        }
                     }

                     vector<double> xnode, ynode, xstep, ystep;
                     xnode.resize(2); ynode.resize(2);
                     xstep.resize(2); ystep.resize(2);
                     if(neighbors || shared_neighbors || neighbors_neighbor)
                     {
                        bool good_intersect = false, match = false;
                        CElem *neighbor_0, *neighbor_1;
                        neighbor_0 = neighbor_vec[0];
                        neighbor_1 = neighbor_vec[1];
                        vector<CNode*> matching_nodes;
                        for(long r=0; r!=(long)neighbor_vec.size(); r=r+2)
                        {
                           match = MSHGetCommonNodes(neighbor_vec[r], neighbor_vec[r+1], matching_nodes);
                           if(match && !good_intersect)
                           {
                              for(int m=0; m!=(int)matching_nodes.size(); ++m)
                              {
                                 //in order to let this function run when there is no deformation
                                 //IMPORTANT NOTE: this has not been tested properly
                                 for(int i=0;i<(int)pcs_vector.size();i++)
                                 {
                                    m_pcs = pcs_vector[i];
                                    if(m_pcs->pcs_type_name.find("DEFORMATION")!=string::npos)
                                    {
                                       xnode[m] = matching_nodes[m]->X_displaced();
                                       ynode[m] = matching_nodes[m]->Y_displaced();
                                    }
                                    else
                                    {
                                       xnode[m] = matching_nodes[m]->X();
                                       ynode[m] = matching_nodes[m]->Y();
                                    }
                                 }
                              }
                                                  //2 steps back
                              xstep[0] = centroid[0]+(step*dx)-3*dx;
                                                  //2 steps forward
                              xstep[1] = centroid[0]+(step*dx);
                                                  //2 steps back
                              ystep[0] = centroid[1]+(step*dy)-3*dy;
                                                  //2 steps forward
                              ystep[1] = centroid[1]+(step*dy);
                              good_intersect = LineSegmentIntersection(xnode, ynode, xstep, ystep, intercept);
                           }
                        }
                        if(!good_intersect)
                        {
                           cout<< "ERROR2 in CalculateFracAperture, element: " <<index<<endl;
                           intercept.push_back( centroid[0]+(step*dx)-1.5*dx );
                           intercept.push_back( centroid[1]+(step*dy)-1.5*dy );
                        }
                     }
                     else
                     {                            //else 3

                        intercept.push_back( centroid[0]+(step*dx)-1.5*dx );
                        intercept.push_back( centroid[1]+(step*dy)-1.5*dy );
                        cout<< "ERROR3 in CalculateFracAperture, element: " <<index<<".  INSTABILITY??"<<endl;

                     }                            //end else 3
                     //  //end else 2

                  }                               //end else 1
               }
               else                               // else 4, point has jumped outside of model boundaries
               {
                  at_frac_bound = true;
                  intercept.push_back( centroid[0]+(step*dx)-1.5*dx );
                  intercept.push_back( centroid[1]+(step*dy)-1.5*dy );
                  cout<< "ERROR5 in CalculateFracAperture, element: " <<index<<"  current_index: "<<current_index<<endl;
                  cout<< "  Intersection with model boundary."<<endl;
               }                                  //end else 4
            }                                     // end while
         }                                        //end big for
         double xtop, xbot, ytop,ybot;
         xtop = intercept[0];
         xbot = intercept[2];
         ytop = intercept[1];
         ybot = intercept[3];

         aperture = sqrt(  pow( (intercept[3]-intercept[1]),2 )  +  pow( (intercept[2]-intercept[0]),2 )  );

         elem->SetAperture(aperture);

         return aperture;
      }
      else                                        // if the aperture has already been set for this timestep
      {
         aperture = elem->GetAperture();
         return aperture;
      }
   }
   else
   {
      cout << "Element "<<index<<"neglected from fracture calculations.\n";
      aperture = -10;
      return aperture;
   }
}
#endif
//------------------------------------------------------------------------
//13. CAPILLARY_PRESSURE_FUNCTION
//------------------------------------------------------------------------
/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/1997 CT Erste Version
09/1998 CT Kurven ueber Schluesselwort #CURVES
08/1999 CT Erweiterung auf n-Phasen begonnen
08/2004 OK Template for MMP implementation
01/2004 OK Mode 3, given saturation
05/2008 WW Brook & Corey
last modification:
ToDo: GetSoilRelPermSatu
**************************************************************************/
double CMediumProperties::CapillaryPressureFunction(long number,double*gp,double theta,\
int phase,double saturation)
{
   //OK411
   theta = theta;
   gp = gp;
   number = number;

   double density_fluid = 0.;
   double van_saturation,van_beta;
   CFluidProperties *m_mfp = NULL;
   // only for liquid phase
   //if (phase == 0)  //YD
   //  return 0.;
   // static int nidx0,nidx1;
   int gueltig;
   double capillary_pressure=0.0;

   //----------------------------------------------------------------------
   switch(capillary_pressure_model)
   {
      case 0:                                     // k = f(x) user-defined function
         capillary_pressure = GetCurveValue((int)capillary_pressure_model_values[0],0,saturation,&gueltig);
         break;
      case 1:                                     // constant
         capillary_pressure = capillary_pressure;
         break;
      case 2:                                     // Lineare Kapillardruck-Saettigungs-Beziehung
         // kap12 steigt linear von 0 auf kap[2] im Bereich satu_water_saturated bis satu_water_residual
         break;
      case 3:                                     // Parabolische Kapillardruck-Saettigungs-Beziehung
         // kap12 steigt parabolisch von 0 auf kap[2] im Bereich satu_water_saturated bis satu_water_residual mit Exponent mat.kap[3]
         break;
      case 4:                                     // Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 Page 894 Equation 21
         // YD Advances in Water Resource 19 (1995) 25-38
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         m_mfp = MFPGet("LIQUID");                //YD
         density_fluid = m_mfp->Density();
         van_beta = 1./(1-saturation_exp[phase]);
         van_saturation =(saturation - saturation_res[phase])/(saturation_max[phase]-saturation_res[phase]);
         capillary_pressure =density_fluid*gravity_constant/capillary_pressure_model_values[0] \
            *pow(pow(van_saturation,-1./saturation_exp[phase])-1.,1./van_beta);
         /*
         //WW
         if (capillary_pressure > (9.0e4 - MKleinsteZahl))
              capillary_pressure = 9.0e4;  // Mehr als Vollsaettigung mit Wasser
         */
         if (capillary_pressure < (0. + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
               capillary_pressure = 0. + MKleinsteZahl;
         //if (capillary_pressure < 200)
         //cout<<"c"<<capillary_pressure;

         break;
      case 5:                                     // Haverkamp Problem: aus SOIL SCI. SOC. AM. J., VOL. 41, 1977, Pages 285ff
         break;
      case 6:                                     // Brooks & Corey. 22.05.2008. WW
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
            saturation = saturation_res[phase] + MKleinsteZahl;
         van_saturation =(saturation - saturation_res[0])/(1.0-saturation_res[0]-saturation_res[1]);
         if(van_saturation<1.0e-9)
            van_saturation = 1.0e-9;
         capillary_pressure =capillary_pressure_model_values[0]*pow(van_saturation, -1.0/saturation_exp[0]);
#ifdef CUTOFF
         // PCH putting Pmax constraint as in TOUGH2
         if(capillary_pressure > PMAX)
            capillary_pressure = PMAX;
#endif
         break;
      case 7:                                     // Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 Page 894 Equation 21
         break;
      case 8:                                     // 3Phasen ueber Kurven
         break;
      case 16:                                    // JOD
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         m_mfp = MFPGet("LIQUID");                //YD
         density_fluid = m_mfp->Density();
         van_beta = 1./(1-saturation_exp[phase]);
         van_saturation =(saturation - saturation_res[phase])/(saturation_max[phase]-saturation_res[phase]);
         capillary_pressure =density_fluid*gravity_constant/saturation_alpha[phase] \
            *pow(pow(van_saturation,-1./saturation_exp[phase])-1.,1./van_beta);

         if (capillary_pressure < (0. + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
               capillary_pressure = 0. + MKleinsteZahl;
         break;
      default:
         cout << "Error in CFluidProperties::CapillaryPressure: no valid material model" << endl;
         break;
   }
   return capillary_pressure;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2002 CT CECalcSatuFromCap
        SB Extensions case 4 and case 5
02/2005 OK CMediumProperties function
08/2005 WW Remove interploation
Last modified:
**************************************************************************/
double CMediumProperties::SaturationCapillaryPressureFunction
(const double capillary_pressure, const int phase)
{
   static double saturation;
   int gueltig;
   double van_beta = 0.0;
   double density_fluid = 0.;
   CFluidProperties *m_mfp = NULL;
   //----------------------------------------------------------------------
   switch(capillary_pressure_model)
   {
      case 0:                                     // k = f(x) user-defined function
         saturation = GetCurveValueInverse((int)capillary_pressure_model_values[0],0,capillary_pressure,&gueltig);
         //WW 07.07.2008; changed 27.10.2010. BG + SB, version 5.0.08
         if (saturation > (saturation_max[0] - MKleinsteZahl))
            saturation = saturation_max[0] - MKleinsteZahl;
         if (saturation < (saturation_res[0] + MKleinsteZahl))
            saturation = saturation_res[0] + MKleinsteZahl;
         break;
      case 1:                                     // constant
         saturation = 1.0;                        //MX test for DECOVALEX
         break;
      case 2:                                     // Lineare Kapillardruck-Saettigungs-Beziehung
         // kap12 steigt linear von 0 auf kap[2] im Bereich satu_water_saturated bis satu_water_residual
         break;
      case 3:                                     // Parabolische Kapillardruck-Saettigungs-Beziehung
         // kap12 steigt parabolisch von 0 auf kap[2] im Bereich satu_water_saturated bis satu_water_residual mit Exponent mat.kap[3]
         break;
      case 4:                                     // Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 Page 894 Equation 21
         m_mfp = mfp_vector[phase];
         density_fluid = m_mfp->Density();
         van_beta = 1/(1-saturation_exp[phase]) ;
         if(capillary_pressure < MKleinsteZahl)
            saturation = saturation_max[phase];
         else
            saturation = pow((pow(capillary_pressure*capillary_pressure_model_values[0]\
               /density_fluid/gravity_constant,van_beta)+1),-1.*saturation_exp[phase]) \
               *(saturation_max[phase]-saturation_res[phase]) + saturation_res[phase];
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;

         /*OK
               a = kap[2];             // Test 0.005
               m = kap[3];             // Test 0.5
               n = kap[4];             // Test 2.
              saturation = pow(1.0 + pow(a*cap[1],n),-m) * (kap[1] - kap[0]) + kap[0];
         */
         break;
      case 5:                                     // Haverkamp Problem: aus SOIL SCI. SOC. AM. J., VOL. 41, 1977, Pages 285ff
         /*
               a = kap[2];             // 1.611e6 fuer Haverkamp
               b = kap[3];             // 3.96 fuer Haverkamp
               // SB: Umrechnen von cap in p_cap [cm]
               p_cap = cap[1] / GetFluidDensity(1, element, 0., 0., 0., 0.) / gravity_constant* 100.0;
               satu_return[1] = (a * (kap[1] - kap[0]))/(a + pow(p_cap,b)) + kap[0];
         */
         break;
      case 6:                                     // Brooks & Corey. 22.05.2008. WW
         if(capillary_pressure < MKleinsteZahl)
            saturation = saturation_max[phase];
         else
            saturation = (1.0-saturation_res[0]-saturation_res[1])
               *pow(capillary_pressure/capillary_pressure_model_values[0], -saturation_exp[0])
               + saturation_res[0];
         if (saturation > (saturation_max[0] - MKleinsteZahl))
            saturation = saturation_max[0] - MKleinsteZahl;
         if (saturation < (saturation_res[0] + MKleinsteZahl))
            saturation = saturation_res[0] + MKleinsteZahl;
         break;
      case 7:                                     // Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 Page 894 Equation 21
         break;
      case 8:                                     // 3Phasen ueber Kurven
         break;
      case 16:                                    // VanGenuchten/Mualem for light oil: R.B. Thoms MastersThesis 2003    JOD
         m_mfp = mfp_vector[phase];
         density_fluid = m_mfp->Density();
         van_beta = 1/(1-saturation_exp[phase]) ;
         if(capillary_pressure < MKleinsteZahl)
            saturation = saturation_max[phase];
         else
            saturation = pow((pow(capillary_pressure * saturation_alpha[phase]\
               /density_fluid/gravity_constant,van_beta)+1),-1.*saturation_exp[phase]) \
               *(saturation_max[phase]-saturation_res[phase]) + saturation_res[phase];
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         break;
      default:
         cout << "Error in CMediumProperties::SaturationCapillaryPressureFunction: no valid material model" << endl;
         break;
   }
   return saturation;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2002 CT CECalcSatuFromCap
        SB Extensions case 4 and case 5
02/2005 OK CMediumProperties function
05/2008 WW Brooks & Corey.
Last modified:
**************************************************************************/
double CMediumProperties::SaturationCapillaryPressureFunction(long number,double*gp,double theta,int phase)
{
   //OK411
   theta = theta;
   gp = gp;
   number = number;

   // static int nidx0,nidx1;
   static double saturation;
   int gueltig;
   double capillary_pressure=0.0;
   double van_beta = 0.0;
   double density_fluid = 0.;
   CFluidProperties *m_mfp = NULL;
   //---------------------------------------------------------------------
   if(mode==2)
   {
      capillary_pressure = argument;
   }
   else
   {
      /*OK411
          nidx0 = PCSGetNODValueIndex("PRESSURE_CAP",0);
          nidx1 = PCSGetNODValueIndex("PRESSURE_CAP",1);
          if(mode==0){ // Gauss point values
            capillary_pressure = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
                               + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
          }
          else{ // Node values
            capillary_pressure = (1.-theta)*GetNodeVal(number,nidx0) \
                               + theta*GetNodeVal(number,nidx1);
          }
      */
   }
   //----------------------------------------------------------------------
   switch(capillary_pressure_model)
   {
      case 0:                                     // k = f(x) user-defined function
         saturation = GetCurveValueInverse((int)capillary_pressure_model_values[0],0,capillary_pressure,&gueltig);
         break;
      case 1:                                     // constant
         break;
      case 2:                                     // Lineare Kapillardruck-Saettigungs-Beziehung
         // kap12 steigt linear von 0 auf kap[2] im Bereich satu_water_saturated bis satu_water_residual
         break;
      case 3:                                     // Parabolische Kapillardruck-Saettigungs-Beziehung
         // kap12 steigt parabolisch von 0 auf kap[2] im Bereich satu_water_saturated bis satu_water_residual mit Exponent mat.kap[3]
         break;
      case 4:                                     // Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 Page 894 Equation 21
         m_mfp = mfp_vector[phase];
         density_fluid = m_mfp->Density();
         van_beta = 1/(1-saturation_exp[phase]) ;
         if(capillary_pressure < MKleinsteZahl)
            saturation = saturation_max[phase];
         else
            saturation = pow((pow(capillary_pressure*capillary_pressure_model_values[0]/density_fluid/gravity_constant,van_beta)+1),-1.*saturation_exp[phase]) \
               *(saturation_max[phase]-saturation_res[phase]) + saturation_res[phase];
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;

         /*OK
               a = kap[2];             // Test 0.005
               m = kap[3];             // Test 0.5
               n = kap[4];             // Test 2.
              saturation = pow(1.0 + pow(a*cap[1],n),-m) * (kap[1] - kap[0]) + kap[0];
         */
         break;
      case 5:                                     // Haverkamp Problem: aus SOIL SCI. SOC. AM. J., VOL. 41, 1977, Pages 285ff
         /*
               a = kap[2];             // 1.611e6 fuer Haverkamp
               b = kap[3];             // 3.96 fuer Haverkamp
               // SB: Umrechnen von cap in p_cap [cm]
               p_cap = cap[1] / GetFluidDensity(1, element, 0., 0., 0., 0.) / gravity_constant* 100.0;
               satu_return[1] = (a * (kap[1] - kap[0]))/(a + pow(p_cap,b)) + kap[0];
         */
         break;
      case 6:                                     // Brooks & Corey. 22.05.2008. WW
         if(capillary_pressure < MKleinsteZahl)
            saturation = saturation_max[phase];
         else
            saturation = (1.0-saturation_res[0]-saturation_res[1])
               *pow(capillary_pressure/capillary_pressure_model_values[0], -saturation_exp[0])
               + saturation_res[0];
         if (saturation > (saturation_max[0] - MKleinsteZahl))
            saturation = saturation_max[0] - MKleinsteZahl;
         if (saturation < (saturation_res[0] + MKleinsteZahl))
            saturation = saturation_res[0] + MKleinsteZahl;
         break;
      case 7:                                     // Van Genuchten: Wasser/Gas aus SOIL SIC. SOC. AM. J. VOL. 44, 1980 Page 894 Equation 21
         break;
      case 8:                                     // 3Phasen ueber Kurven
         break;
      case 16:                                    // JOD not needed??
         m_mfp = mfp_vector[phase];
         density_fluid = m_mfp->Density();
         van_beta = 1/(1-saturation_exp[phase]) ;
         if(capillary_pressure < MKleinsteZahl)
            saturation = saturation_max[phase];
         else
            saturation = pow((pow(capillary_pressure*saturation_alpha[phase]/density_fluid/gravity_constant,van_beta)+1),-1.*saturation_exp[phase]) \
               *(saturation_max[phase]-saturation_res[phase]) + saturation_res[phase];
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
                                                  /* Mehr als Vollsaettigung mit Wasser */
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
                                                  /* Weniger als Residualsaettigung Wasser */
            saturation = saturation_res[phase] + MKleinsteZahl;
         break;
      default:
         cout << "Error in CMediumProperties::SaturationCapillaryPressureFunction: no valid material model" << endl;
         break;
   }
   return saturation;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
08/1999 CT Erste Version (MMTM0699GetSaturationPressureDependency)
05/2001 OK Verallgemeinerung
02/2005 OK CMediumProperties function
08/2005 WW Remove interpolation
03/2007 WW Analytical solution:
02/2008 PCH Brooks-Corey dPc/dSw added
Last modified:
**************************************************************************/
double CMediumProperties::SaturationPressureDependency
(double saturation, double density_fluid, double theta)
//(int phase, long index, double r, double s, double t, double theta)
{

   static double capillary_pressure,capillary_pressure1,capillary_pressure2;
   static double saturation1,saturation2;
   static double dS,dS_dp,dpc;
   double S_e, m, n, alpha;
   int phase = 0;
   S_e=m=n=alpha=0.0;
   //----------------------------------------------------------------------
   // Vollsaettigung?
   mode = 2;
   phase = 0;                                     //WW (int)mfp_vector.size()-1;
   // Althogh p_c is known, we need to calculate it again just to fit the range of p_c>0.0
   // for this derivative calculation
   capillary_pressure = CapillaryPressureFunction(number,NULL,theta,phase, saturation);
   if(capillary_pressure < MKleinsteZahl)
      return 0.;
   switch(capillary_pressure_model)               // 01.3.2007 WW
   {
      default:                                    // k = f(x) user-defined function
         //----------------------------------------------------------------------
         // Wenn wir nah an der Vollsaettigung, ggf. Schrittweite verkleinern
         dS = 1.e-2;
         do
         {
            dS /= 10.;
            saturation1 = saturation - dS;
            capillary_pressure1 = CapillaryPressureFunction(number,NULL,theta,phase,saturation1);
            saturation2 = saturation + dS;
            capillary_pressure2 = CapillaryPressureFunction(number,NULL,theta,phase,saturation2);
                                                  //OK4105
            dpc = capillary_pressure1 - capillary_pressure2;
         }
         while((dS > MKleinsteZahl) && (capillary_pressure2 < MKleinsteZahl / 100.));
         if ( ((capillary_pressure1 > MKleinsteZahl)||(capillary_pressure2 > MKleinsteZahl)) \
            && (dpc > MKleinsteZahl) )
            dS_dp = 2. * dS / (capillary_pressure1 - capillary_pressure2);
         else
         {
            dS_dp = 0.;
            //cout << "Warning in CMediumProperties::SaturationPressureDependency: dS_dp = 0" << endl; //OK4105
         }
         break;
      case 44:                                    // Van Genuchten: 01.3.2007 WW
         m = saturation_exp[phase];
         n = 1./(1.0-m);
         alpha = capillary_pressure_model_values[0]/(density_fluid*gravity_constant);
         // saturation is p_c
         // dS_dp = dS_d/p_c * dp_c/dp = -dS_d/p_c
         dS_dp = m*n*(saturation_max[phase]-saturation_res[phase])*
            pow(1.0+pow(alpha*capillary_pressure,n),-m-1.0)
            *pow(alpha, n)*pow(capillary_pressure, n-1.0);
         break;
   }
   //----------------------------------------------------------------------
   mode = 0;
   return dS_dp;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2009 PCH Brooks-Corey dPc/dSw added
Last modified:
**************************************************************************/
double CMediumProperties::PressureSaturationDependency
(double saturation, double density_fluid, double theta)
{

   density_fluid = density_fluid;                 //OK411

   static double capillary_pressure,capillary_pressure1,capillary_pressure2;
   static double saturation1,saturation2;
   static double dS,dpc;
   // static double dS_dp;
   // 01.3.2007 WW
   double S_e, m, n, alpha, dPcdSe;
   int phase = 0;
   S_e=m=n=alpha=0.0;

   phase = 0;                                     //WW (int)mfp_vector.size()-1;
   // Althogh p_c is known, we need to calculate it again just to fit the range of p_c>0.0
   // for this derivative calculation
   capillary_pressure = CapillaryPressureFunction(number,NULL,theta,phase, saturation);
   //	if(capillary_pressure < MKleinsteZahl)
   //		return 0.;
   switch(capillary_pressure_model)               // 01.3.2007 WW
   {
      case 0:                                     // k = f(x) user-defined function
         dS = 1.e-2;
         do
         {
            dS /= 10.;
            saturation1 = saturation - dS;
            capillary_pressure1 = CapillaryPressureFunction(number,NULL,theta,phase,saturation1);
            saturation2 = saturation + dS;
            capillary_pressure2 = CapillaryPressureFunction(number,NULL,theta,phase,saturation2);
                                                  //OK4105
            dpc = capillary_pressure1 - capillary_pressure2;
         }
         while((dS > MKleinsteZahl) && (capillary_pressure2 < MKleinsteZahl / 100.));

         dPcdSe = (capillary_pressure1 - capillary_pressure2)/(2. * dS);
         break;
      default:                                    // k = f(x) user-defined function
         dS = 1.e-2;
         do
         {
            dS /= 10.;
            saturation1 = saturation - dS;
            capillary_pressure1 = CapillaryPressureFunction(number,NULL,theta,phase,saturation1);
            saturation2 = saturation + dS;
            capillary_pressure2 = CapillaryPressureFunction(number,NULL,theta,phase,saturation2);
                                                  //OK4105
            dpc = capillary_pressure1 - capillary_pressure2;
         }
         while((dS > MKleinsteZahl) && (capillary_pressure2 < MKleinsteZahl / 100.));

         dPcdSe = (capillary_pressure1 - capillary_pressure2)/(2. * dS);
         break;

      case 6:                                     // Brooks-Corey: 02.2008 PCH
         if (saturation > (saturation_max[phase] - MKleinsteZahl))
            saturation = saturation_max[phase] - MKleinsteZahl;
         if (saturation < (saturation_res[phase] + MKleinsteZahl))
            saturation = saturation_res[phase] + MKleinsteZahl;

         S_e = (saturation-saturation_res[phase])/(saturation_max[phase]-saturation_res[phase]);

         dPcdSe=-capillary_pressure_model_values[0]/(saturation_exp[0]*S_e) *
            pow(S_e, -1.0/saturation_exp[0]);
         break;
   }

   return dPcdSe;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
08/1999 CT Erste Version (MMTM0699GetSaturationPressureDependency)
05/2001 OK Verallgemeinerung
02/2005 OK CMediumProperties function
Last modified:
**************************************************************************/
double CMediumProperties::SaturationPressureDependency(long number,double*gp,double theta)
//(int phase, long index, double r, double s, double t, double theta)
{
   //OK411
   gp = gp;

   static double capillary_pressure,capillary_pressure1,capillary_pressure2;
   static double saturation,saturation1,saturation2;
   static double dS,dS_dp,dpc;
   //int nidx0,nidx1;
   int phase = (int)mfp_vector.size()-1;
   //----------------------------------------------------------------------
   if(mode==2)
   {
      saturation = argument;
   }
   else
   {
      /*OK411
          nidx0 = PCSGetNODValueIndex("SATURATION1",0);
          nidx1 = PCSGetNODValueIndex("SATURATION1",1);
          if(mode==0){ // Gauss point values
            saturation = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
                     + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
          }
          else // Node values
          {
            saturation = SaturationCapillaryPressureFunction(number,NULL,theta,phase);
            saturation = (1.-theta)*GetNodeVal(number,nidx0) \
      + theta*GetNodeVal(number,nidx1);
      }
      */
   }
   //----------------------------------------------------------------------
   // Vollsaettigung?
   mode = 2;
   capillary_pressure = CapillaryPressureFunction(number,NULL,theta,phase,saturation);
   if(capillary_pressure < MKleinsteZahl)
      return 0.;
   //----------------------------------------------------------------------
   // Wenn wir nah an der Vollsaettigung, ggf. Schrittweite verkleinern
   dS = 1.e-2;
   do
   {
      dS /= 10.;
      saturation1 = saturation - dS;
      capillary_pressure1 = CapillaryPressureFunction(number,NULL,theta,phase,saturation1);
      saturation2 = saturation + dS;
      capillary_pressure2 = CapillaryPressureFunction(number,NULL,theta,phase,saturation2);
                                                  //OK4105
      dpc = capillary_pressure1 - capillary_pressure2;
   }
   while((dS > MKleinsteZahl) && (capillary_pressure2 < MKleinsteZahl / 100.));
   if ( ((capillary_pressure1 > MKleinsteZahl)||(capillary_pressure2 > MKleinsteZahl)) \
      && (dpc > MKleinsteZahl) )
      dS_dp = 2. * dS / (capillary_pressure1 - capillary_pressure2);
   else
   {
      dS_dp = 0.;
      //cout << "Warning in CMediumProperties::SaturationPressureDependency: dS_dp = 0" << endl; //OK4105
   }
   //----------------------------------------------------------------------
   mode = 0;
   return dS_dp;
}


/*************************************************************************************************************
14.$MASSDISPERSION
************************************************************************************************************/

/*************************************************************************************************************
Density
************************************************************************************************************/
/**************************************************************************
FEMLib-Method:
Task:
Programing:
11/2004 OK Implementation
**************************************************************************/
double CMediumProperties::Density(long element,double*gp,double theta)
{
   //OK411
   gp = gp;

   int no_phases = (int)mfp_vector.size();
   double density = 0.0;
   int i;
   CFluidProperties* m_mfp = NULL;
   //OK411 CSolidProperties* m_msp = NULL;
   char saturation_name[15];
   if(no_phases==1)
   {
      m_mfp = mfp_vector[0];
      density = Porosity(element,theta) * m_mfp->Density();
   }
   else
   {
      for(i=0;i<no_phases;i++)
      {
         m_mfp = mfp_vector[i];
         sprintf(saturation_name,"SATURATION%i",i+1);
         //OK411 density += Porosity(element,theta) * m_mfp->Density() * PCSGetELEValue(element,gp,theta,saturation_name);
      }
   }
   /*OK411
     long group = ElGetElementGroupNumber(element);
     m_msp = msp_vector[group];
     density += (1.-Porosity(element,theta))*fabs(m_msp->Density());
   */
   return density;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
last modified:
**************************************************************************/
void MMPDelete()
{
   long i;
   int no_mmp =(int)mmp_vector.size();
   for(i=0;i<no_mmp;i++)
   {
      delete mmp_vector[i];
   }
   mmp_vector.clear();
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
06/2005 OK Implementation
07/2005 WW Change due to geometry element object applied
last modified:
**************************************************************************/
void MMP2PCSRelation(CRFProcess*m_pcs)
{
   CMediumProperties*m_mmp = NULL;
   Mesh_Group::CElem* m_ele = NULL;
   if(m_pcs->m_msh)
   {
      for(long i=0;i<(long)m_pcs->m_msh->ele_vector.size();i++)
      {
         m_ele = m_pcs->m_msh->ele_vector[i];
         if(m_ele->GetPatchIndex()<(int)mmp_vector.size())
         {
            m_mmp = mmp_vector[m_ele->GetPatchIndex()];
            m_mmp->m_pcs = m_pcs;
         }
      }
   }
   else
   {
      for(int j=0;j<(int)mmp_vector.size();j++)
      {
         m_mmp = mmp_vector[j];
         m_mmp->m_pcs = m_pcs;
      }
   }
}


/**************************************************************************
PCSLib-Function
Liest zu jedem Knoten einen Wert der Permeabilität ein.
Identifikation über Koordinaten
Programing:
10/2003     SB  First Version
01/2004     SB  2. Version
01/2005 OK Check MAT groups //OK41
06/2005 MB msh, loop over mmp groups
09/2005 MB EleClass
//SB/MB ? member function of CFEMesh / CMediumProperties
11/2005 OK GEOMETRY_AREA
04/2006 CMCD Constant area
**************************************************************************/
//MMPGetHeterogeneousFields
void GetHeterogeneousFields()
{
   //OK411 int ok=0;
   //OK411 char* name_file=NULL;
   CMediumProperties *m_mmp = NULL;
   //----------------------------------------------------------------------
   // File handling
   string file_path;
   string file_path_base_ext;
   CGSProject* m_gsp = NULL;
   m_gsp = GSPGetMember("mmp");
   if(m_gsp)
      file_path = m_gsp->path;

   //----------------------------------------------------------------------
   // Tests
   if(mmp_vector.size()==0)
      return;
   //----------------------------------------------------------------------
   //Schleife über alle Gruppen
   for(int i=0;i<(int)mmp_vector.size();i++)
   {
      m_mmp = mmp_vector[i];
      //....................................................................
      // For Permeability
      if(m_mmp->permeability_file.size()>0)
      {
         //OK name_file = (char *) m_mmp->permeability_file.data();
         //OK if(name_file != NULL)
         //OK ok = FctReadHeterogeneousFields(name_file,m_mmp);

         //WW file_path_base_ext = file_path + m_mmp->permeability_file;
                                                  //WW
         m_mmp->SetDistributedELEProperties(m_mmp->permeability_file);
         m_mmp->WriteTecplotDistributedProperties();
      }
      //....................................................................
      // For Porosity
      if(m_mmp->porosity_file.size()>0)
      {
         //CB name_file = (char *) m_mmp->porosity_file.data();
         //CB if(name_file != NULL)
         //CB  ok = FctReadHeterogeneousFields(name_file,m_mmp);
         //file_path_base_ext = file_path + m_mmp->porosity_file;
         //m_mmp->SetDistributedELEProperties(file_path_base_ext); // CB Removed bugs in this function
                                                  // CB Removed bugs in this function
         m_mmp->SetDistributedELEProperties(m_mmp->porosity_file);
         m_mmp->WriteTecplotDistributedProperties();
      }
      //....................................................................
      // GEOMETRY_AREA
      if(m_mmp->geo_area_file.size()>0)
      {
         file_path_base_ext = file_path + m_mmp->geo_area_file;
         m_mmp->SetDistributedELEProperties(file_path_base_ext);
         m_mmp->WriteTecplotDistributedProperties();
      }
      //NW    else m_mmp->SetConstantELEarea(m_mmp->geo_area,i);
      //....................................................................
   }
   //----------------------------------------------------------------------
}


/**************************************************************************
PCSLib-Method:
Programing:
04/2006 CMCD Implementation
**************************************************************************/
void CMediumProperties::SetConstantELEarea(double area, int group)
{
   long i,j, ele_group;
   long no_ele;
   int no_processes = (int) pcs_vector.size();
   if (area != 1.0)
   {
      for (i=0; i<no_processes; i++)
      {
         m_msh = FEMGet(convertProcessTypeToString (pcs_vector[i]->getProcessType()));
         if(!m_msh) return;                       //WW
         no_ele = (long) m_msh->ele_vector.size();
         for(j=0;j<no_ele;j++)
         {
            ele_group = m_msh->ele_vector[j]->GetPatchIndex();
            if (ele_group == group) m_msh->ele_vector[j]->SetFluxArea(area);
         }
      }
   }
}


/**************************************************************************
PCSLib-Method:
Programing:
11/2005 OK Implementation
**************************************************************************/
void CMediumProperties::SetDistributedELEProperties(string file_name)
{
   string line_string, line1;
   string mmp_property_name;
   string mmp_property_dis_type;
   string mmp_property_mesh;
   CElem* m_ele_geo = NULL;
   bool element_area = false;
   long i, j, ihet;
   double mmp_property_value;
   int mat_vector_size=0;                         // Init WW
   double ddummy, conversion_factor=1.0;;         //init WW
   vector <double> xvals, yvals, zvals, mmpvals;
   vector<double>temp_store;
   int c_vals;
   double x, y, z, mmpv;
   std::stringstream in;
   //CB
   vector<double> garage;
   int mat_vec_size=0;
   int por_index = 0;
   int vol_bio_index = 0;
   string outfile;
   int k;

   cout << " SetDistributedELEProperties: " ;
   //----------------------------------------------------------------------
   // File handling
   ifstream mmp_property_file(file_name.data(),ios::in);
   if(!mmp_property_file.good())
   {
      cout << "Warning in CMediumProperties::SetDistributedELEProperties: no MMP property data" << endl;
      return;
   }
   mmp_property_file.clear();
   mmp_property_file.seekg(0,ios::beg);
   //----------------------------------------------------------------------
   line_string = GetLineFromFile1(&mmp_property_file);
   if(!(line_string.find("#MEDIUM_PROPERTIES_DISTRIBUTED")!=string::npos))
   {
      cout << "Keyword #MEDIUM_PROPERTIES_DISTRIBUTED not found" << endl;
      return;
   }
   //----------------------------------------------------------------------
   while(!mmp_property_file.eof())
   {
      line_string = GetLineFromFile1(&mmp_property_file);
      if(line_string.find("STOP")!=string::npos)
         return;
      if(line_string.empty())
      {
         cout << "Error in CMediumProperties::SetDistributedELEProperties - no enough data sets" << endl;
         return;
      }
      //....................................................................
      if(line_string.find("$MSH_TYPE")!=string::npos)
      {
         line_string = GetLineFromFile1(&mmp_property_file);
         mmp_property_mesh = line_string;
         m_msh = FEMGet(line_string);
         if(!m_msh)
         {
            cout << "CMediumProperties::SetDistributedELEProperties: no MSH data" << endl;
            return;
         }
         continue;
      }
      //....................................................................
      if(line_string.find("$MMP_TYPE")!=string::npos)
      {
         element_area = false;
         mmp_property_file >> mmp_property_name;
         cout << mmp_property_name << endl;
         m_msh->mat_names_vector.push_back(mmp_property_name);
         if (mmp_property_name == "GEOMETRY_AREA") element_area = true;
         continue;
      }
      //....................................................................
      if(line_string.find("$DIS_TYPE")!=string::npos)
      {
         mmp_property_file >> mmp_property_dis_type;
         continue;
      }
      //....................................................................
      if(line_string.find("$CONVERSION_FACTOR")!=string::npos)
      {
         mmp_property_file >> conversion_factor;
         continue;
      }
      //....................................................................
      if(line_string.find("$DATA")!=string::npos)
      {
         switch(mmp_property_dis_type[0])
         {
            case 'N':                             // Next neighbour
            case 'G':                             // Geometric mean
               // Read in all values given, store in vectors for x, y, z and value
               i=0;
               while(i==0)
               {
                  line1 = GetLineFromFile1(&mmp_property_file);
                  if(line1.find("STOP")!=string::npos) break;
                  in.str((string)line1);
                  in >> x >> y >> z >> mmpv;
                  in.clear();
                  mmpv *= conversion_factor;      // convert values
                  xvals.push_back(x);
                  yvals.push_back(y);
                  zvals.push_back(z);
                  mmpvals.push_back(mmpv);
               }
               // sort values to mesh
               for(i=0;i<(long)m_msh->ele_vector.size();i++)
               {
                  m_ele_geo = m_msh->ele_vector[i];
                  mat_vector_size = m_ele_geo->mat_vector.Size();
                  // CB Store old values as they are set to zero after resizing
                  for(j=0;j<mat_vector_size;j++)
                     garage.push_back(m_ele_geo->mat_vector(j)) ;
                  m_ele_geo->mat_vector.resize(mat_vector_size+1);
                  // CB Refill old values as they were set to zero after resizing
                  for(j=0;j<mat_vector_size;j++)
                     m_ele_geo->mat_vector(j) = garage[j];
                  garage.clear();
                  if(mmp_property_dis_type[0] == 'N')
                  {
                     // Search for all elements of the mesh, which is the nearest given value in the input file
                     // Return value ihet is the index of the het. val in the mmpval-vector
                     ihet = GetNearestHetVal2(i, m_msh, xvals, yvals, zvals, mmpvals);
                     m_ele_geo->mat_vector(mat_vector_size) = mmpvals[ihet];
                  }
                  if(mmp_property_dis_type[0] == 'G')
                  {
                     mmpv = GetAverageHetVal2(i, m_msh, xvals, yvals, zvals, mmpvals);
                     m_ele_geo->mat_vector(mat_vector_size) = mmpv;
                  }
               }
               break;
            case 'E':                             // Element data
               for(i=0;i<(long)m_msh->ele_vector.size();i++)
               {
                  m_ele_geo = m_msh->ele_vector[i];
                  mmp_property_file >> ddummy >> mmp_property_value;
                  mat_vector_size = m_ele_geo->mat_vector.Size();
                  if (mat_vector_size > 0)
                  {
                     for (c_vals = 0; c_vals < mat_vector_size; c_vals++)
                        temp_store.push_back(m_ele_geo->mat_vector(c_vals));
                     m_ele_geo->mat_vector.resize(mat_vector_size+1);
                     for (c_vals = 0; c_vals < mat_vector_size; c_vals++)
                        m_ele_geo->mat_vector(c_vals)=temp_store[c_vals];
                     m_ele_geo->mat_vector(mat_vector_size) = mmp_property_value;
                     temp_store.clear();
                  }
                  else
                  {
                     m_ele_geo->mat_vector.resize(mat_vector_size+1);
                     m_ele_geo->mat_vector(mat_vector_size) = mmp_property_value;
                  }
                  if (element_area) m_msh->ele_vector[i]->SetFluxArea(mmp_property_value);
                  if(line_string.empty())
                  {
                     cout << "Error in CMediumProperties::SetDistributedELEProperties - not enough data sets" << endl;
                     return;
                  }
               }
               break;
            default:
               cout << " Unknown interpolation option for the values!" << endl;
               break;
         }
         continue;
      }
      //....................................................................
   }
   // CB now set VOL_MAT & VOL_BIO as heterogeneous values, if defined as model 2 and het Porosity
   if( (mmp_property_name == "POROSITY") && (this->vol_bio_model == 2) )
   {
      m_msh->mat_names_vector.push_back("VOL_BIO");
      for(i=0;i<(long)m_msh->ele_vector.size();i++)
      {
         m_ele_geo = m_msh->ele_vector[i];        // Get the element
         mat_vec_size = m_ele_geo->mat_vector.Size();
         // CB Store old values as they are set to zero after resizing
         for(j=0;j<mat_vec_size;j++)
            garage.push_back(m_ele_geo->mat_vector(j)) ;
         m_ele_geo->mat_vector.resize(mat_vec_size+1);
         // CB Refill old values as they were set to zero after resizing
         for(j=0;j<mat_vec_size;j++)
            m_ele_geo->mat_vector(j) = garage[j];
         garage.clear();
         // Set the VOL_BIO value from mmp file input
         m_ele_geo->mat_vector(mat_vec_size) = this->vol_bio;
      }
   }
   if( (mmp_property_name == "POROSITY") && (this->vol_mat_model == 2) )
   {
      m_msh->mat_names_vector.push_back("VOL_MAT");
      // Get the porosity index
      for(por_index=0;por_index<(int)m_msh->mat_names_vector.size();por_index++)
         if(m_msh->mat_names_vector[por_index].compare("POROSITY")==0)
            break;
      // Get the vol_bio index
      for(vol_bio_index=0;vol_bio_index<(int)m_msh->mat_names_vector.size();vol_bio_index++)
         if(m_msh->mat_names_vector[vol_bio_index].compare("VOL_BIO")==0)
            break;
      for(i=0;i<(long)m_msh->ele_vector.size();i++)
      {
         m_ele_geo = m_msh->ele_vector[i];        // Get the element
         mat_vec_size = m_ele_geo->mat_vector.Size();
         // CB Store old values as they are set to zero after resizing
         for(j=0;j<mat_vec_size;j++)
            garage.push_back(m_ele_geo->mat_vector(j)) ;
         m_ele_geo->mat_vector.resize(mat_vec_size+1);
         // CB Refill old values as they were set to zero after resizing
         for(j=0;j<mat_vec_size;j++)
            m_ele_geo->mat_vector(j) = garage[j];
         garage.clear();
         // Set the VOL_MAT value from (1-POROSITY-VOL_BIO)
         m_ele_geo->mat_vector(mat_vec_size) = 1 - m_ele_geo->mat_vector(por_index) - m_ele_geo->mat_vector(vol_bio_index) ;
      }
   }
   //----------------------------------------------------------------------
   //Write sorted output file
   //----------------------------------------------------------------------
   // File handling

                                                  // CB
   for(k=0;k<(int)m_msh->mat_names_vector.size();k++)
   {
      //file_name +="_sorted";
      outfile = m_msh->mat_names_vector[k] + "_sorted";
      ofstream mmp_property_file_out(outfile.data());
      if(!mmp_property_file_out.good())
      {
         cout << "Warning in CMediumProperties::WriteDistributedELEProperties: no MMP property data file to write to" << endl;
         return;
      }
      mmp_property_file_out << "#MEDIUM_PROPERTIES_DISTRIBUTED" << endl;
      mmp_property_file_out << "$MSH_TYPE" << endl << "  " << mmp_property_mesh << endl;
      //mmp_property_file_out << "$MSH_TYPE" << endl << "  " << mmp_property_mesh << endl;
      //mmp_property_file_out << "$MMP_TYPE" << endl << "  " << "PERMEABILITY" << endl;
      mmp_property_file_out << "$MMP_TYPE" << endl << "  " << m_msh->mat_names_vector[k] << endl;
      mmp_property_file_out << "$DIS_TYPE" << endl << "  " << "ELEMENT" << endl;
      mmp_property_file_out << "$DATA" << endl ;
      for(i=0;i<(long)m_msh->ele_vector.size();i++)
      {
         m_ele_geo = m_msh->ele_vector[i];
         mmp_property_file_out << i << "  " << m_ele_geo->mat_vector(k) << endl;
      }
      mmp_property_file_out << "#STOP" << endl;
      mmp_property_file_out.close();
      //----------------------------------------------------------------------
   }

}


/**************************************************************************
PCSLib-Method:
Programing:
11/2005 OK Implementation
**************************************************************************/
void CMediumProperties::WriteTecplotDistributedProperties()
{
   int j, k;
   long i;
   string element_type;
   string m_string = "MAT";
   double m_mat_prop_nod;
   //----------------------------------------------------------------------
   // Path
   string path;
   CGSProject* m_gsp = NULL;
   m_gsp = GSPGetMember("msh");
   if (m_gsp)
   {
      path = m_gsp->path;
   }
   //--------------------------------------------------------------------
   // MSH
   CNode* m_nod = NULL;
   CElem* m_ele = NULL;
   if (!m_msh)
      return;
   //--------------------------------------------------------------------
   // File handling
   string mat_file_name = path + name + "_" + m_msh->pcs_name + "_PROPERTIES"
      + TEC_FILE_EXTENSION;
   fstream mat_file(mat_file_name.data(), ios::trunc | ios::out);
   mat_file.setf(ios::scientific, ios::floatfield);
   mat_file.precision(12);
   if (!mat_file.good())
      return;
   mat_file.seekg(0L, ios::beg);
   //--------------------------------------------------------------------
   if ((long) m_msh->ele_vector.size() > 0)
   {
      m_ele = m_msh->ele_vector[0];
      switch (m_ele->GetElementType())
      {
         case MshElemType::LINE:
            element_type = "QUADRILATERAL";
            break;
         case MshElemType::QUAD:
            element_type = "QUADRILATERAL";
            break;
         case MshElemType::HEXAHEDRON:
            element_type = "BRICK";
            break;
         case MshElemType::TRIANGLE:
            element_type = "TRIANGLE";
            break;
         case MshElemType::TETRAHEDRON:
            element_type = "TETRAHEDRON";
            break;
         case MshElemType::PRISM:
            element_type = "BRICK";
            break;
         default:
            std::cerr
               << "CMediumProperties::WriteTecplotDistributedProperties MshElemType not handled"
               << std::endl;
      }
   }
   //--------------------------------------------------------------------
   // Header
   mat_file << "VARIABLES = X,Y,Z";
   for (j = 0; j < (int) m_msh->mat_names_vector.size(); j++)
   {
      mat_file << "," << m_msh->mat_names_vector[j];
   }
   mat_file << endl;
   mat_file << "ZONE T = " << name << ", " << "N = "
      << (long) m_msh->nod_vector.size() << ", " << "E = "
      << (long) m_msh->ele_vector.size() << ", " << "F = FEPOINT" << ", "
      << "ET = " << element_type << endl;
   //--------------------------------------------------------------------
   // Nodes
   for (i = 0; i < (long) m_msh->nod_vector.size(); i++)
   {
      m_nod = m_msh->nod_vector[i];
      mat_file << m_nod->X() << " " << m_nod->Y() << " " << m_nod->Z();
      for (j = 0; j < (int) m_msh->mat_names_vector.size(); j++)
      {
         m_mat_prop_nod = 0.0;
         for (k = 0; k < (int) m_nod->connected_elements.size(); k++)
         {
            m_ele = m_msh->ele_vector[m_nod->connected_elements[k]];
            m_mat_prop_nod += m_ele->mat_vector(j);
         }
         m_mat_prop_nod /= (int) m_nod->connected_elements.size();
         mat_file << " " << m_mat_prop_nod;
      }
      mat_file << endl;
   }
   //--------------------------------------------------------------------
   // Elements
   for (i = 0; i < (long) m_msh->ele_vector.size(); i++)
   {
      m_ele = m_msh->ele_vector[i];
      //OK if(m_ele->GetPatchIndex()==number) {
      switch (m_ele->GetElementType())
      {
         case MshElemType::LINE:
            mat_file << m_ele->nodes_index[0] + 1 << " "
               << m_ele->nodes_index[1] + 1 << " "
               << m_ele->nodes_index[1] + 1 << " "
               << m_ele->nodes_index[0] + 1 << endl;
            element_type = "QUADRILATERAL";
            break;
         case MshElemType::QUAD:
            mat_file << m_ele->nodes_index[0] + 1 << " "
               << m_ele->nodes_index[1] + 1 << " "
               << m_ele->nodes_index[2] + 1 << " "
               << m_ele->nodes_index[3] + 1 << endl;
            element_type = "QUADRILATERAL";
            break;
         case MshElemType::HEXAHEDRON:
            mat_file << m_ele->nodes_index[0] + 1 << " "
               << m_ele->nodes_index[1] + 1 << " "
               << m_ele->nodes_index[2] + 1 << " "
               << m_ele->nodes_index[3] + 1 << " "
               << m_ele->nodes_index[4] + 1 << " "
               << m_ele->nodes_index[5] + 1 << " "
               << m_ele->nodes_index[6] + 1 << " "
               << m_ele->nodes_index[7] + 1 << endl;
            element_type = "BRICK";
            break;
         case MshElemType::TRIANGLE:
            mat_file << m_ele->nodes_index[0] + 1 << " "
               << m_ele->nodes_index[1] + 1 << " "
               << m_ele->nodes_index[2] + 1 << endl;
            element_type = "TRIANGLE";
            break;
         case MshElemType::TETRAHEDRON:
            mat_file << m_ele->nodes_index[0] + 1 << " "
               << m_ele->nodes_index[1] + 1 << " "
               << m_ele->nodes_index[2] + 1 << " "
               << m_ele->nodes_index[3] + 1 << endl;
            element_type = "TETRAHEDRON";
            break;
         case MshElemType::PRISM:
            mat_file << m_ele->nodes_index[0] + 1 << " "
               << m_ele->nodes_index[0] + 1 << " "
               << m_ele->nodes_index[1] + 1 << " "
               << m_ele->nodes_index[2] + 1 << " "
               << m_ele->nodes_index[3] + 1 << " "
               << m_ele->nodes_index[3] + 1 << " "
               << m_ele->nodes_index[4] + 1 << " "
               << m_ele->nodes_index[5] + 1 << endl;
            element_type = "BRICK";
            break;
         default:
            std::cerr
               << "CMediumProperties::WriteTecplotDistributedProperties MshElemType not handled"
               << std::endl;
      }
   }
}


/**************************************************************************
MSHLib-Method: GetNearestHetVal2
Task:
Programing:
0?/2004 SB Implementation
09/2005 MB EleClass
01/2006 SB ReImplementation with new concept by Olaf, moved here
**************************************************************************/
long GetNearestHetVal2(long EleIndex, CFEMesh *m_msh, vector <double> xvals,  vector <double> yvals,  vector <double> zvals,  vector <double> mmpvals)
{

   (void)mmpvals;
   long i, nextele, no_values;
   double ex, ey, ez, dist, dist1, dist2;
   double x, y, z;
   double* center = NULL;
   Mesh_Group::CElem* m_ele = NULL;
   no_values = (long) xvals.size();

   x=0.0; y=0.0; z=0.0;
   dist = 10000000.0;                             //Startwert
   dist2 = 0.01;                                  // Abstand zwischen eingelesenen Knoten und Geometrieknoten-RF;
   // Achtung, doppelbelegung möglich bei kleinen Gitterabständen
   nextele = -1;

   //Get element data
   m_ele = m_msh->ele_vector[EleIndex];
   center = m_ele->GetGravityCenter();
   x = center[0];
   y = center[1];
   z = center[2];

   //Calculate distances
   for(i=0;i<no_values;i++)
   {
      ex=xvals[i];
      ey=yvals[i];
      ez=zvals[i];
      dist1 = (ex-x)*(ex-x)+(ey-y)*(ey-y)+(ez-z)*(ez-z);
      if(dist1<dist)
      {
         dist = dist1;
         nextele = i;
      }
   }

   return nextele;
}


/**************************************************************************
MSHLib-Method: GetAverageHetVal2
Task:
Programing:
06/2005 MB Implementation
01/2006 SB Adapted to new structure
**************************************************************************/
double GetAverageHetVal2(long EleIndex, CFEMesh *m_msh, vector <double> xvals,  vector <double> yvals,  vector <double> zvals,  vector <double> mmpvals)
{

   long i, j, ihet;
   double average;
   double xp[3],yp[3];
   double value;
   double NumberOfValues;
   double InvNumberOfValues;
   CGLPoint *m_point = NULL;
   Mesh_Group::CElem* m_ele = NULL;
   long   no_values = (long) xvals.size();

   j = 0;                                         //only for 1 value

   //-----------------------------------------------------------------------
   //Get element data
   m_ele = m_msh->ele_vector[EleIndex];
   for(j=0;j<3; j++)
   {
      xp[j] = m_ele->GetNode(j)->X();
      yp[j] = m_ele->GetNode(j)->Y();
      //zp[j] = 0.0;
   }

   //-----------------------------------------------------------------------
   //Find data points in the element
   NumberOfValues = 0;
   InvNumberOfValues = 0;
   m_point = new CGLPoint;

   average = -1;
   value = 0;

   for(i=0;i<no_values;i++)
   {
      if(mmpvals[i]!= -999999.0)                  //Data point not within an element yet
      {
         m_point->x = xvals[i];
         m_point->y = yvals[i];
         m_point->z = 0.0;

         //....................................................................
         //Calculate the product of values in element
                                                  //CC 10/05
         if(m_point->IsInTriangleXYProjection(xp,yp))
         {
            value = value + zvals[i];
            NumberOfValues ++;
            mmpvals[i] = -999999.0;               //used as marker
         }
      }
   }                                              //end for
   //........................................................................
   if(NumberOfValues == 0)                        //if no data points in element --> get neares value
   {
      ihet = GetNearestHetVal2(EleIndex, m_msh, xvals, yvals, zvals, mmpvals);
      if(ihet<0)
         DisplayMsgLn(" Error getting nearest het_value location");
      else
      {
         average = mmpvals[ihet];
      }
   }
   //........................................................................
   else                                           //if data points in element --> Calculate arithmetic mean
   {
      average = value / NumberOfValues;
   }
   delete m_point;
   return average;
}


/**************************************************************************
FEMLib-Method:
Task: set MMP group member
Programing:
01/2006 YD Implementation
**************************************************************************/
void CMediumPropertiesGroup::Set(CRFProcess* m_pcs)
{
   long j,k;
   CFEMesh* m_msh = m_pcs->m_msh;
   CMediumProperties* m_mmp = NULL;
   CElem* elem = NULL;
   //----------------------------------------------------------------------
   // Tests //
   if(!m_msh)
   {
      cout << "Warning in CSourceTermGroup::Set - no MSH data" << endl;
      //return;
   }
   //----------------------------------------------------------------------
   long no_mmp =(long)mmp_vector.size();
   for(j=0;j<no_mmp;j++)
   {
      m_mmp = mmp_vector[j];
      //====================================================================
      if(m_mmp->pcs_type_name.compare(pcs_type_name)==0)
      {
         m_mmp = mmp_vector[j];
         for(k=0;k<(long)m_msh->ele_vector.size();k++)
         {
            elem = m_msh->ele_vector[k];
            elem->SetPatchIndex(j);
         }
      }
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2007 WW Implementation
**************************************************************************/
void CMediumProperties::CalStressPermeabilityFactor(double *kfac, const double T)
{
   switch(permeability_stress_mode)
   {
      case 2:
         CalStressPermeabilityFactor2(kfac, T);
         break;
      case 3:
         CalStressPermeabilityFactor3(kfac);
         break;
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2007 WW Implementation
**************************************************************************/
void CMediumProperties::CalStressPermeabilityFactor2(double *kfac, const double T)
{
   int i, ia, ib, ele_index;
   double xyz[3], sig[3], b[3], b0;
   ElementValue_DM *e_valDM = NULL;
   ele_index = Fem_Ele_Std->Index;
   Fem_Ele_Std->RealCoordinates(xyz);
   b0 = pow(6.0*permeability_tensor[0]/c_coefficient[5], 1.0/3.0);
   if((long)ele_value_dm.size()>0)                // Deformation process coupled
   {
      e_valDM = ele_value_dm[ele_index];
      for(i=0; i<3; i++)
      {
         sig[i] = (*e_valDM->Stress)(i, Fem_Ele_Std->GetGPindex());
         b[i] = c_coefficient[2+i]  + c_coefficient[0]*exp(c_coefficient[1]*sig[i]/(T*GAS_CONSTANT));
         //
      }
   }
   else
   {
      for(i=0; i<3; i++)
      {
         sig[i] = c_coefficient[6+i*4]+c_coefficient[7+i*4]*xyz[0]
            +c_coefficient[8+i*4]*xyz[1]+c_coefficient[9+i*4]*xyz[2];
         b[i] = c_coefficient[2+i] + c_coefficient[0]*exp(c_coefficient[1]*sig[i]/(T*GAS_CONSTANT));
      }
   }
   for(i=0; i<3; i++)
   {
      ia = (i+1)%3;
      ib = (i+2)%3;
      kfac[i] = 0.5*(pow(b[ia],3.0)+pow(b[ib],3.0))/pow(b0,3.0);
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2007 WW Implementation
**************************************************************************/
void CMediumProperties::CalStressPermeabilityFactor3(double *kfac)
{
   int i, ele_index;
   double am, pG, fy;
   double xyz[3], sig[3], ah[3];
   double JRC = c_coefficient[0];
   double a01 = c_coefficient[1];
   double Kn = c_coefficient[2];
   //
   ElementValue_DM *e_valDM = NULL;
   ele_index = Fem_Ele_Std->Index;
   pG = Fem_Ele_Std->PG;
   if(pG<0.0)  pG = 0.0;
   Fem_Ele_Std->RealCoordinates(xyz);
   if((long)ele_value_dm.size()>0)                // Deformation process coupled
   {
      e_valDM = ele_value_dm[ele_index];
      for(i=0; i<3; i++)
         sig[i] = 1.e-6*((*e_valDM->Stress)(i, Fem_Ele_Std->GetGPindex())-max(pG,0.0));
   }
   else
   {
      for(i=0; i<3; i++)
      {
         sig[i] = 1.0e-6*(c_coefficient[9+i*4]+c_coefficient[10+i*4]*xyz[0]
            +c_coefficient[11+i*4]*xyz[1]+c_coefficient[12+i*4]*xyz[2]-max(pG,0.0));
      }
   }
   // am at 100
   double am0_h = a01-(c_coefficient[7]/(-Kn+ c_coefficient[7]/c_coefficient[3])+c_coefficient[4]+c_coefficient[5]+c_coefficient[6]);
   double am0_H = a01-(c_coefficient[8]/(-Kn+ c_coefficient[8]/c_coefficient[3])+c_coefficient[4]+c_coefficient[5]+c_coefficient[6]);
   double ah0_h = am0_h*am0_h;
   double ah0_H = am0_H*am0_H;
   if(ah0_h>am0_h)  ah0_h=am0_h;
   if(ah0_H>am0_H)  ah0_H=am0_h;
   for(i=0; i<3; i++)
   {
      am =  a01-(sig[i]/(-Kn+ sig[i]/c_coefficient[3])+c_coefficient[4]+c_coefficient[5]+c_coefficient[6]);
      ah[i] = am*am;
      if(ah[i]>am) ah[i] = am;
      //
      c_coefficient[21+i] = ah[i]/pow(JRC,2.5);
   }
   kfac[0] = ah[0]*ah[0]/(ah0_h*ah0_h);
   if(geo_dimension==2)
   {
      fy = ah[2]*ah[2]/(ah0_H*ah0_H);
      kfac[1] = 0.5*(fy+kfac[0]);
   }
   else if (geo_dimension==3)
   {
      fy = ah[1]*ah[1]/(ah0_H*ah0_H);
      kfac[1] = fy;
      kfac[2] = 0.5*(kfac[0]+kfac[1]);
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
06/2007 WW Implementation
**************************************************************************/
void CMediumProperties::CalStressPermeabilityFactor3_Coef()
{
   int i;
   double d_max=0.0;
   double a[4], delta[4], Kni[4];
   double A[] = {-0.296, -0.1005, -0.1031, -0.1031 };
   double B[] = {-0.0056, -0.0073, -0.0074, -0.0074 };
   double C[] = {2.241, 1.0082, 1.135, 1.135 };
   double D[] = {-0.245, -0.23, -0.251, -0.251 };
   //
   double C1[] = {84.77, 44.37, 31.38, 20 };
   double C2[] = {0.02, 0.01, 0.01, 0.01 };
   //
   double JRC = c_coefficient[0];
   double JCS = c_coefficient[1];
   double UCS = c_coefficient[2];
   double a0 = 0.2*JRC*(0.2*UCS/JCS-0.1);
   //
   for(i=0; i<4; i++)
   {
      a[i] = a0;
      d_max = A[i]+B[i]*JRC+C[i]*pow(JCS/a0,D[i]);
      Kni[i] = 0.0178*JCS/a[i]+1.748*JRC-7.155;
      delta[i] = 0.01*(C1[i]-C2[i]*JCS/a0)*d_max;
      a0 -=  delta[i];
   }
   //
   //
   c_coefficient[1] = a[0];
   c_coefficient[2] = Kni[3];
   c_coefficient[3] = d_max;
   for(i=4; i<7; i++)
      c_coefficient[i] = delta[i-4];
   // Unit of stresses is MPa
   c_coefficient[7] *= 1.0e-6;
   c_coefficient[8] *= 1.0e-6;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2006 YD Implementation
last modification:
**************************************************************************/
CMediumPropertiesGroup* MMPGetGroup(const string &pcs_type_name)
{
   CMediumPropertiesGroup *m_mmp_group = NULL;
   list<CMediumPropertiesGroup*>::const_iterator p_mmp_group = mmp_group_list.begin();
   while(p_mmp_group!=mmp_group_list.end())
   {
      m_mmp_group = *p_mmp_group;
      if(m_mmp_group->pcs_type_name.compare(pcs_type_name)==0)
         return m_mmp_group;
      ++p_mmp_group;
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2006 YD Implementation
last modified:
**************************************************************************/
void MMPGroupDelete(/*string pcs_type_name*/)
{
   CMediumPropertiesGroup* m_mmp_group = NULL;
   list<CMediumPropertiesGroup*>::const_iterator p=mmp_group_list.begin();
   while (p!=mmp_group_list.end())
   {
      m_mmp_group = *p;
      // if(m_mmp_group->pcs_type_name.compare(pcs_type_name)==0){
      delete m_mmp_group;
      // mmp_group_list.remove(m_mmp_group);
      // return;
      ++p;
   }
   mmp_group_list.clear();
}


/**************************************************************************
FEMLib-Method:
08/2007 OK Implementation
**************************************************************************/
bool MMPExist(ifstream *mmp_file)
{
   string line_string;
   std::stringstream in;
   string name;
   //======================================================================
   while(!mmp_file->eof())
   {
      line_string = GetLineFromFile1(mmp_file);
      //--------------------------------------------------------------------
      //NAME
      if(line_string.find("$NAME")!=string::npos)
      {
         in.str(GetLineFromFile1(mmp_file));
         in >> name;                              //sub_line
         // Test whether MMP already exists
         for(int i=0;i<(int)mmp_vector.size();i++)
         {
            if(name.compare(mmp_vector[i]->name)==0)
            {
               return true;
            }
         }
         in.clear();
         continue;
      }
      //--------------------------------------------------------------------
   }
   //======================================================================
   return false;
}


/**************************************************************************
FEMLib-Method:
08/2007 OK Implementation
06/2009 OK Bug fix
**************************************************************************/
bool MMPExist()
{
   for(int i=0;i<(int)mmp_vector.size();i++)
   {
      for(int j=0;j<(int)mmp_vector.size();j++)
      {
         if(i==j)
            continue;
         if(mmp_vector[i]->name.compare(mmp_vector[j]->name)==0)
         {
            return true;
         }
      }
   }
   return false;
}


/**************************************************************************
FEMLib-Method:
Task: retrun the new permeability based on original permeability
      and old/new porosity
Programing:
11/2008 HS Implementation
last modification:
**************************************************************************/
double CMediumProperties::KozenyCarman(double k0, double n0, double n)
{
   double rt = 0.0;

   // TODO: here it should be the min_porosity and max_porosity instead of 0 and 1
   if (k0 < 1.0e-20 || k0 > 1.0 || n0 <=0 || n0 >= 1 || n <=0 || n >= 1 )
      return k0;
   else
   {
      rt = k0 ;

      rt *=pow( n / n0 , 3 );
   }

   return rt;
}


/**************************************************************************
FEMLib-Method:
Task: retrun the new permeability based on original permeability
      and old/new porosity (Koseny-Carman normalized)
Programing:
11/2008 HS Implementation
last modification:
**************************************************************************/
double CMediumProperties::KozenyCarman_normalized(double k0, double n0, double n)
{
   double rt = 0.0;

   // TODO: here it should be the min_porosity and max_porosity instead of 0 and 1
   if (k0 < 1.0e-20 || k0 > 1.0 || n0 <=0 || n0 >= 1 || n <=0 || n >= 1 )
      return k0;
   else
   {
      rt = k0 ;

      rt *=pow( n / n0 , 3 );

      rt *=pow( (1 - n0) / (1 - n) , 2);
   }

   return rt;
}


///////////////////////////////////////////////////////////////////////////
// old data structure functions //OK411

/**************************************************************************
//3.1 Subfunction of Porosity
 ROCKFLOW - Funktion: MATCalcSoilPorosityMethod1

 Aufgabe:
   Takes the porosity from a Geomechanical model, calulcates the effective stress
   and then takes the value of porosity from a curve.

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E SOIL_PROPERTIES *sp: Zeiger auf eine Instanz vom Typ SOIL_PROPERTIES.

Ergebnis:
- s.o. -
Programming Example
Programmaenderungen:
03/2004  CMCD First version

**************************************************************************/
double CMediumProperties::PorosityEffectiveStress(long index, double element_pressure)
{
   //OK411
   element_pressure = element_pressure;
   index = index;

   /*OK411
       int  nn, i, idummy, Type;
       long *element_nodes;
       double p;
      double znodes[8],ynodes[8],xnodes[8];
       double zelnodes[8],yelnodes[8],xelnodes[8];
      double coords[3];
      double Pie,angle;
      double sigma1,sigma2,sigma3;
      double a1,a2,a3,b1,b2,b3;
      double normx,normy,normz,normlen;
   double dircosl, dircosm, dircosn;
   double tot_norm_stress, eff_norm_stress;
   int material_group;
   double x_mid, y_mid, z_mid;

   // Normal stress calculated according to the orientation of the fracture element

   Type=ElGetElementType(index);
   material_group = ElGetElementGroupNumber(index);

   dircosl = dircosm = dircosn = 0.0;//Initialise variable

   if (Type == 2||Type == 3||Type == 4)  //Function defined for square, triangular and cubic elements
   {
   nn = ElNumberOfNodes[Type - 1];
   element_nodes = ElGetElementNodes(index);

   // Calculate directional cosins, note that this is currently set up
   // Sigma1 is in the y direction
   // Sigma2 is in the z direction
   // Sigma3 is in the x direction
   // This correspondes approximately to the KTB site conditions

   for (i=0;i<nn;i++)
   {
   zelnodes[i]=GetNodeZ(element_nodes[i]);
   yelnodes[i]=GetNodeY(element_nodes[i]);
   xelnodes[i]=GetNodeX(element_nodes[i]);
   }

   // Coordinate transformation to match sigma max direction
   // y direction matches the north direction
   Pie=3.141592654;
   angle=(porosity_model_values[3]*Pie)/180.;
   for (i=0;i<nn;i++)
   {
   znodes[i]=zelnodes[i];
   xnodes[i]=xelnodes[i]*cos(angle)-yelnodes[i]*sin(angle);
   ynodes[i]=xelnodes[i]*sin(angle)+yelnodes[i]*cos(angle);

   }

   if (Type == 2) //Square
   {
   a1=xnodes[2]-xnodes[0];
   a2=ynodes[2]-ynodes[0];
   a3=znodes[2]-znodes[0];
   b1=xnodes[3]-xnodes[1];
   b2=ynodes[3]-ynodes[1];
   b3=znodes[3]-znodes[1];
   normx=a2*b3-a3*b2;
   normy=a3*b1-a1*b3;
   normz=a1*b2-a2*b1;
   normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
   dircosl= normy/normlen;
   dircosm= normz/normlen;
   dircosn= normx/normlen;
   }

   if (Type == 3) //Cube
   {
   a1=xnodes[2]-xnodes[0];
   a2=ynodes[2]-ynodes[0];
   a3=znodes[2]-znodes[0];
   b1=xnodes[3]-xnodes[1];
   b2=ynodes[3]-ynodes[1];
   b3=znodes[3]-znodes[1];
   normx=a2*b3-a3*b2;
   normy=a3*b1-a1*b3;
   normz=a1*b2-a2*b1;
   normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
   dircosl= normy/normlen;
   dircosm= normz/normlen;
   dircosn= normx/normlen;
   }

   if (Type == 4) //Triangle
   {
   a1=xnodes[1]-xnodes[0];
   a2=ynodes[1]-ynodes[0];
   a3=znodes[1]-znodes[0];
   b1=xnodes[2]-xnodes[1];
   b2=ynodes[2]-ynodes[1];
   b3=znodes[2]-znodes[1];
   normx=a2*b3-a3*b2;
   normy=a3*b1-a1*b3;
   normz=a1*b2-a2*b1;
   normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
   dircosl= normy/normlen;
   dircosm= normz/normlen;
   dircosn= normx/normlen;
   }

   // Calculation of average location of element
   CalculateSimpleMiddelPointElement(index,coords);
   x_mid = coords[0];
   y_mid = coords[1];
   z_mid = coords[2];

   //Calculate fluid pressure in element
   p = element_pressure;

   //Calcualtion of stress according to Ito & Zoback 2000 for KTB hole

   sigma1=z_mid*0.045*-1e6;
   sigma2=z_mid*0.028*-1e6;
   sigma3=z_mid*0.02*-1e6;

   //Calculate total normal stress on element
   //Note in this case sigma2 corresponds to the vertical stress
   tot_norm_stress=sigma1*dircosl*dircosl+sigma2*dircosm*dircosm+sigma3*dircosn*dircosn;

   //Calculate normal effective stress
   eff_norm_stress=tot_norm_stress-p;

   //Take value of storage from a curve
   porosity=GetCurveValue((int) porosity_model_values[0], 0, eff_norm_stress, &idummy);

   }
   // default value if element type is not included in the method for calculating the normal
   else porosity=porosity_model_values[0];
   */
   return porosity;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 MX Implementation
last modification:
**************************************************************************/
double CMediumProperties::PorosityVolumetricFreeSwellingConstantIonicstrength(long index,double saturation,double temperature)
{
   //OK411
   index = index;
   saturation = saturation;
   temperature = temperature;

   /*OK411
     double mat_mp_m, beta;
     double porosity = 0.0;
     static double theta;
     static double porosity_n, porosity_IL, d_porosity, \
                   density_rock, fmon;
     static double S_0;
     static double epsilon;
     static double satu_0=0.20;
     //  static double porosity_min=0.05;
     static double ion_strength;
   static double F_const=96484.6, epsilon_0=8.854e-12;
   static double R=8.314510, psi=1.0;
   theta = 1.0;
   //--------------------------------------------------------------------
   // MMP medium properties
   S_0 =      porosity_model_values[1];      // Specific surface area [m^2/g]
   fmon =     porosity_model_values[2];     // Anteil quelfaehige mineral [-]
   mat_mp_m = porosity_model_values[3]; // Schichtenanzahl eines quellfähigen Partikels [-]
   beta =     porosity_model_values[6];     // modifications coefficient (z.B. free swelling Weimar beta=3.0)
   //--------------------------------------------------------------------
   // MSP solid properties
   CSolidProperties *m_msp = NULL;
   long group = ElGetElementGroupNumber(index);
   m_msp = msp_vector[group];
   density_rock  = m_msp->Density(1);
   //--------------------------------------------------------------------
   // State properties
   ion_strength = porosity_model_values[4]; // Ionic strength [M]
   satu_0 =       porosity_model_values[5];       // Initial saturation, if the sample is homogenous [-]
   //--------------------------------------------------------------------
   // Interlayer porosity calculation
   if (abs(temperature)>1.0e10) temperature=298.0;  //TODO MX
   epsilon = 87.0 + exp(-0.00456*(temperature-273.0));
   porosity_n = porosity_model_values[0];
   // Maximal inter layer porosity
   porosity_IL = fmon * psi * mat_mp_m * S_0 * (density_rock * 1.0e3) \
   * sqrt(epsilon * epsilon_0 * R * temperature / (2000.0 * F_const * F_const * ion_strength ));
   d_porosity=porosity_IL*(pow(saturation,beta)-pow(satu_0,beta));
   porosity_IL *=pow(saturation,beta);
   ElSetElementVal(index,PCSGetELEValueIndex("POROSITY_IL"),porosity_IL);
   // Total porosity calculation
   porosity = porosity_n+d_porosity;
   ElSetElementVal(index,PCSGetELEValueIndex("POROSITY"),porosity);
   */
   return porosity;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 MX Implementation
last modification:
**************************************************************************/
double CMediumProperties::PorosityEffectiveConstrainedSwelling(long index,double saturation,double temperature, double *porosity_sw)
{
   //OK411
   index = index;
   saturation = saturation;
   temperature = temperature;
   porosity_sw = porosity_sw;

   /*OK411
    // Soil Properties

     double mat_mp_m, beta;
     double porosity = 0.0;
     static double theta;
     double porosity_n, porosity_IL, d_porosity, \
                   n_total, density_rock, fmon;
   //  double n_max, n_min;
     static double S_0;
     // Component Properties
   static double  satu, epsilon;
   static double satu_0=0.20;
   static double porosity_min=0.05;
   static double ion_strength;
   static double F_const=96484.6, epsilon_0=8.854e-12;
   static double R=8.314510, psi=1.0;
   // Fluid properies
   int phase=1;

   //--------------------------------------------------------------------
   // MMP medium properties
   S_0 =      porosity_model_values[1];      // Specific surface area [m^2/g]
   fmon =     porosity_model_values[2];     // Anteil quelfaehige mineral [-]
   mat_mp_m = porosity_model_values[3]; // Schichtenanzahl eines quellfähigen Partikels [-]
   beta =     porosity_model_values[6];     // modifications coefficient (z.B. free swelling Weimar beta=3.0)
   //--------------------------------------------------------------------
   // MSP solid properties
   CSolidProperties *m_msp = NULL;
   //long group = ElGetElementGroupNumber(index);

   long group = m_pcs->m_msh->ele_vector[index]->GetPatchIndex();
   m_msp = msp_vector[group];
   density_rock  = m_msp->Density(1);
   //--------------------------------------------------------------------

   // Component Properties
   ion_strength = MATCalcIonicStrengthNew(index);
   if (ion_strength == 0.0){
   ion_strength = porosity_model_values[4]; //Ionic strength [M]
   }
   satu_0 = porosity_model_values[5];       //Initial saturation, if the sample is homogenous [-]
   porosity_min = porosity_model_values[7]; //minimal porosity after swelling compaction

   //  ion_strength = MATGetSoilIonicStrength(index);

   // Field variables
   //theta = GetNumericalTimeCollocation("TRANSPORT");
   theta = 1.0;
   // T = MATGetTemperatureGP (index,0.0,0.0,0.0,theta);
   //  T=298.0;
   phase=1;
   satu = saturation; //MATGetSaturationGP(phase, index, 0.0, 0.0, 0.0, theta); //only for fluid, phase=1

   //-----------------------------------------------------------------------
   // Interlayer Porositaet berechnen
   if (abs(temperature)>1.0e10) temperature=298.0;  //TODO MX
   epsilon =87.0+exp(-0.00456*(temperature-273.0));
   porosity_n = porosity_model_values[0];

   // Maximal inter layer porosity
   porosity_IL=fmon * psi * mat_mp_m * S_0 * (density_rock * 1.0e3) \
   * sqrt(epsilon * epsilon_0 * R * temperature / (2000.0 * F_const * F_const * ion_strength ));
   d_porosity=porosity_IL*(pow(satu, beta)-pow(satu_0, beta));

   //-----------Interlayer porosity calculation------------------

   //  porosity_IL = porosity_IL*satu;
   porosity_IL  *=pow(satu,beta);
   ElSetElementVal(index,PCSGetELEValueIndex("POROSITY_IL"),porosity_IL);
   // constrained swelling

   //-----------Effective porosity calculation------------------
   porosity = porosity_n-d_porosity;

   if (porosity<porosity_min)
   porosity =porosity_min;

   ElSetElementVal(index,PCSGetELEValueIndex("POROSITY"),porosity);
   //-----------Void ratio for swelling pressure calculation------------------
   //  e = porosity/(1.-porosity);
   //  ElSetElementVal(index,PCSGetELEValueIndex("VoidRatio"),e);
   //-----------Swelling potential calculation------------------
   // constrained swelling
   //  n_total=porosity_n - d_porosity;
   n_total=porosity_n + d_porosity-porosity_min;

   //  if(n_total>porosity_min)
   if(n_total>=0)
   {
   *porosity_sw = n_total;
   }
   else
   //    *porosity_sw=-porosity_IL*(satu-satu_0)+(porosity_n-porosity_min);
   *porosity_sw=d_porosity+(porosity_n-porosity_min);
   ElSetElementVal(index,PCSGetELEValueIndex("POROSITY_SW"),*porosity_sw);

   //-----------Swelling pressure calculation------------------
   // MATCalSoilSwellPressMethod0(index);
   */
   return porosity;
}


/**************************************************************************/
/* ROCKFLOW - Funktion: CalculateSoilPorosityMethod1
 */
/* Aufgabe:
   Berechnet die Porositaet in Abhaengigkeit von der Konzentration
   (Salzloesungsmodell)
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E SOIL_PROPERTIES *sp: Zeiger auf eine Instanz vom Typ SOIL_PROPERTIES.
 */
/* Ergebnis:
   - s.o. -
 */
/* Programmaenderungen:
   02/2000     RK     Erste Version
   11/2001     AH     Warnung entfernt

 */
/**************************************************************************/
/*
double CalculateSoilPorosityMethod1(SOIL_PROPERTIES * sp, long index)
{
    // double grainbound_limit = 1.73e-8;
    double porosity_limit = 0.01;
    int porosity_dependence_model;
    double val0, val1;
    static double porosity = 0.0;
    static double rho;
    static double theta;
    static double porosity_n, density_rock;
static double dissolution_rate, solubility_coefficient;
int timelevel = 1;

porosity_dependence_model = get_sp_porosity_dependence_model(sp);
switch (porosity_dependence_model) {
case 1: //Salzloesungsmodell
val0 = get_sp_porosity_model_field_value(sp, 0);
val1 = get_sp_porosity_model_field_value(sp, 1);

theta = GetNumericalTimeCollocation("TRANSPORT");
density_rock = GetSolidDensity(index);
dissolution_rate = GetTracerDissolutionRate(index, 0, 0);
porosity_n = GetElementSoilPorosity(index);
rho = GetFluidDensity(0, index, 0., 0., 0., theta);
if (GetTracerSolubilityModel(index, 0, 0) == 1) {
//if (GetRFProcessHeatReactModel())
solubility_coefficient = CalcTracerSolubilityCoefficient(index,0,0,1,PCSGetNODValueIndex("CONCENTRATION1",timelevel),1,PCSGetNODValueIndex("CONCENTRATION1",timelevel));
//else solubility_coefficient = CalcTracerSolubilityCoefficient(index,0,0,1,PCSGetNODValueIndex("CONCENTRATION1",timelevel),0,PCSGetNODValueIndex("CONCENTRATION1",timelevel));
} else {
solubility_coefficient = GetTracerSolubilityCoefficient(index, 0, 0);
}

porosity = porosity_n + 2 * porosity_n * dissolution_rate * val1 * rho
* (solubility_coefficient - val0) * dt / density_rock;

if (porosity > porosity_limit)
porosity = porosity_limit;

break;

default:
DisplayMsgLn("Unknown porosity dependence model!");
break;

}

return porosity;
}
*/

/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 MX Implementation
last modification:
**************************************************************************/
double CMediumProperties::PorosityVolumetricFreeSwelling(long index,double saturation,double temperature)
{
   //OK411
   index = index;
   saturation = saturation;
   temperature = temperature;

   /*OK411
     double mat_mp_m, beta;
     double porosity = 0.0;
     static double theta;
     static double porosity_n, porosity_IL, d_porosity, \
                   density_rock, fmon;
     static double S_0;
     static double epsilon;
     static double satu_0=0.20;
     //  static double porosity_min=0.05;
     static double ion_strength;
   static double F_const=96484.6, epsilon_0=8.854e-12;
   static double R=8.314510, psi=1.0;
   theta = 1.0;
   //--------------------------------------------------------------------
   // MMP medium properties
   S_0 =      porosity_model_values[1];      // Specific surface area [m^2/g]
   fmon =     porosity_model_values[2];     // Anteil quelfaehige mineral [-]
   mat_mp_m = porosity_model_values[3]; // Schichtenanzahl eines quellfähigen Partikels [-]
   beta =     porosity_model_values[6];     // modifications coefficient (z.B. free swelling Weimar beta=3.0)
   //--------------------------------------------------------------------
   // MSP solid properties
   CSolidProperties *m_msp = NULL;
   long group = ElGetElementGroupNumber(index);
   m_msp = msp_vector[group];
   density_rock  = m_msp->Density(1);
   //--------------------------------------------------------------------
   // State properties
   ion_strength = MATCalcIonicStrengthNew(index);
   if (ion_strength == 0.0){
   ion_strength = porosity_model_values[4]; // Ionic strength [M]
   }
   satu_0 =       porosity_model_values[5];       // Initial saturation, if the sample is homogenous [-]
   //--------------------------------------------------------------------
   // Interlayer porosity calculation
   epsilon = 87.0 + exp(-0.00456*(temperature-273.0));
   porosity_n = porosity_model_values[0];
   // Maximal inter layer porosity
   porosity_IL = fmon * psi * mat_mp_m * S_0 * (density_rock * 1.0e3) \
   * sqrt(epsilon * epsilon_0 * R * temperature / (2000.0 * F_const * F_const * ion_strength ));
   d_porosity=porosity_IL*(pow(saturation,beta)-pow(satu_0,beta));
   porosity_IL *=pow(saturation,beta);
   ElSetElementVal(index,PCSGetELEValueIndex("POROSITY_IL"),porosity_IL);
   // Total porosity calculation
   porosity = porosity_n+d_porosity;
   ElSetElementVal(index,PCSGetELEValueIndex("POROSITY"),porosity);
   */
   return porosity;
}


/**************************************************************************
FEMLib-Method: PorosityVolumetricChemicalReaction
Task: Porosity variation owing to chemical reactions
Programing:
05/2005 MX Implementation
last modification:
**************************************************************************/
double CMediumProperties::PorosityVolumetricChemicalReaction(long index)
{
   //OK411
   index = index;

   /*OK411
     int n=0, i, timelevel, nn=0, m=0;
     static long *element_nodes;
   //  long component;
     double porosity=0.0, tot_mineral_volume=0.0, tot_mineral_volume0=0.0;
     double conc[100];
     double mineral_volume[100],molar_volume[100];
   //  double *ergebnis=NULL;
    // REACTION_MODEL *rcml=NULL;

   //  rcml=GETReactionModel(index);
   REACT *rc = NULL; //SB
   rc->GetREACT();

   // MMP Medium Properties
   porosity = porosity_model_values[0];
   if (!rc) return porosity;

   //  m=rcml->number_of_equi_phases;
   m = rc->rcml_number_of_equi_phases;
   if (m ==0 || ElGetElement(index)==NULL) // wenn solid phases and Element existiert
   return porosity;

   // mineral molar volume abholen
   for (i=0; i<m; i++){
   molar_volume[i]=porosity_model_values[i+1];
   }

   tot_mineral_volume0=porosity_model_values[m+1]; //initial total volume of the minerals

   //  if (!rcml) {return porosity;}

   // calculate the concentrations of each solid phases (in mol/l soil) as element value from adjecent node values
   //  n=rcml->number_of_master_species;
   n = rc->rcml_number_of_master_species;
   //  n = get_rcml_number_of_master_species(rcml);
   for (int component=0; component<m+2+n+rc->rcml_number_of_ion_exchanges; component++) {
   conc[component] =0.0;
   // Not used: CompProperties *m_cp = cp_vec[component];
   //    int z_i = m_cp->valence;
   //	 m_cp->compname; //What's this for
   if (component>=n+2 && component<n+2+m){

   if (ElGetElementActiveState(index)){
   nn = ElGetElementNodesNumber(index);
   element_nodes = ElGetElementNodes(index);
   for (int j=0;j<nn;j++) {
   timelevel=1;
   conc[component] += PCSGetNODConcentration(element_nodes[j],component,timelevel);
   }
   conc[component] /= (double)nn;
   element_nodes = NULL;
   }

   // calculate the solid phase: volume =v_mi*Ci
   timelevel=0;
   //      conc[i]=CalcElementMeanConcentration (index, i, timelevel, ergebnis);
   mineral_volume[component-n-2] = conc[component]*molar_volume[component-n-2];
   tot_mineral_volume += mineral_volume[component-n-2];
   }
   } //for

   porosity += tot_mineral_volume0 - tot_mineral_volume;
   //  ElSetElementVal(index,"POROSITY",porosity);
   ElSetElementVal(index,PCSGetELEValueIndex("POROSITY"),porosity);
   */
   return porosity;
}


/**************************************************************************
FEMLib-Method: TortuosityFunction
Task:
Programing:
05/2007 PCH Diffusion tensor is handled by tortuosity tensor as in
      permeability. Agreed with GK
last modification:
**************************************************************************/
double CMediumProperties::TortuosityFunction(long number, double *gp, double theta,CFiniteElementStd* assem)
{
   //OK411
   theta = theta;
   gp = gp;
   number = number;
   assem = assem;

   // static int nidx0,nidx1;
   //----------------------------------------------------------------------
   /*OK411
     int count_nodes;
     long* element_nodes = NULL;
     double primary_variable[10];		//OK To Do
     int i;
     int no_pcs_names =(int)pcs_name_vector.size();
     for(i=0;i<no_pcs_names;i++){
       nidx0 = PCSGetNODValueIndex(pcs_name_vector[i],0);
       nidx1 = PCSGetNODValueIndex(pcs_name_vector[i],1);
       if(mode==0){ // Gauss point values
         primary_variable[i] = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
   + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
   }
   else if(mode==1){ // Node values
   primary_variable[i] = (1.-theta)*GetNodeVal(number,nidx0) \
   + theta*GetNodeVal(number,nidx1);
   }
   else if(mode==2){ // Element average value
   count_nodes = ElNumberOfNodes[ElGetElementType(number) - 1];
   element_nodes = ElGetElementNodes(number);
   for (i = 0; i < count_nodes; i++)
   primary_variable[i] += GetNodeVal(element_nodes[i],nidx1);
   primary_variable[i]/= count_nodes;
   }
   } //For Loop
   */
   switch (tortuosity_model)
   {
      case 0:                                     // Tortuosity is read from a curve
         //To do
         break;
      case 1:                                     // Constant Tortuosity
         tortuosity = tortuosity_model_values[0];
         break;
      default:
         DisplayMsgLn("Unknown tortuosisty model!");
         break;
   }
   return (tortuosity);
}


/**************************************************************************/
/* ROCKFLOW - Funktion: NonlinearFlowFunction
 */
/* Aufgabe:
   Berechnung der relativen Permeabilitaet
   fuer nichtlineare Fliessgesetze
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
 */
/* Ergebnis:
   0 bei Fehler, sonst 1 (dient nur zum Abbrechen der Funktion)
 */
/* Programmaenderungen:
   06/1999   OK   Implementierung
   09/2004   CMCD In GeoSys 4
 */
/**************************************************************************/
double CMediumProperties::NonlinearFlowFunction(long index, double *gp, double theta)
{
   double k_rel = 1.0;
   //OK411
   theta = theta;
   gp = gp;
   index = index;

   /*OK411
      //Pressure variable name (PRESSURE 1) not used as yet in this function
      int i;
      //Pointer to fluid properties
   //	int no_phases =(int)mfp_vector.size();
      CFluidProperties *m_mfp = NULL;
      int phase = (int)mfp_vector.size()-1;
      m_mfp = mfp_vector[phase];
       // Knotendaten
       int nn,Type;
       long *element_nodes;
   double p_element_node[4], h_element_node[4], z_element_node[4];
   // Materialdaten
   double alpha;
   double g, rho, mu;
   // Elementdaten
   double grad_h[2],grad_x,grad_y;
   double mult[2];
   double detjac, *invjac, invjac2d[4];
   double grad_omega[8];
   double grad_h_min = MKleinsteZahl;
   double xgt[3],ygt[3],zgt[3];				//CMCD Global x,y,z coordinates of traingular element
   double xt[3],yt[3];							//CMCD Local x,y coordinates of traingular element
   double pt_element_node[4], ht_element_node[4], zt_element_node[4]; //CMCD Pressure, depth head traingular element
   double dN_dx[3],dN_dy[3], area;				//CMCD Shape function derivates for triangles
   double isotropicgradient,porosity, Re, Lambda, apperture, hyd_radius, perm;
   double linear_q,turbulent_q,linear_grad, turbulent_grad, flow_rate, temp;
   double dircos[6];							//CMCD 04 2004
   gp[0]= gp[1] = gp[2] = 0.0;					//Gaus points, triangular interpretation value not relevant

   element_nodes = ElGetElementNodes(index);
   g = gravity_constant;
   //OK4104 theta = GetNumericalTimeCollocation("PRESSURE1");
   Type = ElGetElementType(index);
   //--------------------------------------------------------------------
   // MMP medium properties
   CMediumProperties *m_mmp = NULL;
   long group = ElGetElementGroupNumber(index);
   m_mmp = mmp_vector[group];
   //--------------------------------------------------------------------
   switch (flowlinearity_model){
   case 1://Element Type Dependent
   alpha = flowlinearity_model_values[0];
   rho = m_mfp->Density();

   switch (Type){
   case 1:
   nn = ElNumberOfNodes[0];
   for (i = 0; i < nn; i++) {
   p_element_node[i] = GetNodeVal(element_nodes[i],1);
   z_element_node[i] = GetNode(element_nodes[i])->z;
   h_element_node[i] = (p_element_node[i]) / (g * rho) + z_element_node[i];
   }
   invjac = GEOGetELEJacobianMatrix(index, &detjac);
   //          if(fabs(h_element_node[1]-h_element_node[0])>MKleinsteZahl)
   if (fabs(h_element_node[1] - h_element_node[0]) > grad_h_min) {
   k_rel = \
   pow(fabs(0.5 * sqrt(MSkalarprodukt(invjac, invjac, 3))), (alpha - 1.0)) * \
   pow(fabs(h_element_node[1] - h_element_node[0]), (alpha - 1.0));
   if (k_rel > 1)
   k_rel = 1.;
   } else
   k_rel = 1.0;

   break;

   case 2:
   nn = ElNumberOfNodes[1];        // Knotenanzahl nn muss 4 sein !
   for (i = 0; i < nn; i++) {
   p_element_node[i] = GetNodeVal(element_nodes[i],1);
   z_element_node[i] = GetNode(element_nodes[i])->z;
   h_element_node[i] = p_element_node[i] / (g * rho) + z_element_node[i];
   }
   Calc2DElementJacobiMatrix(index, 0.0, 0.0, invjac2d, &detjac);
   MGradOmega2D(grad_omega, 0, 0);         // Gradientenmatrix
   MMultMatVec(grad_omega, 2, 4, h_element_node, 4, mult, 2);
   MMultVecMat(mult, 2, invjac2d, 2, 2, grad_h, 2);
   //          if( (fabs(grad_h[0])>MKleinsteZahl)||(fabs(grad_h[1])>MKleinsteZahl) ) )
   if ((fabs(grad_h[0]) > grad_h_min) || (fabs(grad_h[1]) > grad_h_min)) {
   k_rel = \
   pow(fabs(sqrt((grad_h[0]) * (grad_h[0]) + (grad_h[1]) * (grad_h[1]))), \
   (alpha - 1.0));
   } else {
   k_rel = 1.0;
   }
   break;
   case 3:
   k_rel = 1.;
   break;
   }

   case 2://Equivalent Fractured Media represented by triangles   CMCD April 2004

   //Geometry
   nn = ElNumberOfNodes[Type - 1];
   element_nodes = ElGetElementNodes(index);

   for (i = 0; i < nn; i++)
   {
   xgt[i] = GetNodeX(element_nodes[i]);
   ygt[i] = GetNodeY(element_nodes[i]);
   zgt[i] = GetNodeZ(element_nodes[i]);
   }
   //Input parameters
   porosity = CMediumProperties::Porosity(index,theta);
   alpha = flowlinearity_model_values[0];
   apperture = porosity / flowlinearity_model_values[1]; //Change equivalent porosity to individual fracture porosity
   Re = flowlinearity_model_values[2];

   //Fluid properties
   rho = m_mfp->Density();
   mu  = m_mfp->Viscosity();
   Re = 0.0; //Reynolds number for turbulent flow CMCD 04. 2004
   Lambda = 0.0; //Frictional Resistence

   //Flow status
   hyd_radius = 2*apperture;
   perm = (apperture * apperture)/12.;
   linear_q = (Re*mu)/(hyd_radius*rho); //max linear q

   //Fluid pressure at each node
   for (i = 0; i < nn; i++) {
   pt_element_node[i] = GetNodeVal(element_nodes[i],1);
   zt_element_node[i] = GetNode(element_nodes[i])->z;
   ht_element_node[i] = pt_element_node[i] + ((g * rho) * zt_element_node[i]); //In Pascal
   }
   Calc2DElementCoordinatesTriangle(index,xt,yt,dircos); //CMCD included 03/2004
   area = ElGetElementVolume(index)/m_mmp->geo_area; //CMCD March 2004 removed wrong area in  code
   //Shape function derivatives
   dN_dx[0] = (yt[1] - yt[2]) / (2. * area);
   dN_dx[1] = (yt[2] - yt[0]) / (2. * area);
   dN_dx[2] = (yt[0] - yt[1]) / (2. * area);
   dN_dy[0] = (xt[2] - xt[1]) / (2. * area);
   dN_dy[1] = (xt[0] - xt[2]) / (2. * area);
   dN_dy[2] = (xt[1] - xt[0]) / (2. * area);
   grad_x = MSkalarprodukt(dN_dx, ht_element_node, nn);
   grad_y = MSkalarprodukt(dN_dy, ht_element_node, nn);
   //v2[0] = MSkalarprodukt(dN_dx, zg, nn);
   //v2[1] = MSkalarprodukt(dN_dy, zg, nn);
   //Assume isotropic nonlinear flow (p268 Kolditz 2001)
   linear_q = linear_q/3.0; // Here the whole element is considered hence 4* to remove avereaging effects
   isotropicgradient = pow((grad_x*grad_x+grad_y*grad_y),0.5);
   flow_rate = (perm * isotropicgradient)/mu;
   if (flow_rate > linear_q){
   turbulent_q = flow_rate-linear_q;
   linear_grad = (linear_q *apperture * mu)/perm;
   turbulent_grad = isotropicgradient - linear_grad;
   temp = pow((turbulent_grad/(rho*g)),1-alpha)/(turbulent_grad/(rho*g));
   k_rel = ((linear_grad*1.0) +(turbulent_grad*temp))/isotropicgradient;
   }
   else {
   k_rel = 1.0;
   }

   //velovec[0] = (-k_x * k_rel_grad_p * k_rel_S*k_rel / mu)* (v1[0] + (rho * g * v2[0]));
   //velovec[1] = (-k_y * k_rel_grad_p * k_rel_S*k_rel / mu) * (v1[1] + (rho * g * v2[1]));

   // special stop CMCD

   //		if (index == 7021){
   //			printf("\n");
   //			printf("Element 4516 k_rel = %g\n",k_rel);
   //			}
   //
   //		break;
   }
   */
   return k_rel;
}


/**************************************************************************
 ROCKFLOW - Funktion: Storage Function

 Aufgabe:
   Berechnet Speicherkoeffizienten

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long index: Elementnummer

 Ergebnis:
   rel. Permeabilitaet

Programmaenderungen:
1/2001   C.Thorenz   Erste Version
7/2007   C.Thorenz   Div. Druckabhaengigkeiten
06/2003 OK/CM case 4: Storage as function of effective stress read from curve
11/2003   CMCD		Geothermal methods added (20)

Case overview
0		Curve
1		Constant
2		Funktion der effektiven Spannung und des Drucks in Elementmitte
3		Funktion der effektiven Spannung und des Drucks in Elementmitte ueber Kurve
4		Storage as function of effective stress read from curve
5     Storage as normal stress in element in stress field defined by KTB stress field.
6		Storage as normal stress in element in stress field defined by KTB stress
field, function to increase storage with distance from borehole.
**************************************************************************/
double CMediumProperties::StorageFunction(long index,double *gp,double theta)
{
   //OK411
   theta = theta;
   gp = gp;
   index = index;

   //int nn, i, Type;
   //int idummy;
   //double p; //WW, sigma, z[8];
   //int phase;
   //double density_solid, stress_eff,S;
   //double coords[3];
   //double znodes[8],ynodes[8],xnodes[8];
   //double zelnodes[8],yelnodes[8],xelnodes[8];
   //double Pie,angle;
   //double sigma1,sigma2,sigma3;
   //double a1,a2,a3,b1,b2,b3;
   //double normx,normy,normz,normlen;
   double dircosl, dircosm, dircosn;
   //double tot_norm_stress, eff_norm_stress;
   //int material_group;
   //double x_mid, y_mid, z_mid, x_bore, y_bore, z_bore, distance;
   dircosl = dircosm = dircosn = 0.0;             //Initialise variable
   // static int nidx0,nidx1;
   //double primary_variable[10];		//OK To Do
   //int count_nodes;
   /*OK411
      long* element_nodes = NULL;
      int no_pcs_names =(int)pcs_name_vector.size();
      for(i=0;i<no_pcs_names;i++)
       {
        nidx0 = PCSGetNODValueIndex(pcs_name_vector[i],0);
        nidx1 = PCSGetNODValueIndex(pcs_name_vector[i],1);
        if(mode==0) // Gauss point values
         {
          primary_variable[i] = (1.-theta)*InterpolValue(index,nidx0,gp[0],gp[1],gp[2]) \
                        + theta*InterpolValue(index,nidx1,gp[0],gp[1],gp[2]);		}
   else if(mode==1) // Node values
   {
   primary_variable[i] = (1.-theta)*GetNodeVal(index,nidx0) \
   + theta*GetNodeVal(index,nidx1);
   }
   else if(mode==2) // Element average value
   {
   //MX		count_nodes = ElNumberOfNodes[ElGetElementType(number) - 1];
   count_nodes = ElNumberOfNodes[ElGetElementType(index) - 1];
   //MX		element_nodes = ElGetElementNodes(number);
   element_nodes = ElGetElementNodes(index);
   for (i = 0; i < count_nodes; i++)
   primary_variable[i] += GetNodeVal(element_nodes[i],nidx1);
   primary_variable[i]/= count_nodes;
   }
   }
   */
   switch (storage_model)
   {

      case 0:

         //OK411 storage = GetCurveValue((int) storage_model_values[0], 0, primary_variable[0], &idummy);

         break;

      case 1:
         // Konstanter Wert
         storage = storage_model_values[0];
         break;

      case 2:
         // Funktion der effektiven Spannung und des Drucks in Elementmitte
#ifdef obsolete                             //WW. 06.11.2008
         // Den Druck holen
         p = primary_variable[0];

         /* Mittlere Tiefe */
         nn = ElNumberOfNodes[ElGetElementType(index) - 1];
         element_nodes = ElGetElementNodes(index);

         for (i = 0; i < nn; i++)
            z[i] = GetNodeZ(element_nodes[i]);

         /* Spannung = sigma(z0) + d_sigma/d_z*z */
         //OKsigma = storage_model_values[2] + storage_model_values[3] * InterpolValueVector(index, z, 0., 0., 0.);
         sigma = storage_model_values[2] + storage_model_values[3] * InterpolValueVector(ElGetElementType(index), z, 0., 0., 0.);

         /* Auf effektive Spannung umrechnen */
         sigma -= p;

         storage = exp(storage_model_values[0] - storage_model_values[1] * log(sigma));
#endif                                      //#ifdef obsolete //WW. 06.11.2008
         break;

      case 3:
         /* Funktion der effektiven Spannung und des Drucks in Elementmitte
                    ueber Kurve */
#ifdef obsolete                             //WW. 06.11.2008
         /* Den Druck holen */
         p = primary_variable[0];

         /* Mittlere Tiefe */
         nn = ElNumberOfNodes[ElGetElementType(index) - 1];
         element_nodes = ElGetElementNodes(index);

         for (i = 0; i < nn; i++)
            z[i] = GetNodeZ(element_nodes[i]);

         /* Spannung = sigma(z0) + d_sigma/d_z*z */
         sigma = storage_model_values[1] + storage_model_values[2] * InterpolValueVector(ElGetElementType(index), z, 0., 0., 0.);

         /* Auf effektive Spannung umrechnen */
         sigma -= p;

         storage = GetCurveValue((int)storage_model_values[0],0,sigma,&i);
#endif
         break;

      case 4:                                     /* McD Storage as function of effective stress read from curve */
         /*OK411
               CalculateSimpleMiddelPointElement(index, coords);
                 p = primary_variable[0];
                density_solid = storage_model_values[2];
               stress_eff = (fabs(coords[2])*gravity_constant*density_solid) - p;
                 storage =GetCurveValue((int) storage_model_values[0], 0, stress_eff, &idummy);
               break;
         */
      case 5:                                     /* Stroage : Normal stress calculated according to the orientation of the fracture element*/
         /*OK411
               Type=ElGetElementType(index);
               material_group = ElGetElementGroupNumber(index);

                  if (Type == 2||Type == 3||Type == 4)  //Function defined for square, triangular and cubic elements
                  {
                     nn = ElNumberOfNodes[Type - 1];
                     element_nodes = ElGetElementNodes(index);

                  for (i=0;i<nn;i++)
                  {
         zelnodes[i]=GetNodeZ(element_nodes[i]);
         yelnodes[i]=GetNodeY(element_nodes[i]);
         xelnodes[i]=GetNodeX(element_nodes[i]);
         }

         Pie=3.141592654;
         angle=(storage_model_values[3]*Pie)/180.;
         for (i=0;i<nn;i++)
         {
         znodes[i]=zelnodes[i];
         xnodes[i]=xelnodes[i]*cos(angle)-yelnodes[i]*sin(angle);
         ynodes[i]=xelnodes[i]*sin(angle)+yelnodes[i]*cos(angle);

         }

         if (Type == 2) //Square
         {
         a1=xnodes[2]-xnodes[0];
         a2=ynodes[2]-ynodes[0];
         a3=znodes[2]-znodes[0];
         b1=xnodes[3]-xnodes[1];
         b2=ynodes[3]-ynodes[1];
         b3=znodes[3]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
         }

         if (Type == 3) //Cube
         {
         a1=xnodes[2]-xnodes[0];
         a2=ynodes[2]-ynodes[0];
         a3=znodes[2]-znodes[0];
         b1=xnodes[3]-xnodes[1];
         b2=ynodes[3]-ynodes[1];
         b3=znodes[3]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
         }

         if (Type == 4) //Triangle
         {
         a1=xnodes[1]-xnodes[0];
         a2=ynodes[1]-ynodes[0];
         a3=znodes[1]-znodes[0];
         b1=xnodes[2]-xnodes[1];
         b2=ynodes[2]-ynodes[1];
         b3=znodes[2]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
         }
         // Calculate average location of the element
         CalculateSimpleMiddelPointElement(index, coords);
         x_mid = coords[0];
         y_mid = coords[1];
         z_mid = coords[2];

         p = primary_variable[0];

         //Calcualtion of stress according to Ito & Zoback 2000 for KTB hole

         sigma1=z_mid*0.045*-1e6;
         sigma2=z_mid*0.028*-1e6;
         sigma3=z_mid*0.02*-1e6;

         ///Calculate total normal stress on element
         //Note in this case sigma2 corresponds to the vertical stress
         tot_norm_stress=sigma1*dircosl*dircosl+sigma2*dircosm*dircosm+sigma3*dircosn*dircosn;

         //Calculate normal effective stress
         eff_norm_stress=tot_norm_stress-p;

         // special stop CMCD
         if (eff_norm_stress>220000000.){
         phase=0;
         }

         //Take value of storage from a curve
         S=GetCurveValue((int) storage_model_values[0], 0, eff_norm_stress, &idummy);
         }

         // default value if element type is not included in the method for calculating the normal
         else S=storage_model_values[2];
         storage = S;
         break;
         */
      case 6:                                     /* Normal stress calculated according to the orientation of the fracture element*/
         /*OK411
               Type=ElGetElementType(index);
               material_group = ElGetElementGroupNumber(index);

               if (material_group == 0)
               {
                  if (Type == 2||Type == 3||Type == 4)
                  {
                     nn = ElNumberOfNodes[Type - 1];
                     element_nodes = ElGetElementNodes(index);

         for (i=0;i<nn;i++)
         {
         zelnodes[i]=GetNodeZ(element_nodes[i]);
         yelnodes[i]=GetNodeY(element_nodes[i]);
         xelnodes[i]=GetNodeX(element_nodes[i]);
         }

         Pie=3.141592654;
         angle=(storage_model_values[3]*Pie)/180.;
         for (i=0;i<nn;i++)
         {
         znodes[i]=zelnodes[i];
         xnodes[i]=xelnodes[i]*cos(angle)-yelnodes[i]*sin(angle);
         ynodes[i]=xelnodes[i]*sin(angle)+yelnodes[i]*cos(angle);

         }

         if (Type == 2)
         {
         a1=xnodes[2]-xnodes[0];
         a2=ynodes[2]-ynodes[0];
         a3=znodes[2]-znodes[0];
         b1=xnodes[3]-xnodes[1];
         b2=ynodes[3]-ynodes[1];
         b3=znodes[3]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
         }

         if (Type == 3)
         {
         a1=xnodes[2]-xnodes[0];
         a2=ynodes[2]-ynodes[0];
         a3=znodes[2]-znodes[0];
         b1=xnodes[3]-xnodes[1];
         b2=ynodes[3]-ynodes[1];
         b3=znodes[3]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
         }

         if (Type == 4)
         {
         a1=xnodes[1]-xnodes[0];
         a2=ynodes[1]-ynodes[0];
         a3=znodes[1]-znodes[0];
         b1=xnodes[2]-xnodes[1];
         b2=ynodes[2]-ynodes[1];
         b3=znodes[2]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
         }

         CalculateSimpleMiddelPointElement(index,coords);
         x_mid = coords[0];
         y_mid = coords[1];
         z_mid = coords[2];

         p = primary_variable[0];

         x_bore = storage_model_values[4];
         y_bore = storage_model_values[5];
         z_bore = storage_model_values[6];

         distance = pow(pow((x_mid-x_bore),2.0) + pow((y_mid-y_bore),2.0) + pow((z_mid-z_bore),2.0),0.5);

         sigma1=z_mid*0.045*-1e6;
         sigma2=z_mid*0.028*-1e6;
         sigma3=z_mid*0.02*-1e6;

         ///Calculate total normal stress on element
         tot_norm_stress=sigma1*dircosl*dircosl+sigma2*dircosm*dircosm+sigma3*dircosn*dircosn;

         eff_norm_stress=tot_norm_stress-p;

         S=GetCurveValue((int) storage_model_values[0], 0, eff_norm_stress, &idummy);
         if (distance > storage_model_values[7]){
         distance = storage_model_values[7];
         }
         if (distance > 2) S = S + (S * (distance-2) * storage_model_values[8]);
         }

         else S=storage_model_values[2];
         }

         else S=storage_model_values[2];
         storage = S;
         */
         break;
      case 7:                                     // poroelasticity RW
         storage = storage_model_values[1];
         break;
      default:
         storage = 0.0;                           //OK DisplayMsgLn("The requested storativity model is unknown!!!");
         break;
   }
   return storage;
}


double CMediumProperties::PermeabilityPressureFunction(long index,double *gp,double theta)
{
   double k_rel=0.0;
   //OK411
   theta = theta;
   gp = gp;
   index = index;

#ifdef obsolete                                //OK411
   int nn, i, idummy,p_idx1;
   long *element_nodes;
   static double gh,  p, eins_durch_rho_g, sigma, z[8], h[8], grad_h[3];
   static double invjac[8], detjac, grad_omega[8];
   static double mult[3];
   double x_mid,y_mid,z_mid,coords[3];
   double density_solid, stress_eff;
   //	double porosity, factora, factorb,valuelogk;
   CFluidProperties* m_mfp=NULL;

   //Collect primary variables
   static int nidx0,nidx1;
   double primary_variable[10];                   //OK To Do
   int count_nodes;
   int no_pcs_names =(int)pcs_name_vector.size();
   for(i=0;i<no_pcs_names;i++)
   {
      nidx0 = PCSGetNODValueIndex(pcs_name_vector[i],0);
      nidx1 = PCSGetNODValueIndex(pcs_name_vector[i],1);
      if(mode==0)                                 // Gauss point values
      {
         primary_variable[i] = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
            + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
      }
      else if(mode==1)                            // Node values
      {
         primary_variable[i] = (1.-theta)*GetNodeVal(number,nidx0) \
            + theta*GetNodeVal(number,nidx1);
      }
      else if(mode==2)                            // Element average value
      {
         count_nodes = ElNumberOfNodes[ElGetElementType(number) - 1];
         element_nodes = ElGetElementNodes(number);
         for (i = 0; i < count_nodes; i++)
            primary_variable[i] += GetNodeVal(element_nodes[i],nidx1);
         primary_variable[i]/= count_nodes;
      }
   }
   switch (permeability_pressure_model)
   {
      case 0:                                     //Curve function
         k_rel = 1.0;
         break;
      case 1:                                     //No functional dependence
         k_rel = 1.0;
         break;
         /* Funktion der effektiven Spannung */
      case 2:
#ifdef obsolete                             //WW. 06.11.2008
         /* Den Druck holen */
         p = primary_variable[0];
         /* Mittlere Tiefe */
         nn = ElNumberOfNodes[ElGetElementType(index) - 1];
         element_nodes = ElGetElementNodes(index);
         for (i = 0; i < nn; i++)
            z[i] = GetNodeZ(element_nodes[i]);
         /* Spannung = sigma(z0) + d_sigma/d_z*z */
         sigma = permeability_pressure_model_values[2] + permeability_pressure_model_values[3] * InterpolValueVector(ElGetElementType(index), z, 0., 0., 0.);
         /* Auf effektive Spannung umrechnen */
         sigma -= p;
         k_rel = exp(permeability_pressure_model_values[0] - permeability_pressure_model_values[1] * log(sigma));
#endif
         break;
      case 3:                                     /* Turbulentes Fliessen */
         /* k_rel = max(min((grad(h)*alpha1)^(alpha2-1), alpha4),alpha3) */
         m_mfp = MFPGet("LIQUID");
                                                  // YDGetFluidDensity(0, index, 0., 0., 0., 1.);
         eins_durch_rho_g = 1./gravity_constant/m_mfp->Density();
         nn = ElNumberOfNodes[ElGetElementType(index) - 1];
         element_nodes = ElGetElementNodes(index);
         p_idx1 = PCSGetNODValueIndex("PRESSURE1",1);
         for (i=0;i<nn;i++)
            h[i] = GetNodeVal(element_nodes[i], p_idx1)*eins_durch_rho_g
               + GetNodeZ(element_nodes[i]);
         switch (ElGetElementType(index))
         {
            default:
               DisplayMsgLn("Error in GetSoilRelPermPress!");
               DisplayMsgLn("  Nonlinear permeability not available!");
               abort();
            case 2:
               Calc2DElementJacobiMatrix(index,0.0,0.0,invjac,&detjac);
               MGradOmega2D(grad_omega,0,0);      /* Gradientenmatrix */
               MMultMatVec(grad_omega,2,4,h,4,mult,2);
               MMultVecMat(mult,2,invjac,2,2,grad_h,2);
               gh = sqrt(grad_h[0]*grad_h[0]+grad_h[1]*grad_h[1])*permeability_pressure_model_values[0];
               if(gh<MKleinsteZahl)
                  k_rel = 1.0;
               else
                  k_rel = max(min(pow(gh,permeability_pressure_model_values[1]-1.0),permeability_pressure_model_values[3]),permeability_pressure_model_values[2]);
         }
      case 4:
#ifdef obsolete                             //WW. 06.11.2008
         /* Funktion der effektiven Spannung ueber Kurve */
         /* Den Druck holen */
         p = primary_variable[0];
         /* Mittlere Tiefe */
         nn = ElNumberOfNodes[ElGetElementType(index) - 1];
         element_nodes = ElGetElementNodes(index);
         for (i = 0; i < nn; i++)
            z[i] = GetNodeZ(element_nodes[i]);
         /* Spannung = sigma(z0) + d_sigma/d_z*z */
         sigma = permeability_pressure_model_values[1] + permeability_pressure_model_values[2] * InterpolValueVector(ElGetElementType(index), z, 0., 0., 0.);
         /* Auf effektive Spannung umrechnen */
         sigma -= p;
         k_rel = GetCurveValue((int)permeability_pressure_model_values[0], 0, sigma, &i);
#endif
         break;
      case 5:
         /* Funktion der effektiven Spannung ueber Kurve CMCD 26.06.2003*/
         /* Average depth */
         CalculateSimpleMiddelPointElement(index,coords);
         x_mid = coords[0];
         y_mid = coords[1];
         z_mid = coords[2];
         p = primary_variable[0];
         density_solid = permeability_pressure_model_values[2];
         stress_eff = (fabs(z_mid)*gravity_constant*density_solid) - p;
         k_rel = GetCurveValue((int) permeability_pressure_model_values[0], 0, stress_eff, &idummy);
         break;
      case 6:
         k_rel = PermeabilityPressureFunctionMethod1(index,primary_variable[0]);
         break;
      case 7:
         k_rel = PermeabilityPressureFunctionMethod2(index,primary_variable[0]);
         break;
      case 8:
         k_rel = PermeabilityPressureFunctionMethod3(index,primary_variable[0]);
         break;
      case 9:
         k_rel = PermeabilityPressureFunctionMethod4(index,primary_variable[0], primary_variable[1]);
         break;
      default:                                    // CMCD Einbau
         k_rel = 1.0;                             // CMCD Einbau
         break;                                   // CMCD Einbau
   }
#endif
   return k_rel;
}


/**************************************************************************
12.(ii)a Subfunction of Permeability_Function_Pressure
ROCKFLOW - Funktion: MATCalcPressurePermeabilityMethod1
Application to KTB
 Aufgabe:
   Calculates relative permeability from
   the normal stress according to orientation of the fractures in the KTB site system
   converts the normal stress to effective stress, reads from a curve what the permeability
   of a fracture under the given effective stress is.

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
E: long index: Elementnummer, double primary_variable pressure
R: relative permeability

Ergebnis:
rel. Permeabilitaet

Programmaenderungen:
09/2004   CMCD Inclusion in GeoSys 4
**************************************************************************/
double CMediumProperties::PermeabilityPressureFunctionMethod1(long index,double pressure)
{
   //OK411
   index = index;
   pressure = pressure;

#ifdef obsolete                                //OK411
   double R,Pie,p;
   double znodes[8],ynodes[8],xnodes[8], angle;
   double x_mid, y_mid, z_mid, coords[3];
   double zelnodes[8],yelnodes[8],xelnodes[8];
   double sigma1,sigma2,sigma3;
   double a1,a2,a3,b1,b2,b3;
   double normx,normy,normz,normlen;
   double dircosl, dircosm, dircosn;
   double tot_norm_stress, eff_norm_stress;
   int material_group,Type;
   int phase,nn,i,idummy;
   long *element_nodes;
   dircosl = dircosm = dircosn =0.0;
   Type=ElGetElementType(index);
   material_group = ElGetElementGroupNumber(index);

   /* special stop CMCD*/
   if (index == 6590)
   {
      phase=0;
      /* Stop here*/
   }

   if (Type == 2||Type == 3||Type == 4)           /*Function defined for square, triangular and cubic elements*/
   {
      nn = ElNumberOfNodes[Type - 1];
      element_nodes = ElGetElementNodes(index);

      /* Calculate directional cosins, note that this is currently set up*/
      /* Sigma1 is in the y direction*/
      /* Sigma2 is in the z direction*/
      /* Sigma3 is in the x direction*/
      /* This correspondes approximately to the KTB site conditions*/

      for (i=0;i<nn;i++)
      {
         zelnodes[i]=GetNodeZ(element_nodes[i]);
         yelnodes[i]=GetNodeY(element_nodes[i]);
         xelnodes[i]=GetNodeX(element_nodes[i]);
      }

      /*Coordinate transformation to match sigma max direction*/
      /* y direction matches the north direction*/

      Pie=3.141592654;
      for (i=0;i<nn;i++)
      {
         znodes[i]=zelnodes[i];
         angle=(permeability_pressure_model_values[3]*Pie)/180.;
         xnodes[i]=xelnodes[i]*cos(angle)-yelnodes[i]*sin(angle);
         ynodes[i]=xelnodes[i]*sin(angle)+yelnodes[i]*cos(angle);

      }

      if (Type == 2)                              /*Square*/
      {
         a1=xnodes[2]-xnodes[0];
         a2=ynodes[2]-ynodes[0];
         a3=znodes[2]-znodes[0];
         b1=xnodes[3]-xnodes[1];
         b2=ynodes[3]-ynodes[1];
         b3=znodes[3]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
      }

      if (Type == 3)                              /*Cube*/
      {
         a1=xnodes[2]-xnodes[0];
         a2=ynodes[2]-ynodes[0];
         a3=znodes[2]-znodes[0];
         b1=xnodes[3]-xnodes[1];
         b2=ynodes[3]-ynodes[1];
         b3=znodes[3]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
      }

      if (Type == 4)                              /*Triangle*/
      {
         a1=xnodes[1]-xnodes[0];
         a2=ynodes[1]-ynodes[0];
         a3=znodes[1]-znodes[0];
         b1=xnodes[2]-xnodes[1];
         b2=ynodes[2]-ynodes[1];
         b3=znodes[2]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
      }

      /* Calculation of average location of element*/
      CalculateSimpleMiddelPointElement(index,coords);
      x_mid = coords[0];
      y_mid = coords[1];
      z_mid = coords[2];

      /*Calculate fluid pressure in element*/
      p = pressure;
      /*Calcualtion of stress according to Ito & Zoback 2000 for KTB hole*/

      sigma1=z_mid*0.045*-1e6;
      sigma2=z_mid*0.028*-1e6;
      sigma3=z_mid*0.02*-1e6;

      /*Calculate total normal stress on element
      Note in this case sigma2 corresponds to the vertical stress*/
      tot_norm_stress=sigma1*dircosl*dircosl+sigma2*dircosm*dircosm+sigma3*dircosn*dircosn;

      /*Calculate normal effective stress*/
      eff_norm_stress=tot_norm_stress-p;

      /*Take value of storage from a curve*/
      R=GetCurveValue((int) permeability_pressure_model_values[0], 0, eff_norm_stress, &idummy);

      return R;

   }
   /* default value if element type is not included in the method for calculating the normal */
   else R=permeability_pressure_model_values[2];
#endif
   return R;
}


/**************************************************************************
 12.(ii)b Subfunction of Permeability_Function_Pressure
 Function: MATCalcPressurePermeabilityMethod2
 Application to KTB
 Aufgabe:
   Calculates relative permeability from
   the normal stress according to orientation of the fractures in the KTB site system
   converts the normal stress to effective stress, reads from a curve what the permeability
   of a fracture under the given effective stress is. Permeability is then related to
   the distance from the borehole.

Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
E: long index: Elementnummer, double primary_variable pressure
R: relative permeability

Ergebnis:
rel. Permeabilitaet

Programmaenderungen:
09/2004   CMCD Inclusion in GeoSys 4
**************************************************************************/
double CMediumProperties::PermeabilityPressureFunctionMethod2(long index,double pressure)
{
   //OK411
   index = index;
   pressure = pressure;

#ifdef obsolete                                //OK411
   double R,Pie,p;
   double znodes[8],ynodes[8],xnodes[8], angle;
   double x_mid, y_mid, z_mid, coords[3];
   double zelnodes[8],yelnodes[8],xelnodes[8];
   double sigma1,sigma2,sigma3;
   double a1,a2,a3,b1,b2,b3;
   double normx,normy,normz,normlen;
   double dircosl, dircosm, dircosn;
   double tot_norm_stress, eff_norm_stress;
   int material_group,Type;
   double x_bore, y_bore, z_bore, distance;
   int nn,i,idummy;
   long *element_nodes;
   dircosl = dircosm = dircosn =0.0;
   Type=ElGetElementType(index);
   material_group = ElGetElementGroupNumber(index);

   /* special stop CMCD*/
   if (index == 4516)
   {
      /* Stop here */
   }

   if (Type == 2||Type == 3||Type == 4)           /*Function defined for square, triangular and cubic elements*/
   {
      nn = ElNumberOfNodes[Type - 1];
      element_nodes = ElGetElementNodes(index);

      /* Calculate directional cosins, note that this is currently set up*/
      /* Sigma1 is in the y direction*/
      /* Sigma2 is in the z direction*/
      /* Sigma3 is in the x direction*/
      /* This correspondes approximately to the KTB site conditions*/

      for (i=0;i<nn;i++)
      {
         zelnodes[i]=GetNodeZ(element_nodes[i]);
         yelnodes[i]=GetNodeY(element_nodes[i]);
         xelnodes[i]=GetNodeX(element_nodes[i]);
      }

      Pie=3.141592654;
      angle=(permeability_pressure_model_values[3]*Pie)/180.;
      for (i=0;i<nn;i++)
      {
         znodes[i]=zelnodes[i];
         xnodes[i]=xelnodes[i]*cos(angle)-yelnodes[i]*sin(angle);
         ynodes[i]=xelnodes[i]*sin(angle)+yelnodes[i]*cos(angle);

      }

      if (Type == 2)                              /*Square*/
      {
         a1=xnodes[2]-xnodes[0];
         a2=ynodes[2]-ynodes[0];
         a3=znodes[2]-znodes[0];
         b1=xnodes[3]-xnodes[1];
         b2=ynodes[3]-ynodes[1];
         b3=znodes[3]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
      }

      if (Type == 3)                              /*Cube*/
      {
         a1=xnodes[2]-xnodes[0];
         a2=ynodes[2]-ynodes[0];
         a3=znodes[2]-znodes[0];
         b1=xnodes[3]-xnodes[1];
         b2=ynodes[3]-ynodes[1];
         b3=znodes[3]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
      }

      if (Type == 4)                              /*Triangle*/
      {
         a1=xnodes[1]-xnodes[0];
         a2=ynodes[1]-ynodes[0];
         a3=znodes[1]-znodes[0];
         b1=xnodes[2]-xnodes[1];
         b2=ynodes[2]-ynodes[1];
         b3=znodes[2]-znodes[1];
         normx=a2*b3-a3*b2;
         normy=a3*b1-a1*b3;
         normz=a1*b2-a2*b1;
         normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
         dircosl= normy/normlen;
         dircosm= normz/normlen;
         dircosn= normx/normlen;
      }

      CalculateSimpleMiddelPointElement(index,coords);
      x_mid = coords[0];
      y_mid = coords[1];
      z_mid = coords[2];

      /*Calculate fluid pressure in element*/
      p = pressure;

      /*Calculate distance from borehole*/
      x_bore = permeability_pressure_model_values[4];
      y_bore = permeability_pressure_model_values[5];
      z_bore = permeability_pressure_model_values[6];

      distance = pow(pow((x_mid-x_bore),2.0) + pow((y_mid-y_bore),2.0) + pow((z_mid-z_bore),2.0),0.5);

      /*Calcualtion of stress according to Ito & Zoback 2000 for KTB hole*/

      sigma1=z_mid*0.045*-1e6;
      sigma2=z_mid*0.028*-1e6;
      sigma3=z_mid*0.02*-1e6;

      ///Calculate total normal stress on element
      /*Note in this case sigma2 corresponds to the vertical stress*/
      tot_norm_stress=sigma1*dircosl*dircosl+sigma2*dircosm*dircosm+sigma3*dircosn*dircosn;

      /*Calculate normal effective stress*/
      eff_norm_stress=tot_norm_stress-p;

      /*Take value of storage from a curve*/
      R=GetCurveValue((int) permeability_pressure_model_values[0], 0, eff_norm_stress, &idummy);
      if (distance > permeability_pressure_model_values[7])
      {
         distance = permeability_pressure_model_values[7];

      }
      if (distance > 2) R = R + (R * (distance-2) * permeability_pressure_model_values[8]);

      return R;
   }
   /* default value if element type is not included in the method for calculating the normal */
   else R=permeability_pressure_model_values[2];
#endif
   return R;
}


/**************************************************************************
 12.(ii)c Subfunction of Permeability_Function_Pressure
 Function: MATCalcPressurePermeabilityMethod3
 Application to Urach

 Aufgabe:
   The normal stress across the fractures is calculated according to an approximate formulation
   from the relationship of stress with depth. This normal stress is then converted into a permeabilty
   by reference to effective stress and a curve.

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
E: long index: Elementnummer, double primary_variable pressure
R: relative permeability

Ergebnis:
rel. Permeabilitaet

Programmaenderungen:
09/2004   CMCD Inclusion in GeoSys 4
**************************************************************************/
double CMediumProperties::PermeabilityPressureFunctionMethod3(long index,double pressure)
{
   //OK411
   index = index;
   pressure = pressure;

   double perm=0.0;
#ifdef obsolete                                //OK411
   /* Funktion der effektiven Spannung ueber Kurve CMCD 26.06.2003*/
   double x_mid, y_mid, z_mid, coords[3];
   double stress_eff,p,perm;
   int idummy;
   /* Average depth */
   CalculateSimpleMiddelPointElement(index,coords);
   x_mid = coords[0];
   y_mid = coords[1];
   z_mid = coords[2];

   /*Calculate fluid pressure in element*/
   p = pressure;
   stress_eff = (fabs(z_mid)*permeability_pressure_model_values[1]*1e6)-p;
   perm = GetCurveValue((int) permeability_pressure_model_values[0], 0, stress_eff, &idummy);
#endif
   return perm;
}


/**************************************************************************
12.(ii)d Subfunction of Permeability_Function_Pressure
 Function: MATCalcPressurePermeabilityMethod4
 Application to Urach

 Aufgabe:
   The normal stress across the fractures is calculated according to an approximate formulation
   from the relationship of stress with depth. This normal stress is then adjusted to take account of
   thermal cooling, and the resulting effective stress across the fracture isconverted into a permeabilty.
   by reference to a curve.

Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
E: long index: Elementnummer, double primary_variable pressure
R: relative permeability

Ergebnis:
rel. Permeabilitaet

Programmaenderungen:
09/2004   CMCD Inclusion in GeoSys 4
**************************************************************************/
double CMediumProperties::PermeabilityPressureFunctionMethod4(long index,double pressure, double temperature)
{
   //OK411
   index = index;
   pressure = pressure;
   temperature = temperature;
   double perm=0.0;

#ifdef obsolete                                //OK411
   /* Funktion der effektiven Spannung ueber Kurve CMCD 26.06.2003*/
   double x_mid, y_mid, z_mid, coords[3];
   double stress_eff,p,perm, thermal_stress;
   int idummy;
   /* Average depth */
   CalculateSimpleMiddelPointElement(index,coords);
   x_mid = coords[0];
   y_mid = coords[1];
   z_mid = coords[2];

   /*Calculate effective stress in the element*/
   p = pressure;
   stress_eff = (fabs(z_mid)*permeability_pressure_model_values[1]*1e6)-p;

   /*Impact of thermal stress*/
   thermal_stress = (permeability_pressure_model_values[2]-temperature)*permeability_pressure_model_values[3]*permeability_pressure_model_values[4];

   stress_eff = stress_eff - thermal_stress;
   if (stress_eff < 0.0)
   {
      stress_eff = 0.0;
   }
   /*Read value effective stress against curve*/
   perm = GetCurveValue((int) permeability_pressure_model_values[0], 0, stress_eff, &idummy);
#endif
   return perm;
}


double CMediumProperties::PermeabilityPorosityFunction(long number,double *gp,double theta)
{
   double k_rel=0.0;
   //OK411
   theta = theta;
   gp = gp;
   number = number;
#ifdef obsolete                                //OK411
   int i;
   long *element_nodes;
   double factora, factorb,valuelogk;

   //Collect primary variables
   static int nidx0,nidx1;
   double primary_variable[10];                   //OK To Do
   int count_nodes;
   int no_pcs_names =(int)pcs_name_vector.size();
   for(i=0;i<no_pcs_names;i++)
   {
      nidx0 = PCSGetNODValueIndex(pcs_name_vector[i],0);
      nidx1 = PCSGetNODValueIndex(pcs_name_vector[i],1);
      if(mode==0)                                 // Gauss point values
      {
         primary_variable[i] = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
            + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
      }
      else if(mode==1)                            // Node values
      {
         primary_variable[i] = (1.-theta)*GetNodeVal(number,nidx0) \
            + theta*GetNodeVal(number,nidx1);
      }
      else if(mode==2)                            // Element average value
      {
         count_nodes = ElNumberOfNodes[ElGetElementType(number) - 1];
         element_nodes = ElGetElementNodes(number);
         for (i = 0; i < count_nodes; i++)
            primary_variable[i] += GetNodeVal(element_nodes[i],nidx1);
         primary_variable[i]/= count_nodes;
      }
   }
   switch (permeability_porosity_model)
   {
      case 0:                                     //Reserved for a curve function
         k_rel=1.0;
         break;
      case 1:                                     //Constant value function
         k_rel=1.0;
         break;
      case 2:                                     //Ming Lian function
         factora=permeability_porosity_model_values[0];
         factorb=permeability_porosity_model_values[1];
         valuelogk=factora+(factorb*porosity);
         k_rel=0.987e-12*pow(10.,valuelogk);
         CMediumProperties *m_mmp = NULL;
         long group = ElGetElementGroupNumber(number);
         m_mmp = mmp_vector[group];
         double* k_ij;
         k_ij = m_mmp->PermeabilityTensor(number);//permeability;
         k_rel /= k_ij[0];
         break;
   }
#endif
   return k_rel;
}
