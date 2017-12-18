/**************************************************************************
FEMLib - Object: MFP Fluid Properties
Task:
Programing:
08/2004 OK Implementation
last modified:
**************************************************************************/
#include "makros.h"
// C++ STL
//#include <math.h>
//#include <fstream>
//#include <iostream>
#include <cfloat>

// FEM-Makros
//#include "mathlib.h"
#include "eos.h"                                  //NB
// GeoSys-GeoLib
#include "files0.h"
// GeoSys-FEMLib
#include "fem_ele_std.h"
//
#include "rf_mfp_new.h"
//#include "rf_mmp_new.h"
extern double InterpolValue(long number,int ndx,double r,double s,double t);
//#include "rf_pcs.h"
#include "rfmat_cp.h"
extern double GetCurveValue(int,int,double,int*);
#include "tools.h"                                //GetLineFromFile

using namespace std;

/* Umrechnungen SI - Amerikanisches System */
//WW #include "steam67.h"
#define PSI2PA 6895.
#define PA2PSI 1.4503263234227701232777374909355e-4
#define GAS_CONSTANT    8314.41
#define COMP_MOL_MASS_AIR    28.96
#define COMP_MOL_MASS_WATER  18.016
#define GAS_CONSTANT_V  461.5                     //WW
#define T_KILVIN_ZERO  273.15                     //AKS
double gravity_constant = 9.81;                   //TEST for FEBEX OK 9.81;

//==========================================================================
std::vector<CFluidProperties*>mfp_vector;

/**************************************************************************
FEMLib-Method:
Task: OBJ constructor
Programing:
08/2004 OK Implementation
**************************************************************************/
CFluidProperties::CFluidProperties(void)
{
   name = "WATER";
   phase = 0;
   // Density
   density_model = 1;
   rho_0 = 1000.;
   drho_dp = 0.;
   drho_dT = 0.;
   drho_dC = 0.;
   // Viscosity
   viscosity_model = 1;
   my_0 = 1e-3;
   dmy_dp = 0.;
   dmy_dT = 0.;
   dmy_dC = 0.;
   // Thermal properties
   heat_capacity_model = 1;
   specific_heat_capacity = 4680.;                //CMCD we should give this as SHC not HC GeoSys 4 9/2004
   heat_conductivity_model = 1;
   heat_conductivity = 0.6;
   // Electrical properties
   // Chemical properties
   diffusion_model = 1;
   diffusion = 2.13e-6;
   // State variables
   p_0 = 101325.;
   T_0 = 293.;
   C_0 = 0.;
   Z = 1.;
   cal_gravity = true;
   // Data
   mode = 0;                                      // Gauss point values
   Fem_Ele_Std = NULL;
   // WW
   molar_mass = COMP_MOL_MASS_AIR;
}


/**************************************************************************
FEMLib-Method:
Task: OBJ deconstructor
Programing:
05/2010 OK/AKS Implementation
**************************************************************************/
CFluidProperties::~CFluidProperties(void)
{

   for (int i=0; i < (int)component_vector.size() ; i++ )
   {
      component_vector[i] = NULL;
   }
   component_vector.clear();
}


/**************************************************************************
FEMLib-Method:
Task: OBJ read function
Programing:
08/2004 OK Implementation
11/2004 SB string streaming
**************************************************************************/
std::ios::pos_type CFluidProperties::Read(std::ifstream *mfp_file)
{
   std::string sub_line;
   std::string line_string;
   std::string delimiter(" ");
   bool new_keyword = false;
   std::string hash("#");
   std::ios::pos_type position;
   std::string sub_string;
   bool new_subkeyword = false;
   std::string dollar("$");
   std::string delimiter_type(":");
   std::stringstream in;
   //========================================================================
   // Schleife ueber alle Phasen bzw. Komponenten
   while (!new_keyword)
   {
      new_subkeyword = false;
      position = mfp_file->tellg();
      //SB    mfp_file->getline(buffer,MAX_ZEILE);
      //SB    line_string = buffer;
      line_string = GetLineFromFile1(mfp_file);
      if(line_string.size() < 1) break;
      if(line_string.find(hash)!=string::npos)
      {
         new_keyword = true;
         break;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$FLUID_TYPE")!=string::npos)
      {
         in.str(GetLineFromFile1(mfp_file));
         in >> name;                              //sub_line
         in.clear();
         continue;
      }
      //....................................................................
                                                  // NB 4.8.01
      if(line_string.find("$FLUID_NAME")!=string::npos)
      {
         in.str(GetLineFromFile1(mfp_file));
         in >> fluid_name;                        //sub_line
         therm_prop (fluid_name);                 //NB 4.9.05 (getting thermophysical constants of specified substance)
         //TODO: add choosing of property functions in input file (NB)
         in.clear();
         continue;
      }
      //....................................................................
                                                  // NB Oct-2009
      if(line_string.find("$COMPRESSIBILITY")!=string::npos)
      {
         in.str(GetLineFromFile1(mfp_file));
         in >> compressibility_model_pressure;    //sub_line 1 for first phase
         in >> compressibility_pressure;          //sub_line 1
         in.clear();
         in >> compressibility_model_temperature; //sub_line 2 for second phase
         in >> compressibility_temperature;       //sub_line 2
         in.clear();

         // available models see CFluidProperties::drhodP and CFluidProperties::drhodT
         // 0 incompressible
         // 1 constant slope
         // 2 slope from fct_table
         // 3 difference quotient
         // 4 analytical derivation

         in.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$DAT_TYPE")!=string::npos)
      {
         in.str(GetLineFromFile1(mfp_file));
         in >> name;                              //sub_line
         in.clear();
         continue;
      }
                                                  //YD/WW subkeyword found
      if(line_string.find("$NON_GRAVITY")!=string::npos)
      {
         cal_gravity = false;
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$DENSITY")!=string::npos)
      {
         new_subkeyword = false;
         in.str(GetLineFromFile1(mfp_file));
         in >> density_model;
         if(density_model==0)                     // rho = f(x)
         {
            in >> rho_fct_name;
         }
         if(density_model==1)                     // rho = const
         {
            in >> rho_0;
         }
         if(density_model==2)                     // rho(p) = rho_0*(1+beta_p*(p-p_0))
         {
            in >> rho_0;
            in >> p_0;
            in >> drho_dp;
            density_pcs_name_vector.push_back("PRESSURE1");
         }
         if(density_model==3)                     // rho(C) = rho_0*(1+beta_C*(C-C_0))
         {
            in >> rho_0;
            in >> C_0;
            in >> drho_dC;
            //        density_pcs_name_vector.push_back("CONCENTRATION1");
                                                  // PCH
            density_pcs_name_vector.push_back("Isochlor");
         }
         if(density_model==4)                     // rho(T) = rho_0*(1+beta_T*(T-T_0))
         {
            in >> rho_0;
            in >> T_0;
            in >> drho_dT;
            density_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(density_model==5)                     // rho(C,T) = rho_0*(1+beta_C*(C-C_0)+beta_T*(T-T_0))
         {
            in >> rho_0;
            in >> C_0;
            in >> drho_dC;
            in >> T_0;
            in >> drho_dT;
            density_pcs_name_vector.push_back("CONCENTRATION1");
            density_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(density_model==6)                     // rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0))
         {
            in >> rho_0;
            in >> p_0;
            in >> drho_dp;
            in >> T_0;
            in >> drho_dT;
            density_pcs_name_vector.push_back("PRESSURE1");
            density_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(density_model==7)                     // rho(p,p_v,T)
         {
            // no input data required
         }
         if(density_model==8)                     // rho(p,T,C)
         {
            in >> C_0;
            density_pcs_name_vector.push_back("PRESSURE1");
            density_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(density_model==9)                     //WW
         {                                        // Molar mass
            in >> molar_mass;
         }

         if((density_model==10)                   //NB 4.8.01  read density from a rho-P-T table
            ||(density_model==11)                 //NB 4.9.05  Peng-Robinson Equation of State
            ||(density_model==12)                 //NB 4.9.05  Redlich-Kwong Equation of State
            ||(density_model==13)                 //NB JUN 09  Fundamental equation
            ||(density_model==14))                //AKS MAY 10  Extended Ideal gas Eq. based on Super compressibility factor
         {
            std::string arg1,arg2,arg3;
            in >> arg1 >> arg2 >> arg3;           //get up to three arguments for density model

            if (isdigit(arg1[0])!=0)              // first argument is reference temperature
            {
               T_0 = atof(arg1.c_str());
               arg1 = arg2;
               arg2 = arg3;
            }

            if (arg1.length()==0)                 // if no arguments are given use standard
            {
               arg1="PRESSURE1";
               arg2="TEMPERATURE1";
            } else
            {
               if (arg2.length()==0)              // if only PRESSURE argument is given
               {
                  arg2="TEMPERATURE1";
               }
            }

            density_pcs_name_vector.push_back(arg1);
            if (T_Process) density_pcs_name_vector.push_back(arg2);
         }
         //      mfp_file->ignore(MAX_ZEILE,'\n');
         in.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$VISCOSITY")!=string::npos)
      {
         in.str(GetLineFromFile1(mfp_file));
         in >> viscosity_model;
         if(viscosity_model==0)                   // my = fct(x)
         {
            in >> my_fct_name;
         }
         if(viscosity_model==1)                   // my = const
         {
            in >> my_0;
         }
         if(viscosity_model==2)                   // my(p) = my_0*(1+gamma_p*(p-p_0))
         {
            in >> my_0;
            in >> p_0;
            in >> dmy_dp;
            viscosity_pcs_name_vector.push_back("PRESSURE1");
         }
         if(viscosity_model==3)                   // my(T), Yaws et al. (1976)
         {
                                                  //OK4704
            viscosity_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(viscosity_model==4)                   // my(T), ???
         {
         }
         if(viscosity_model==5)                   // my(p,T), Reichenberg (1971)
         {
         }
         if(viscosity_model==6)                   // my(C,T),
         {
         }
         if(viscosity_model==7)                   // my(p,T,C)
         {
            in >> C_0;
            viscosity_pcs_name_vector.push_back("PRESSURE1");
            viscosity_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(viscosity_model==9)                   // my(rho,T)
         {

            std::string arg1,arg2;
            in >> arg1 >> arg2;                   //get up to three arguments for density model

            if (arg1.length()==0)                 // if no arguments are given use standard
            {
               arg1="PRESSURE1";
               arg2="TEMPERATURE1";
            } else
            {
               if (arg2.length()==0)              // if only PRESSURE argument is given
               {
                  arg2="TEMPERATURE1";
               }
            }

            viscosity_pcs_name_vector.push_back(arg1);
            if (T_Process) viscosity_pcs_name_vector.push_back(arg2);

         }
         //    mfp_file->ignore(MAX_ZEILE,'\n');
         in.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$SPECIFIC_HEAT_CAPACITY")!=string::npos)
      {
         in.str(GetLineFromFile1(mfp_file));
         in >> heat_capacity_model;
         if(heat_capacity_model==0)               // c = fct(x)
         {
            in >> heat_capacity_fct_name;
         }
         if(heat_capacity_model==1)               // c = const
         {
            in >> specific_heat_capacity;
            specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
            specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(heat_capacity_model==2)               // my(p,T,C)
         {
            in >> C_0;
            specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
            specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(heat_capacity_model==3)               // YD: improved phase change
         {
            in >> T_Latent1;                      // Tmin for phase change
            in >> T_Latent2;                      // Tmax for phase change
            in >> heat_phase_change_curve;
            specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
            specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
            specific_heat_capacity_pcs_name_vector.push_back("SATURATION1");
            enthalpy_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(heat_capacity_model==4)               // YD: improved phase change, function
         {
            in >> T_Latent1;                      // Tmin for phase change
            in >> T_Latent2;                      // Tmax for phase change
            in >> specific_heat_capacity;         // ^c
            in >> latent_heat;                    // L
            specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
            specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
            specific_heat_capacity_pcs_name_vector.push_back("SATURATION1");
            enthalpy_pcs_name_vector.push_back("TEMPERATURE1");
         }
         in.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$HEAT_CONDUCTIVITY")!=string::npos)
      {
         in.str(GetLineFromFile1(mfp_file));
         in >> heat_conductivity_model;
         if(heat_conductivity_model==0)           // my = fct(x)
         {
            in >> heat_conductivity_fct_name;
         }
         if(heat_conductivity_model==1)           // my = const
         {
            in >> heat_conductivity;
         }
         if(heat_conductivity_model==2)           // my = f(p,T,C)
         {
            in >> C_0;
            heat_conductivity_pcs_name_vector.push_back("PRESSURE1");
            heat_conductivity_pcs_name_vector.push_back("TEMPERATURE1");
         }
         if(heat_conductivity_model==3)           // my = f(p,T) NB
         {

            std::string arg1,arg2;
            in >> arg1 >> arg2;                   //get up to three arguments for density model

            if (arg1.length()==0)                 // if no arguments are given use standard
            {
               arg1="PRESSURE1";
               arg2="TEMPERATURE1";
            } else
            {
               if (arg2.length()==0)              // if only PRESSURE argument is given
               {
                  arg2="TEMPERATURE1";
               }
            }

            heat_conductivity_pcs_name_vector.push_back(arg1);
            heat_conductivity_pcs_name_vector.push_back(arg2);
         }
         in.clear();                              //OK
         continue;
      }
                                                  // subkeyword found
      if(line_string.find("$PHASE_DIFFUSION")!=string::npos)
      {
         in.str(GetLineFromFile1(mfp_file));
         in >> diffusion_model;
         if(diffusion_model==0)                   // D = fct(x)
         {
            in >> dif_fct_name;
         }
         if(diffusion_model==1)                   // D = const //MX
         {
            in >> diffusion;
         }
         in.clear();
         continue;
      }
      //....................................................................
                                                  // CMCD outer space version
      if(line_string.find("$GRAVITY")!=string::npos)
      {
         in.str(GetLineFromFile1(mfp_file));
         in >> gravity_constant;
         in.clear();
         continue;
      }
      //....................................................................
   }
   return position;
}


/**************************************************************************
FEMLib-Method:
Task: Master read function
Programing:
08/2004 OK Implementation
01/2005 OK Boolean type
01/2005 OK Destruct before read
**************************************************************************/
bool MFPRead(std::string file_base_name)
{
   //----------------------------------------------------------------------
   //OK  MFPDelete();
   //----------------------------------------------------------------------
   CFluidProperties *m_mfp = NULL;
   char line[MAX_ZEILE];
   std::string sub_line;
   std::string line_string;
   std::ios::pos_type position;
   //========================================================================
   // File handling
   std::string mfp_file_name = file_base_name + MFP_FILE_EXTENSION;
   std::ifstream mfp_file (mfp_file_name.data(),std::ios::in);
   if (!mfp_file.good())
      return false;
   mfp_file.seekg(0L,std::ios::beg);
   //========================================================================
   // Keyword loop
   std::cout << "MFPRead" << std::endl;
   while (!mfp_file.eof())
   {
      mfp_file.getline(line,MAX_ZEILE);
      line_string = line;
      if(line_string.find("#STOP")!=std::string::npos)
         return true;
      //----------------------------------------------------------------------
                                                  // keyword found
      if(line_string.find("#FLUID_PROPERTIES")!=std::string::npos)
      {
         m_mfp = new CFluidProperties();
         m_mfp->file_base_name = file_base_name;
         position = m_mfp->Read(&mfp_file);
         m_mfp->phase = (int)mfp_vector.size();   //OK4108
         mfp_vector.push_back(m_mfp);
         mfp_file.seekg(position,std::ios::beg);
      }                                           // keyword found
   }                                              // eof
   //========================================================================
   // Configuration
   int i;
   int no_fluids =(int)mfp_vector.size();
   if(no_fluids==1)
   {
      m_mfp = mfp_vector[0];
      m_mfp->phase=0;
   }
   else if(no_fluids==2)
   {
      for(i=0;i<no_fluids;i++)
      {
         m_mfp = mfp_vector[i];
         if(m_mfp->name.find("GAS")!=string::npos)
            m_mfp->phase=0;
         else
            m_mfp->phase=1;
      }
   }
   //----------------------------------------------------------------------
   // Test
   if(mfp_vector.size()==0)
   {
      std::cout << "Error in MFPRead: no MFP data" << std::endl;
      abort();
   }
   //----------------------------------------------------------------------
   return true;
}


/**************************************************************************
FEMLib-Method:
Task: write function
Programing:
11/2004 SB Implementation
last modification:
**************************************************************************/
void CFluidProperties::Write(std::ofstream* mfp_file)
{
   //KEYWORD
   *mfp_file  << "#FLUID_PROPERTIES" << std::endl;
   *mfp_file  << " $FLUID_TYPE" << std::endl;
   *mfp_file  << "  " << name << std::endl;
   *mfp_file  << " $DAT_TYPE" << std::endl;
   *mfp_file  << "  " << name << std::endl;
   *mfp_file  << " $DENSITY" << std::endl;
   if(density_model == 0) *mfp_file << "  " << density_model << " " << rho_fct_name << std::endl;
   if(density_model == 1) *mfp_file << "  " << density_model << " " << rho_0 << std::endl;
   //todo
   *mfp_file  << " $VISCOSITY" << endl;
   if(viscosity_model == 0) *mfp_file << "  " << viscosity_model << " " << my_fct_name << std::endl;
   if(viscosity_model == 1) *mfp_file << "  " << viscosity_model << " " << my_0 << std::endl;
   //todo
   *mfp_file  << " $SPECIFIC_HEAT_CAPACITY" << endl;
   if(heat_capacity_model == 0) *mfp_file << "  " << heat_capacity_model << " " << heat_capacity_fct_name << std::endl;
   if(heat_capacity_model == 1) *mfp_file << "  " << heat_capacity_model << " " << specific_heat_capacity << std::endl;
   *mfp_file  << " $SPECIFIC_HEAT_CONDUCTIVITY" << endl;
   if(heat_conductivity_model == 0) *mfp_file << "  " << heat_conductivity_model << " " << heat_conductivity_fct_name << std::endl;
   if(heat_conductivity_model == 1) *mfp_file << "  " << heat_conductivity_model << " " << heat_conductivity << std::endl;
   //--------------------------------------------------------------------
}


/**************************************************************************
FEMLib-Method:
Task: Master write function
Programing:
08/2004 OK Implementation
last modification:
**************************************************************************/
void MFPWrite(std::string base_file_name)
{
   CFluidProperties *m_mfp = NULL;
   string sub_line;
   string line_string;
   ofstream mfp_file;
   //========================================================================
   // File handling
   std::string mfp_file_name = base_file_name + MFP_FILE_EXTENSION;
   mfp_file.open(mfp_file_name.data(),std::ios::trunc|std::ios::out);
   mfp_file.setf(std::ios::scientific,std::ios::floatfield);
   mfp_file.precision(12);
   if (!mfp_file.good()) return;
   //  mfp_file.seekg(0L,ios::beg);
   //========================================================================
   mfp_file << "GeoSys-MFP: Material Fluid Properties -------------" << std::endl;
   //========================================================================
   // OUT vector
   int mfp_vector_size =(int)mfp_vector.size();
   int i;
   for(i=0;i<mfp_vector_size;i++)
   {
      m_mfp = mfp_vector[i];
      m_mfp->Write(&mfp_file);
   }
   mfp_file << "#STOP";
   mfp_file.close();
   //  delete mfp_file;
}


////////////////////////////////////////////////////////////////////////////
// Properties functions
////////////////////////////////////////////////////////////////////////////
/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
09/2005 WW implementation
11/2005 YD modification
11/2005 CMCD Inclusion current and previous time step quantities
05/2007 PCH improvement for density-dependent flow
last modification:
**************************************************************************/
void CFluidProperties::CalPrimaryVariable(std::vector<std::string>& pcs_name_vector)
{
   CRFProcess* m_pcs = NULL;

   int nidx0,nidx1;
   if(!Fem_Ele_Std)                               //OK
      return;

   primary_variable[0] = 0;
   primary_variable[1] = 0;
   primary_variable[2] = 0;

   for(int i=0;i<(int)pcs_name_vector.size();i++)
   {
      //MX  m_pcs = PCSGet("HEAT_TRANSPORT");
      m_pcs = PCSGet(pcs_name_vector[i],true);
      if (!m_pcs) return;                         //MX
      nidx0 = m_pcs->GetNodeValueIndex(pcs_name_vector[i]);
      nidx1 = nidx0+1;

      if(mode==0)                                 // Gauss point values
      {
         primary_variable_t0[i]= Fem_Ele_Std->interpolate(nidx0,m_pcs);
         primary_variable_t1[i]= Fem_Ele_Std->interpolate(nidx1,m_pcs);
         primary_variable[i] = (1.-Fem_Ele_Std->pcs->m_num->ls_theta)*Fem_Ele_Std->interpolate(nidx0,m_pcs)+
            Fem_Ele_Std->pcs->m_num->ls_theta*Fem_Ele_Std->interpolate(nidx1,m_pcs);

      }
      else if(mode==2)                            // Element average value
      {
         primary_variable[i] = (1.-Fem_Ele_Std->pcs->m_num->ls_theta)*Fem_Ele_Std->elemnt_average(nidx0,m_pcs)
            + Fem_Ele_Std->pcs->m_num->ls_theta*Fem_Ele_Std->elemnt_average(nidx1,m_pcs);
         primary_variable_t0[i]= Fem_Ele_Std->elemnt_average(nidx0,m_pcs);
         primary_variable_t1[i]= Fem_Ele_Std->elemnt_average(nidx1,m_pcs);
      }
      if (mode==3)                                //NB, just testing
      {
         primary_variable[i]=Fem_Ele_Std->interpolate(nidx0,m_pcs);

      }
   }
}


////////////////////////////////////////////////////////////////////////////
// Fluid density
/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
08/2004 OK MFP implementation
           based on MATCalcFluidDensity by OK/JdJ,AH,MB
11/2005 YD Modification
05/2008 WW Add an argument: double* variables: P, T, C
last modification:
NB 4.9.05
**************************************************************************/
double CFluidProperties::Density(double* variables)
{
   static double density;
   // static double air_gas_density,vapour_density,vapour_pressure;
   int fct_number = 0;
   int gueltig;
   //----------------------------------------------------------------------
   if(variables)                                  // This condition is added by WW
   {
      //----------------------------------------------------------------------
      // Duplicate the following lines just to enhance computation. WW
      switch(density_model)
      {
         case 0:                                  // rho = f(x)
            density = GetCurveValue(fct_number,0,variables[0],&gueltig);
            break;
         case 1:                                  // rho = const
            density = rho_0;
            break;
         case 2:                                  // rho(p) = rho_0*(1+beta_p*(p-p_0))
            density = rho_0*(1.+drho_dp*(max(variables[0],0.0)-p_0));
            break;
         case 3:                                  // rho(C) = rho_0*(1+beta_C*(C-C_0))
            density = rho_0*(1.+drho_dC*(max(variables[2],0.0)-C_0));
            break;
         case 4:                                  // rho(T) = rho_0*(1+beta_T*(T-T_0))
            density = rho_0*(1.+drho_dT*(max(variables[1],0.0)-T_0));
            break;
         case 5:                                  // rho(C,T) = rho_0*(1+beta_C*(C-C_0)+beta_T*(T-T_0))
            density = rho_0*(1.+drho_dC*(max(variables[2],0.0)-C_0)+drho_dT*(max(variables[1],0.0)-T_0));
            break;
         case 6:                                  // rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0))
            density = rho_0*(1.+drho_dp*(max(variables[0],0.0)-p_0)+drho_dT*(max(variables[1],0.0)-T_0));
            break;
         case 7:                                  // Pefect gas. WW
            density = variables[0]*molar_mass/(GAS_CONSTANT*variables[1]);
            break;
         case 10:                                 // Get density from temperature-pressure values from fct-file	NB 4.8.01
            if(!T_Process) variables[1]=T_0;
            density = GetMatrixValue(variables[1],variables[0],fluid_name,&gueltig);
            break;
         case 11:                                 //Peng-Robinson EOS for different fluids NB 4.9.05
            if(!T_Process) variables[1]=T_0;
            density = rkeos(variables[1],variables[0],fluid_id);

            break;
         case 12:                                 // Redlich-Kwong EOS for different fluids NB 4.9.05
            if(!T_Process) variables[1]=T_0;
                                                  //NB
            density = preos(variables[1],variables[0],fluid_id);
            break;
         case 13:                                 // Helmholtz free Energy NB JUN 09
                                                  //NB
            density = zbrent(variables[1],variables[0],fluid_id,1e-8);
            break;
         case 14:                                 //AKS empiricaly extented Ideal gas Eq for real gas // it has used with fractional mass transport Eq.//
            density = MixtureSubProperity(5, (long)  variables[2], variables[0], variables[1])* variables[0] / (CalCopressibility((long)  variables[2], variables[0], variables[1] )*variables[1] * GAS_CONSTANT) ;
            break;
         default:
            std::cout << "Error in CFluidProperties::Density: no valid model" << std::endl;
            break;
      }
   }
   else
   {
      CalPrimaryVariable(density_pcs_name_vector);
      //----------------------------------------------------------------------
      switch(density_model)
      {
         case 0:                                  // rho = f(x)
            density = GetCurveValue(fct_number,0,primary_variable[0],&gueltig);
            break;
         case 1:                                  // rho = const
            density = rho_0;
            break;
         case 2:                                  // rho(p) = rho_0*(1+beta_p*(p-p_0))
            density = rho_0*(1.+drho_dp*(max(primary_variable[0],0.0)-p_0));
            break;
         case 3:                                  // rho(C) = rho_0*(1+beta_C*(C-C_0))
            density = rho_0*(1.+drho_dC*(max(primary_variable[0],0.0)-C_0));
            break;
         case 4:                                  // rho(T) = rho_0*(1+beta_T*(T-T_0))
            density = rho_0*(1.+drho_dT*(max(primary_variable[0],0.0)-T_0));
            break;
         case 5:                                  // rho(C,T) = rho_0*(1+beta_C*(C-C_0)+beta_T*(T-T_0))
            density = rho_0*(1.+drho_dC*(max(primary_variable[0],0.0)-C_0)+drho_dT*(max(primary_variable[1],0.0)-T_0));
            break;
         case 6:                                  // rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0))
            density = rho_0*(1.+drho_dp*(max(primary_variable[0],0.0)-p_0)+drho_dT*(max(primary_variable[1],0.0)-T_0));
            break;
         case 7:                                  // rho_w^l(p,T) for gas phase
            /* //WW
            vapour_pressure = MFPCalcVapourPressure(primary_variable[0]);
            air_gas_density = (COMP_MOL_MASS_AIR * (primary_variable[1]-vapour_pressure)) / (GAS_CONSTANT*(primary_variable[0]+0.0));
            vapour_density = (COMP_MOL_MASS_WATER*vapour_pressure) / (GAS_CONSTANT*(primary_variable[0]+0.0));
            density = vapour_density + air_gas_density;
             */
            break;
         case 8:                                  // M14 von JdJ
            density = MATCalcFluidDensityMethod8(primary_variable[0],primary_variable[1],primary_variable[2]);
            break;
         case 10:                                 // Get density from temperature-pressure values from fct-file NB
            if(!T_Process) primary_variable[1]=T_0;
            density = GetMatrixValue(primary_variable[1],primary_variable[0],fluid_name,&gueltig);
            break;
         case 11:                                 //Peng-Robinson equation of state NB
            if(!T_Process) primary_variable[1]=T_0;
            density = rkeos(primary_variable[1],primary_variable[0],fluid_id);
            break;
         case 12:                                 //Redlich-Kwong equation of state NB
            if(!T_Process) primary_variable[1]=T_0;
            density = preos(primary_variable[1],primary_variable[0],fluid_id);
            break;
         case 13:                                 // Helmholtz free Energy NB JUN 09
            if(!T_Process) primary_variable[1]=T_0;
                                                  //NB
            density = zbrent(primary_variable[1],primary_variable[0],fluid_id,1e-8);
            break;
         case 14:                                 //AKS empiricaly extented Ideal gas Eq for real gas
            density = MixtureSubProperity(5, (long)  variables[2], variables[0], variables[1])* variables[0] / (CalCopressibility((long)  variables[2], variables[0], variables[1] )*variables[1] * GAS_CONSTANT) ;
            //Replace MixtureSubProperity(5, (long)  variables[2], variables[0], variables[1])* variables[0] by molar_mass in case of no MASS_TRANSPORT proces
         default:
            std::cout << "Error in CFluidProperties::Density: no valid model" << std::endl;
            break;
      }
   }
   return density;
}


/*************************************************************************
 ROCKFLOW - Funktion: MATCalcFluidDensityMethod8

 Task:
   Density of a geothermal fluid as a function of temperature, pressure
   and concentration.

   Fluid-Density function according to IAPWS-IF97

 Programmaenderungen:
   09/2003   CMCD  ARL  First implementation
09/2004   CMCD  Inclusion in GeoSys vs. 4

*************************************************************************/
double CFluidProperties::MATCalcFluidDensityMethod8(double Press, double TempK, double Conc)
{
   Conc = Conc;
   /*int c_idx;*/
   double rho_0;
   double GammaPi, Pressurevar, Tau, pressure_average, temperature_average;
   double Tstar, Pstar,GazConst;
   double L[35],J[35],n[35];
   int i;
   double salinity;

   pressure_average = Press;
   temperature_average = TempK;
   salinity = C_0;
   Tstar = 1386;
   Pstar = 16.53e6;                               // MPa
   GazConst = 0.461526e3;                         //

   n[0] = 0.0;
   n[1] = 0.14632971213167;
   n[2] = -0.84548187169114;
   n[3] = -0.37563603672040e1;
   n[4] = 0.33855169168385e1;
   n[5] = -0.95791963387872;
   n[6] = 0.15772038513228;
   n[7] = -0.16616417199501e-1;
   n[8] = 0.81214629983568e-3;
   n[9] = 0.28319080123804e-3;
   n[10] = -0.60706301565874e-3;
   n[11] = -0.18990068218419e-1;
   n[12] = -0.32529748770505e-1;
   n[13] = -0.21841717175414e-1;
   n[14] = -0.52838357969930e-4;
   n[15] = -0.47184321073267e-3;
   n[16] = -0.30001780793026e-3;
   n[17] = 0.47661393906987e-4;
   n[18] = -0.44141845330846e-5;
   n[19] = -0.72694996297594e-15;
   n[20] = -0.31679644845054e-4;
   n[21] = -0.28270797985312e-5;
   n[22] = -0.85205128120103e-9;
   n[23] = -0.22425281908000e-5;
   n[24] = -0.65171222895601e-6;
   n[25] = -0.14341729937924e-12;
   n[26] = -0.40516996860117e-6;
   n[27] = -0.12734301741641e-8;
   n[28] = -0.17424871230634e-9;
   n[29] = -0.68762131295531e-18;
   n[30] = 0.14478307828521e-19;
   n[31] = 0.26335781662795e-22;
   n[32] = -0.11947622640071e-22;
   n[33] = 0.18228094581404e-23;
   n[34] = -0.93537087292458e-25;

   L[0] =0.;
   L[1] =0.;
   L[2] =0.;
   L[3] =0.;
   L[4] =0.;
   L[5] =0.;
   L[6] =0.;
   L[7] =0.;
   L[8] =0.;
   L[9] =1.;
   L[10] =1.;
   L[11] =1.;
   L[12] =1.;
   L[13] =1.;
   L[14] =1.;
   L[15] =2.;
   L[16] =2.;
   L[17] =2.;
   L[18] =2.;
   L[19] =2.;
   L[20] =3.;
   L[21] =3.;
   L[22] =3.;
   L[23] =4.;
   L[24] =4.;
   L[25] =4.;
   L[26] =5.;
   L[27] =8.;
   L[28] =8.;
   L[29] =21.;
   L[30] =23.;
   L[31] =29.;
   L[32] =30.;
   L[33] =31.;
   L[34] =32.;

   J[0] =-2.;
   J[1] =-2.;
   J[2] =-1.;
   J[3] =0.;
   J[4] =1.;
   J[5] =2.;
   J[6] =3.;
   J[7] =4.;
   J[8] =5.;
   J[9] =-9.;
   J[10] =-7.;
   J[11] =-1.;
   J[12] =0.;
   J[13] =1.;
   J[14] =3.;
   J[15] =-3.;
   J[16] =0.;
   J[17] =1.;
   J[18] =3.;
   J[19] =17.;
   J[20] =-4.;
   J[21] =0.;
   J[22] =6.;
   J[23] =-5.;
   J[24] =-2.;
   J[25] =10.;
   J[26] =-8.;
   J[27] =-11.;
   J[28] =-6.;
   J[29] =-29.;
   J[30] =-31.;
   J[31] =-38.;
   J[32] = -39.;
   J[33] =-40.;
   J[34] =-41.;

   Pressurevar = pressure_average / Pstar;
   Tau = Tstar / temperature_average;

   /*BEGIN:Calculation of GammaPi*/
   GammaPi = 0.;

   for (i=1; i<35; i++)
   {

      GammaPi = GammaPi - (n[i]) * (L[i]) * (pow((7.1-Pressurevar),(L[i] -1))) * (pow((Tau-1.222),J[i]));

   }
   /*END: Calculation of GammaPi*/

   /*BEGIN: Calculation of density*/
   rho_0 = pressure_average / (GazConst * temperature_average * Pressurevar * GammaPi);
   /*END: Calculation of density*/

   /*  return rho_0 + drho_dC * (concentration_average - c0);   */
   /*printf("%f", rho_0 + salinity);*/
   return rho_0 + salinity;
}


////////////////////////////////////////////////////////////////////////////
// Fluid viscosity
/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
08/2004 OK Implementation
11/2005 YD Modification
10/2008 OK Faster data access
03/2009 NB Viscosity depending on Density()
**************************************************************************/
                                                  //OK4709
double CFluidProperties::Viscosity(double* variables)
{
   static double viscosity;
   int fct_number = 0;
   int gueltig;
   double mfp_arguments[2];
   double density;

   // double TTT=0, PPP=0;
   // long Element_Index;
   // int nod_index;

   // string val_name;
   CRFProcess* m_pcs = NULL;

   //----------------------------------------------------------------------
   bool New = false;                              // To be
   if(fem_msh_vector.size()>0) New = true;

   //----------------------------------------------------------------------
   if(variables)                                  //OK4709: faster data access
   {
      primary_variable[0] = variables[0];         //p (single phase)
      primary_variable[1] = variables[1];         //T (temperature)
      primary_variable[2] = variables[2];         //C (salinity)//index in case of AIR_FLOW

   }
   else
   {
      CalPrimaryVariable(viscosity_pcs_name_vector);
   }
   //----------------------------------------------------------------------
   switch(viscosity_model)
   {
      case 0:                                     // rho = f(x)
         viscosity = GetCurveValue(fct_number,0,primary_variable[0],&gueltig);
         break;
      case 1:                                     // my = const
         viscosity = my_0;
         break;
      case 2:                                     // my(p) = my_0*(1+gamma_p*(p-p_0))
         viscosity = my_0*(1.+dmy_dp*(max(primary_variable[0],0.0)-p_0));
         break;
      case 3:                                     // my^l(T), Yaws et al. (1976)
         if(mode==1)                              //OK4704 for nodal output
         {
            m_pcs = PCSGet("HEAT_TRANSPORT");
            //if(!m_pcs) return 0.0;
            primary_variable[1] = m_pcs->GetNodeValue(node,m_pcs->GetNodeValueIndex("TEMPERATURE1")+1);
         }
                                                  //ToDo pcs_name
         viscosity = LiquidViscosity_Yaws_1976(primary_variable[1]);
         break;
      case 4:                                     // my^g(T), Marsily (1986)
         viscosity = LiquidViscosity_Marsily_1986(primary_variable[0]);
         break;
      case 5:                                     // my^g(p,T), Reichenberg (1971)
         viscosity = GasViscosity_Reichenberg_1971(primary_variable[0],primary_variable[1]);
         break;
      case 6:                                     // my(C,T),
         viscosity = LiquidViscosity_NN(primary_variable[0],primary_variable[1]);
         break;
      case 7:                                     // my(p,C,T),
         viscosity = LiquidViscosity_CMCD(primary_variable[0],primary_variable[1],primary_variable[2]);
         break;
      case 9:                                     // viscosity as function of density and temperature, NB

         if(!T_Process) primary_variable[1]=T_0;
         // TODO: switch case different models...
         // Problem.cpp 3 PS_GLOBAL, 1212,1313 pcs_type
         //TODO: default fluid_ID, if not specified
         mfp_arguments[0] = primary_variable[0];  // rescue primary_variable before its destroyed by Density();
         mfp_arguments[1] = primary_variable[1];

         density = Density(mfp_arguments);        //TODO: (NB) store density (and viscosity) as secondary variable
                                                  //NB
         viscosity = Fluid_Viscosity(density,mfp_arguments[1],mfp_arguments[0],fluid_id);

         break;

      case 10:                                    // my(rho, T) for real gases mixture
         viscosity = GasViscosity_Chung_1988((long)  primary_variable[2], primary_variable[0], primary_variable[1]);
         break;
      default:
         cout << "Error in CFluidProperties::Viscosity: no valid model" << endl;
         break;
   }
   //----------------------------------------------------------------------

   return viscosity;

}


/**************************************************************************
FEMLib-Method:
Task:
   Dynamische Gasviskositaet nach Reichenberg (1971)
   als Funktion von Temperatur und Druck
   in Reid et al. (1988), S. 420
Programing:
08/2004 OK MFP implementation based on CalcFluidViscosityMethod7 by OK
last modification:
**************************************************************************/
double CFluidProperties::GasViscosity_Reichenberg_1971(double p,double T)
{
   double my,my0;
   double A,B,C,D;
   double Q,Pr,Tr;
   double alpha1, alpha2, beta1, beta2, gamma1, gamma2, delta1, delta2;
   double a,c,d;
   double pc,Tc;

   alpha1 = 1.9824e-3;
   alpha2 = 5.2683;
   beta1  = 1.6552;
   beta2  = 1.2760;
   gamma1 = 0.1319;
   gamma2 = 3.7035;
   delta1 = 2.9496;
   delta2 = 2.9190;
   a = -0.5767;
   c = -79.8678;
   d = -16.6169;
   Q = 1.0;
   Tc = 126.2;
   Tr = T/Tc;
   pc = 33.9*10000.;                              /* bar->Pascal*/
   Pr = p/pc;

   A = alpha1/Tr * exp(alpha2*pow(Tr,a));
   B = A*(beta1*Tr-beta2);
   C = gamma1/Tr * exp(gamma2*pow(Tr,c));
   D = delta1/Tr * exp(delta2*pow(Tr,d));

                                                  /* Poise->Pa*s */
   my0 = 26.69*sqrt(28.96)*sqrt(T)/(3.7*3.7) * 1.e-6 * 0.1;
   my = my0 * ( 1.0 + Q * (A * pow(Pr,1.5))/(B*Pr+(1/(1+C*pow(Pr,D)))) );
   return my;
}


/**************************************************************************
FEMLib-Method:
Task:
  Dynamic viscosity of real gaseous mixture
   als Funktion von Temperatur und Density
Programing:
05/2010  Implementation  AKS

**************************************************************************/
double CFluidProperties::GasViscosity_Chung_1988(long idx_elem, double p,double T)
{
   double my,my0,myp,myk,Omega;
   double A,B,C,D,E,F,G,H,S,W,M,Vc,T_str,Fc,Tc,Y,G1,G2,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10;
   double dens_arg[3];
   A = 1.16145;
   B = 0.14874;
   C = 0.52487;
   D = 0.77320;
   E = 2.16178;
   F = 2.43787;
   G = -6.435e-4;
   H = 7.27371;
   S = 18.0323;
   W = -0.7683;
   A1 = 6.32403 + 50.41190*MixtureSubProperity(4, idx_elem, p, T);
   A2 = 1.2102E-03 - 1.1536E-03*MixtureSubProperity(4, idx_elem, p, T);
   A3 = 5.28346 + 254.209*MixtureSubProperity(4, idx_elem, p, T);
   A4 = 6.62263 + 38.0957*MixtureSubProperity(4, idx_elem, p, T);
   A5 = 19.74540 + 7.63034*MixtureSubProperity(4, idx_elem, p, T);
   A6 = -1.89992 - 12.53670*MixtureSubProperity(4, idx_elem, p, T);
   A7 = 24.27450 + 3.44945*MixtureSubProperity(4, idx_elem, p, T);
   A8 = 0.79716 + 1.11764*MixtureSubProperity(4, idx_elem, p, T);
   A9 = -0.23816 + 0.067695*MixtureSubProperity(4, idx_elem, p, T);
   A10 = 0.068629 + 0.34793*MixtureSubProperity(4, idx_elem, p, T);

   dens_arg[0] = p;
   dens_arg[1] = T;
   dens_arg[2] = idx_elem;

                                                  // 1.2593*mixture_energy_parameter
   Tc = 1.2593*MixtureSubProperity(3, idx_elem, p, T);
                                                  // 1 - 0.2756*mixture_acentric_factor
   Fc = 1.0 - 0.2756*MixtureSubProperity(4, idx_elem, p, T);
   M = MixtureSubProperity(5, idx_elem, p, T);    // mixture_molecular_weight
   T_str = T/MixtureSubProperity(3, idx_elem, p, T);
                                                  // mixture_critical_volume UNIT: m3/kmole
   Vc = pow(MixtureSubProperity(2, idx_elem, p, T)/0.809, 3.0);
   Y = Density(dens_arg)*Vc/(M*6000.0);
   G1 = (1 - 0.5*Y)*pow((1-Y), -3.0);
   G2 = ((A1*(1 - exp(-A4*Y))/Y) + A2*G1*exp(A5*Y)+ A3*G1)/(A1*A4 + A2 + A3);
   Omega = A*pow(T_str, -B) + C*exp(-T_str*D) + E*exp(-T_str*F) + G*pow(T_str, B)*sin(S*pow(T_str, W) - H);

   my0 = (4.0785e-05)*pow(M*T, 0.5)*Fc/(pow(Vc,0.6666)*Omega);
   myk = my0*(pow(G2, -1.0) + A6*Y);
   myp = ( (3.6344e-05)*pow(M*Tc, 0.5)/pow(Vc, 0.6666) )*A7*Y*Y*G2*exp(A8 + A9*pow(T_str, -1.0) + A10*pow(T_str, -2.0));
   my = 0.1*(myp + myk);                          // g/(cm.s)=0.1 Pa.s
   return my;
}


/**************************************************************************
FEMLib-Method:
Task:
   Dynamische Fluessigkeits-Viskositaet nach Yaws et al. (1976)
   als Funktion von Temperatur
   in Reid et al. (1988), S. 441/455
   Eqn.(3): ln(my) = A + B/T + CT + DT^2
Programing:
08/2004 OK MFP implementation
           based on CalcFluidViscosityMethod8 by OK (06/2001)
last modification:
**************************************************************************/
double CFluidProperties::LiquidViscosity_Yaws_1976(double T)
{
   double ln_my,my;
   double A,B,C,D;

   A = -2.471e+01;
   B =  4.209e+03;
   C =  4.527e-02;
   D = -3.376e-5;

   ln_my = A + B/T + C*T + D*T*T;
   my = exp(ln_my);                               /* in cP */
   my = my * 1.e-3;                               /* in Pa s */
   return my;
}


/**************************************************************************
FEMLib-Method:
Task:
   Fluessigkeits-Viskositaet in Abhaengigkeit von der Temperatur
   (nach Marsily 1986)
Programing:
08/2004 OK MFP implementation
           based on CalcFluidViscosityMethod9 by OK (05/2001)
last modification:
**************************************************************************/
double CFluidProperties::LiquidViscosity_Marsily_1986(double T)
{
   double my;
   my = 2.285e-5 + 1.01e-3*log(T);
   return my;
}


/**************************************************************************
FEMLib-Method:
Task:
   Liefert die Viskositaet in Abhaengigkeit von der Konzentration
   und der Temperatur.
Programing:
02/2000 AH Erste Version
08/2004 OK MFP implementation
last modification:
**************************************************************************/
double CFluidProperties::LiquidViscosity_NN(double c,double T)
{
   double f1, f2, mu0 = 0.001, mu;
   double omega0, omega, sigma0, sigma;

   if (rho_0 <  MKleinsteZahl|| T_0 < MKleinsteZahl)
      return 0.;

   omega = c / rho_0;
   omega0 = C_0 / rho_0;
   sigma = (T - T_0) / T_0;
   sigma0 = 0.;

   f1 = (1. + 1.85 * omega0 - 4.1 * omega0 * omega0 + 44.5 * omega0 * omega0 * omega0) /
      (1. + 1.85 * omega - 4.1 * omega * omega + 44.5 * omega * omega * omega);
   f2 = (1 + 0.7063 * sigma - 0.04832 * sigma * sigma * sigma) /
      (1 + 0.7063 * sigma0 - 0.04832 * sigma0 * sigma * sigma0);
   mu = mu0 / (f1 + f2);
   return mu;
}


////////////////////////////////////////////////////////////////////////////
// Fluid thermal properties
/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
08/2004 OK MFP implementation based on MATCalcFluidHeatCapacity (OK)
10/2005 WW/YD Case 3: simplified phase change
10/2005 YD Case 4: improved phase change
11/2005 CMCD edited cases and expanded case 3 & 4
last modification: NB Jan 09 4.9.05
**************************************************************************/
                                                  //NB
double CFluidProperties::SpecificHeatCapacity(double *variables)
{
   int gueltig = -1;
   double pressure, saturation, temperature;

   if(variables)                                  //NB Jan 09
   {
      primary_variable[0] = variables[0];         //p (single phase)
      primary_variable[1] = variables[1];         //T (temperature)
      primary_variable[2] = variables[2];         //ELE index
   }
   else
   {
      CalPrimaryVariable(specific_heat_capacity_pcs_name_vector);
   }

   pressure = primary_variable[0];
   temperature =  primary_variable[1];
   saturation = primary_variable[2];
   //......................................................................
   //
   switch(heat_capacity_model)
   {
      case 0:                                     // c = f(x)
         specific_heat_capacity = GetCurveValue(0,0,temperature,&gueltig);
         break;
      case 1:                                     // c = const, value already read in to specific_heat_capacity
         break;
      case 2:                                     // c = f(p,T,Conc)
         specific_heat_capacity = MATCalcFluidHeatCapacityMethod2(pressure, temperature, saturation);
         break;
      case 5:
         specific_heat_capacity = GetCurveValue(heat_phase_change_curve,0,temperature_buffer,&gueltig);
         break;
      case 10:                                    //AKS for real gases mixture
                                                  // mixture cp= sum_i sum_j x_i*x_j*intrc* sqrt[cp_i(p,T)*cp_j(p,T)], with interation
         specific_heat_capacity = MixtureSubProperity(7, (long)  primary_variable[2], primary_variable[0], primary_variable[1]);
         break;
   }
   return specific_heat_capacity;
}


/**************************************************************************
FEMLib-Method:
Task: calculate heat capacity for phase change
Programing:
02/2008 JOD moved from CFluidProperties::SpecificHeatCapacity()
last modification:
**************************************************************************/
double CFluidProperties::PhaseChange()
{
   int gueltig = -1;
   double heat_capacity_phase_change = 0;
   double pressure, saturation, temperature;
   double density_vapor, humi, drdT;
   double H1,H0,T0,T_1;                           //T1 defined at 662 in steam67???

   //......................................................................
   CalPrimaryVariable(specific_heat_capacity_pcs_name_vector);
   pressure = primary_variable[0];
   temperature =  primary_variable[1];
   saturation = primary_variable[2];

   if(heat_capacity_model == 3)
   {

      T_1 = primary_variable_t1[1];
      if(T_1 <= T_Latent1 || T_1 >= T_Latent2)
         heat_capacity_phase_change = GetCurveValue(heat_phase_change_curve,0,temperature_buffer,&gueltig);
      else
      {
         heat_capacity_model = 5;                 // ??? JOD
         H1 = CalcEnthalpy(T_1);
         T0 = primary_variable_t0[1];
         if(fabs(T_1-T0)<1.0e-8)
            T_1 +=1.0e-8;
         H0 = CalcEnthalpy(T0);
         heat_capacity_phase_change = (H1-H0)/(T_1-T_0);
         heat_capacity_model = 3;
      }

   }
   else if(heat_capacity_model == 4)
   {

      T_1 = primary_variable_t1[1];
      if(T_1 <= T_Latent1 || T_1 >= T_Latent2)
      {

         humi = exp( pressure /( GAS_CONSTANT_V * temperature_buffer * Density() ) );
         density_vapor = humi * Density();
         drdT = ( vaporDensity_derivative( temperature_buffer )* humi \
            - density_vapor * pressure / ( GAS_CONSTANT_V * Density() * pow( temperature_buffer, 2.0 ) ) ) / Density();
         H1 =  latent_heat + specific_heat_capacity * ( temperature_buffer - T_Latent1);
         heat_capacity_phase_change = H1*drdT;
      }
      else
      {
         heat_capacity_model = 5;
         H1 = CalcEnthalpy(T_1);
         T0 = primary_variable_t0[1];
         if(fabs(T_1-T0)<1.0e-8)
            T_1 +=1.0e-8;
         H0 = CalcEnthalpy(T0);
         heat_capacity_phase_change = (H1-H0)/(T_1-T0);
         heat_capacity_model = 4;
      }

   }

   return heat_capacity_phase_change;

}


/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
02/2007 WW MFP implementation based on MATCalcFluidHeatCapacity (OK)
**************************************************************************/
double MFPCalcFluidsHeatCapacity(CFiniteElementStd* assem)
{
   double heat_capacity_fluids=0.0;
   double PG=0.0, Sw = 0.0,TG,rhow,rho_gw,p_gw,dens_aug[3],rho_g;
   double dens_arg[3];                            //AKS
   //
   CFluidProperties *m_mfp = NULL;
   CRFProcess* m_pcs = assem->cpl_pcs;
   //AKS
                                                  //rho(p,T) //AKS
   if(assem->FluidProp->density_model==14 && assem->MediaProp->heat_diffusion_model==273 && assem->cpl_pcs )
   {
                                                  // pressure
      dens_arg[0]=assem->interpolate(assem->NodalValC1);
                                                  // temperature
      dens_arg[1]=assem->interpolate(assem->NodalVal1)+T_KILVIN_ZERO;
      dens_arg[2]=assem->Index;                   //ELE index
      heat_capacity_fluids = assem->FluidProp->Density(dens_arg)*assem->FluidProp->SpecificHeatCapacity(dens_arg);
   }
   else
   {
      //
      //if (m_pcs->pcs_type_name.find("MULTI_PHASE_FLOW")!=string::npos)
      if (m_pcs->type==1212)                      // non-isothermal multi-phase flow
      {
                                                  // Capillary pressure
         PG = assem->interpolate(assem->NodalValC1);
         Sw = assem->MediaProp->SaturationCapillaryPressureFunction(PG,0);
         double PG2 = assem->interpolate(assem->NodalVal_p2);
         TG=assem->interpolate(assem->NodalVal1)+T_KILVIN_ZERO;
         rhow=assem->FluidProp->Density();
         rho_gw = assem->FluidProp->vaporDensity(TG)*exp(-PG*COMP_MOL_MASS_WATER/(rhow*GAS_CONSTANT*TG));
         p_gw = rho_gw*GAS_CONSTANT*TG/COMP_MOL_MASS_WATER;
         dens_aug[0] = PG2-p_gw;
         dens_aug[1] = TG;
         m_mfp = mfp_vector[1];
                                                  // 2 Dec 2010 AKS
         rho_g = rho_gw + m_mfp->Density(dens_aug);
         //double rho_g = PG2*COMP_MOL_MASS_AIR/(GAS_CONSTANT*(assem->TG+273.15));\\WW
         //
         m_mfp = mfp_vector[0];
         heat_capacity_fluids = Sw * m_mfp->Density() * m_mfp->SpecificHeatCapacity();
         m_mfp = mfp_vector[1];
         heat_capacity_fluids += (1.0-Sw) * rho_g * m_mfp->SpecificHeatCapacity();
      }

      else
      {
         heat_capacity_fluids = assem->FluidProp->Density() * assem->FluidProp->SpecificHeatCapacity();

         if(m_pcs->type != 1)                     // neither liquid nor ground water flow
         {

                                                  //  pressure
            PG = assem->interpolate(assem->NodalValC1);

            if(PG < 0.0)
            {
               Sw = assem->MediaProp->SaturationCapillaryPressureFunction(-PG,0);
               heat_capacity_fluids *= Sw;
               if( assem->GasProp != 0)
                  heat_capacity_fluids += (1.-Sw) * assem->GasProp->Density() * assem->GasProp->SpecificHeatCapacity();
               heat_capacity_fluids += (1.-Sw) * assem->FluidProp->PhaseChange();
            }

         }

      }
   }
   return heat_capacity_fluids;
}


/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
08/2004 OK MFP implementation based on MATCalcFluidHeatCapacity (OK)
10/2005 YD/OK: general concept for heat capacity
overloaded function, see above, taken out, CMCD
**************************************************************************/
/*double MFPCalcFluidsHeatCapacity(double temperature, CFiniteElementStd* assem)
{
  double saturation_phase;
  double heat_capacity_fluids=0.0;
  int nidx0,nidx1;
  //--------------------------------------------------------------------
  // MMP medium properties
  bool New = false; // To be removed. WW
  if(fem_msh_vector.size()>0) New = true;
  //----------------------------------------------------------------------
  CFluidProperties *m_mfp = NULL;
int no_fluids =(int)mfp_vector.size();
switch(no_fluids){
//....................................................................
case 1:
//m_mfp = mfp_vector[0];
if(New) //WW
{
heat_capacity_fluids = assem->FluidProp->Density() \
* assem->FluidProp->SpecificHeatCapacity();
}
else
heat_capacity_fluids = assem->FluidProp->Density() \
* assem->FluidProp->SpecificHeatCapacity();

break;
//....................................................................
case 2:
if(New) //WW
{
nidx0 = GetNodeValueIndex("SATURATION1");
nidx1 = nidx0+1;
saturation_phase = (1.-assem->pcs->m_num->ls_theta)*assem->interpolate(nidx0)
+ assem->pcs->m_num->ls_theta*assem->interpolate(nidx1);
}
else
{
nidx0 = PCSGetNODValueIndex("SATURATION1",0);
nidx1 = PCSGetNODValueIndex("SATURATION1",1);
saturation_phase = (1.-assem->pcs->m_num->ls_theta)*assem->interpolate(nidx0) \
+ assem->pcs->m_num->ls_theta*assem->interpolate(nidx1);
}
m_mfp = mfp_vector[0];
heat_capacity_fluids = saturation_phase \
* assem->FluidProp->Density() \
* assem->FluidProp->SpecificHeatCapacity();
m_mfp = mfp_vector[1];
heat_capacity_fluids += (1.0-saturation_phase) \
* assem->FluidProp->Density() \
* assem->FluidProp->SpecificHeatCapacity();
break;
//....................................................................
case 3: // Entropy based
break;
//....................................................................
default:
cout << "Error in MFPCalcFluidsHeatCapacity: no fluid phase data" << endl;
}
return heat_capacity_fluids;
}*/

/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
08/2004 OK MFP implementation based on MATCalcFluidHeatCapacity (OK)
11/2005 YD Modification
last modification:
NB 4.9.05
**************************************************************************/
                                                  //NB Dec 08 4.9.05
double CFluidProperties::HeatConductivity(double *variables)
{
   int fct_number = 0;
   int gueltig;

   if(variables)                                  //NB Dec 08
   {
      primary_variable[0] = variables[0];         //p (single phase)
      primary_variable[1] = variables[1];         //T (temperature)
      primary_variable[2] = variables[2];         //ELE index
   }
   else
   {
      CalPrimaryVariable(heat_conductivity_pcs_name_vector);
   }

   switch(heat_conductivity_model)
   {
      case 0:                                     // rho = f(x)
         heat_conductivity = GetCurveValue(fct_number,0,primary_variable[0],&gueltig);
         break;
      case 1:                                     // c = const
         heat_conductivity = heat_conductivity;
         break;
      case 2:
         heat_conductivity = MATCalcHeatConductivityMethod2(primary_variable[0],primary_variable[1], primary_variable[2]);
         break;
      case 3:                                     // NB
         heat_conductivity = Fluid_Heat_Conductivity (Density(),primary_variable[1],fluid_id);
         break;
      case 10:                                    // AKS for real gases mixture
                                                  // mixture k= sum_i sum_j x_i*x_j*intrc*sqrt[k_i(p,T)*k_j(p,T)]
         heat_conductivity = MixtureSubProperity(6, (long)  primary_variable[2], primary_variable[0], primary_variable[1]);
         break;
   }
   return heat_conductivity;
}


/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programing:
08/2004 OK MFP implementation based on MATCalcFluidHeatCapacity (OK)
last modification:
10/2010 TF changed access to process type
**************************************************************************/
double MFPCalcFluidsHeatConductivity(long index,double*gp,double theta, CFiniteElementStd* assem)
{
   gp = gp;                                       //OK411
   index = index;
   double saturation_phase = 0.0;                 //OK411
   double heat_conductivity_fluids=0.0;
   int nidx0,nidx1;
   bool New = false;                              // To be removed. WW
   if(fem_msh_vector.size()>0) New = true;

   //--------------------------------------------------------------------
   //----------------------------------------------------------------------
   CFluidProperties *m_mfp = NULL;
   int no_fluids =(int)mfp_vector.size();
   //YD-----------
   CRFProcess* m_pcs = NULL;
   for(int i=0;i<(int)pcs_vector.size();i++)
   {
      m_pcs = pcs_vector[i];
      //    if(m_pcs->pcs_type_name.find("RICHARDS_FLOW"))
      if(m_pcs->getProcessType () == RICHARDS_FLOW)
         no_fluids =1;
   }
   //YD-----------
   switch(no_fluids)
   {
      case 1:
         m_mfp = mfp_vector[0];
         heat_conductivity_fluids = m_mfp->HeatConductivity();
         break;
      case 2:
         m_mfp = mfp_vector[0];
         if(New)                                  //WW
         {
            nidx0 = m_pcs->GetNodeValueIndex("SATURATION1");
            nidx1 = nidx0+1;
            saturation_phase = (1.-theta)*assem->interpolate(nidx0,m_pcs)
               + theta*assem->interpolate(nidx1,m_pcs);

         }
         heat_conductivity_fluids = saturation_phase
            * m_mfp->HeatConductivity();
         m_mfp = mfp_vector[1];
         heat_conductivity_fluids += (1.0-saturation_phase )
            * m_mfp->HeatConductivity();
         break;
      default:
         cout << "Error in MFPCalcFluidsHeatConductivity: no fluid phase data" << endl;
   }
   return heat_conductivity_fluids;
}


////////////////////////////////////////////////////////////////////////////
// Fluid phase change properties
#ifdef obsolete                                   //WW
/**************************************************************************
FEMLib-Method:
Task: Vapour pressure from table
Programing:
03/2002 OK/JdJ Implementation
last modification:
**************************************************************************/
double MFPCalcVapourPressure(double temperature)
{
   double temperature_F;
   double pressure;
   double quality;
   double weight;
   double enthalpy;
   double entropy;
   double saturation_temperature;
   double saturation_pressure;
   double degrees_superheat;
   double degrees_subcooling;
   double viscosity;
   double critical_velocity;
   int action=0;
   double pressure_vapour;

   /*
     double vapour_pressure;
     double vapour_enthalpy = 2258.0; kJ/kg
     double potenz;
     potenz = ((1./temperature_ref)-(1./(*temperature))) * \
                     ((vapour_enthalpy*comp_mol_mass_water)/gas_constant);
     vapour_pressure = pressure_ref * exp(potenz);
   */
   pressure = 1.e-3;                              /*Vorgabe eines vernueftigen Wertes*/
   pressure *= PA2PSI;                            /* Umrechnung Pa in psia */
   temperature -= 273.15;                         /* Kelvin -> Celsius */
   temperature_F = temperature*1.8+32.;           /* Umrechnung Celsius in Fahrenheit */
   steam67 (&temperature_F,
      &pressure,
      &quality,
      &weight,
      &enthalpy,
      &entropy,
      &saturation_temperature,
      &saturation_pressure,
      &degrees_superheat,
      &degrees_subcooling,
      &viscosity,
      &critical_velocity,
      action);
   pressure_vapour = saturation_pressure * PSI2PA;/* Umrechnung psia in Pa */
   /*
   density_vapour = 0.062427962/weight; Dichte in kg/m^3
   enthalpy_vapour = enthalpy * 1055. / 0.454; Umrechnung von btu/lbm in J/kg
   */
   return pressure_vapour;
}


/**************************************************************************
FEMLib-Method:
Task: Get phase and species related enthalpies
Programing:
03/2002 OK/JdJ Implementation
08/2004 OK MFP implementation
last modification:
**************************************************************************/
double CFluidProperties::Enthalpy(int comp,double temperature)
{
   double temperature_F;
   double pressure;
   double quality;
   double weight;
   double enthalpy=0.0;
   double entropy;
   double saturation_temperature;
   double saturation_pressure;
   double degrees_superheat;
   double degrees_subcooling;
   double viscosity;
   double critical_velocity;
   int action=0;

   if((phase==0)&&(comp==0))
   {
      enthalpy = 733.0*temperature + (GAS_CONSTANT*(temperature+0.0))/COMP_MOL_MASS_AIR;
   }
   else if((phase==0)&&(comp==1))                 /* h_w^g: water species in gaseous phase */
   {
      pressure = 1.e-3;                           /*Vorgabe eines vernuenftigen Wertes */
      pressure *= PA2PSI;                         /* Umrechnung Pa in psia */
      temperature -= 273.15;                      /* Kelvin -> Celsius */
      temperature_F = temperature*1.8+32.;        /* Umrechnung Celsius in Fahrenheit*/
      steam67 (&temperature_F,
         &pressure,
         &quality,
         &weight,
         &enthalpy,
         &entropy,
         &saturation_temperature,
         &saturation_pressure,
         &degrees_superheat,
         &degrees_subcooling,
         &viscosity,
         &critical_velocity,
         action);
      enthalpy = enthalpy * 1055. / 0.454;        /* Umrechnung von btu/lbm in J/kg */
   }
   else if((phase==1)&&(comp==0))                 /* h_a^l: air species in liquid phase */
   {
   }
   else if((phase==1)&&(comp==1))                 /* h_w^l: water species in liquid phase */
   {
   }
   return enthalpy;
}


/**************************************************************************
FEMLib-Method:
Task: Get phase and species related enthalpies
Programing:
03/2002 OK/JdJ Implementation
08/2004 OK MFP implementation
last modification:
**************************************************************************/
double CFluidProperties::EnthalpyPhase(long number,int comp,double*gp,double theta)
{
   double temperature;
   double enthalpy=0.0;
   double mass_fraction_air,enthalpy_air,mass_fraction_water,enthalpy_water;
   int nidx0,nidx1;

   nidx0 = PCSGetNODValueIndex("TEMPERATURE1",0);
   nidx1 = PCSGetNODValueIndex("TEMPERATURE1",1);
   temperature = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
      + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);

   if(phase==0)
   {
      if(comp<0)
      {
         comp=0;
         mass_fraction_air = MassFraction(number,comp,gp,theta);
         comp=1;
         mass_fraction_water = MassFraction(number,comp,gp,theta);
         comp=0;
         enthalpy_air = Enthalpy(comp,temperature);
         comp=1;
         enthalpy_water = Enthalpy(comp,temperature);

         enthalpy = mass_fraction_air * enthalpy_air + \
            mass_fraction_water * enthalpy_water;
      }
      else if (comp>=0)
         enthalpy = Enthalpy(comp,temperature);
   }
   else if(phase==1)
   {
                                                  //
      enthalpy = SpecificHeatCapacity() * temperature;
   }

   return enthalpy;
}


/**************************************************************************
FEMLib-Method:
Task:
   Calculate Henry constant
Programing:
02/2002 JdJ First implementation
04/2003 JdJ temperature conversion from celsius to kelvin
last modification:
**************************************************************************/
double MFPCalcHenryConstant(double temperature)
{
   double HenryConstant;
   HenryConstant =  (0.8942 +1.47*exp(-0.04394*(temperature-273.14)))*0.0000000001;
   return HenryConstant;
}


/**************************************************************************
FEMLib-Method:
Task:
   Calculate mass fractions
   of gas phase    X^g_a, X^g_w according to Claudius-Clapeyron law
   of liquid phase X^l_a, X^l_w according to Henry law
Programing:
01/2002 OK/JdJ First implementation
03/2002 JdJ Gas phase in argument for pressure in mass fraction.
04/2003 JdJ Dichteberechnung (setzt voraus, das Phase 0 gas ist).
04/2003 JdJ Druckberechnung (setzt voraus, das Phase 0 gas ist).
04/2004 JdJ Bugfix Gauss Points
08/2004 OK MFP implementation
last modification:
**************************************************************************/
double CFluidProperties::MassFraction(long number,int comp,double*gp,double theta,CFiniteElementStd* assem)
{
   double mass_fraction=0.0;
   double mass_fraction_air_in_gas,mass_fraction_air_in_liquid;
   double gas_density=0.0;
   double vapour_pressure;
   double temperature;
   double henry_constant;
   int nidx0,nidx1;
   double gas_pressure;
   /*--------------------------------------------------------------------------*/
   /* Get and calc independent variables */
   nidx0 = PCSGetNODValueIndex("PRESSURE1",0);
   nidx1 = PCSGetNODValueIndex("PRESSURE1",1);
   if(mode==0)                                    // Gauss point values
   {
      gas_pressure = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
         + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
   }
   else                                           // Node values
   {
      gas_pressure = (1.-theta)*GetNodeVal(number,nidx0) \
         + theta*GetNodeVal(number,nidx1);
   }
   nidx0 = PCSGetNODValueIndex("TEMPERATURE1",0);
   nidx1 = PCSGetNODValueIndex("TEMPERATURE1",1);
   if(mode==0)                                    // Gauss point values
   {
      temperature = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
         + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
   }
   else                                           // Node values
   {
      temperature = (1.-theta)*GetNodeVal(number,nidx0) \
         + theta*GetNodeVal(number,nidx1);
   }
   gas_density = Density();
   vapour_pressure = MFPCalcVapourPressure(temperature);
   /*--------------------------------------------------------------------------*/
   /* Calc mass fractions */
   switch (phase)
   {
      case 0:                                     /* gas phase */
         mass_fraction_air_in_gas = \
            ((gas_pressure-vapour_pressure)*COMP_MOL_MASS_AIR) \
            / (GAS_CONSTANT*(temperature+0.0)*gas_density);
         mass_fraction_air_in_gas = MRange(0.0,mass_fraction_air_in_gas,1.0);
         if(comp==0)                              /* air specie */
         {
            mass_fraction = mass_fraction_air_in_gas;
         }
         if(comp==1)                              /* water specie */
         {
            mass_fraction = 1.0 - mass_fraction_air_in_gas;
         }
         break;
      case 1:                                     /* liquid phase */
         henry_constant = MFPCalcHenryConstant(temperature);
         mass_fraction_air_in_liquid = \
            COMP_MOL_MASS_AIR / (COMP_MOL_MASS_AIR \
            - COMP_MOL_MASS_WATER * (1.0-1.0/(henry_constant*(gas_pressure-vapour_pressure))));
         mass_fraction_air_in_liquid = MRange(0.0,mass_fraction_air_in_liquid,1.0);
         if(comp==0)                              /* air specie */
         {
            mass_fraction = mass_fraction_air_in_liquid;
         }
         if(comp==1)                              /* water specie X_w^l = 1-X_a^l */
         {
            mass_fraction = 1.0 - mass_fraction_air_in_liquid;
         }
         break;
   }
   /*--------------------------------------------------------------------------*/
   return mass_fraction;
}


/**************************************************************************
FEMLib-Method:
Task:
   Calculate Henry constant
Programing:
02/2002 OK/JdJ Implementation
08/2004 OK MFP implementation
last modification:
**************************************************************************/
double CFluidProperties::InternalEnergy(long number,double*gp,double theta)
{
   double energy=0.0;
   int nidx0,nidx1;
   double temperature,pressure;

   nidx0 = PCSGetNODValueIndex("TEMPERATURE1",0);
   nidx1 = PCSGetNODValueIndex("TEMPERATURE1",1);
   temperature = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
      + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
   energy = SpecificHeatCapacity() * temperature; //YD
   //energy = HeatCapacity(temperature,m_ele_fem_std) * temperature;

   //if(name.find("GAS")!=string::npos){
   if(phase==0)
   {
      nidx0 = PCSGetNODValueIndex("PRESSURE1",0);
      nidx1 = PCSGetNODValueIndex("PRESSURE1",1);
      pressure = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
         + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
      energy += pressure / Density();             //YD
   }
   return energy;
}


/**************************************************************************
FEMLib-Method:
Task:
   Calculate Henry constant
   temperature in Kelvin
Programing:
01/2003 OK/JdJ Implementation
08/2004 OK MFP implementation
last modification:
**************************************************************************/
double CFluidProperties::DensityTemperatureDependence(long number,int comp,double*gp,double theta)
{
   double vapour_pressure;
   double dvapour_pressure_dT;
   double drho_dT;
   double temperature;
   int nidx0,nidx1;
   //----------------------------------------------------------------------
   // State functions
   nidx0 = PCSGetNODValueIndex("TEMPERATURE1",0);
   nidx1 = PCSGetNODValueIndex("TEMPERATURE1",1);
   if(mode==0)                                    // Gauss point values
   {
      temperature = (1.-theta)*InterpolValue(number,nidx0,gp[0],gp[1],gp[2]) \
         + theta*InterpolValue(number,nidx1,gp[0],gp[1],gp[2]);
   }
   else                                           // Node values
   {
      temperature = (1.-theta)*GetNodeVal(number,nidx0) \
         + theta*GetNodeVal(number,nidx1);
   }
   //----------------------------------------------------------------------
   // Vapour
   vapour_pressure = MFPCalcVapourPressure(temperature);
   dvapour_pressure_dT = COMP_MOL_MASS_WATER * Enthalpy(comp,temperature) \
      / (GAS_CONSTANT*temperature*temperature) \
      * vapour_pressure;

   drho_dT = -1.0 * COMP_MOL_MASS_WATER / (GAS_CONSTANT*temperature) \
      * (dvapour_pressure_dT - vapour_pressure/temperature);
   //----------------------------------------------------------------------
   // Test
   if((phase>0)||(comp==0))
   {
      DisplayMsgLn("MATCalcFluidDensityTemperatureDependence: Incorrect use of function");
      abort();
   }
   return drho_dT;
}
#endif                                            // if define obsolete. WW

/**************************************************************************
FEMLib-Method:
Task:

Programing:

last modification:
**************************************************************************/
double CFluidProperties::LiquidViscosity_CMCD(double Press,double TempK,double C)
{
   C = C;
   /*CMcD variables for 20 ALR*/
   double A1,A2,A3,A4,A5,A6,A7,A8;                /*constants*/
   double TempC,TempF, Pbar,Salinity;             /* Temperature [K], Temperature [F], Pressure [bar]*/
   double my_Zero,PsatBar, PsatKPa;               /*my pure water, my saline water, Saturation pressure [bar], Saturation pressure [KPa]*/
                                                  /*intermediate values*/
   double sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8, exponent;
   /*CMcD end variables for 20 ALR*/
   /* //Prepared for introduction of solute transport in PCS version
      //Average Concentration
      comp=0; // nur fuer Einkomponenten-Systeme
     timelevel=1;
     concentration_average = 0.0;
     for (i = 0; i < count_nodes; i++)
       concentration_average += GetNodeVal(element_nodes[i], c_idx);
     concentration_average /= count_nodes;
     fp->salinity=fp->rho_0+concentration_average*fp->drho_dC;
               */
   //Link to function from ALR
   Salinity=C_0/1000.;
   /***constant values*******/
   A1 = -7.419242;
   A2 = -0.29721;
   A3 = -0.1155286;
   A4 = -0.008685635;
   A5 = 0.001094098;
   A6 = 0.00439993;
   A7 = 0.002520658;
   A8 = 0.0005218684;

   /*Unit conversions*/
   TempC = TempK-273.15;
   TempF = TempC*1.8 + 32.0;
   Pbar = Press/100000.0;
   /*end of units conversion*/

   /*Calculation of the saturation pressure*/
   sum1=pow((0.65-0.01*TempK),0)*A1;
   sum2=pow((0.65-0.01*TempK),1)*A2;
   sum3=pow((0.65-0.01*TempK),2)*A3;
   sum4=pow((0.65-0.01*TempK),3)*A4;
   sum5=pow((0.65-0.01*TempK),4)*A5;
   sum6=pow((0.65-0.01*TempK),5)*A6;
   sum7=pow((0.65-0.01*TempK),6)*A7;
   sum8=pow((0.65-0.01*TempK),7)*A8;

                                                  /*intermediate value*/
   exponent = sum1+sum2+sum3+sum4+sum5+sum6+sum7+sum8;
   exponent= exponent*(374.136-TempC)/TempK;      /*intermediate value*/

   PsatKPa = exp(exponent)*22088;                 /*saturation pressure in kPa*/
   PsatBar = PsatKPa/(1000*100000);               /*Saturation pressure in bar*/

   /*Viscosity of pure water in Pa-S*/
   my_Zero = 243.18e-7 * (pow(10.,(247.8/(TempK-140)))) * (1+(Pbar-PsatBar)*1.0467e-6 * (TempK-305));

   /*Viscosity of saline water in Pa-S*/
   viscosity = my_Zero * (1-0.00187* (pow(Salinity,0.5)) + 0.000218* (pow(Salinity,2.5))+(pow(TempF,0.5)-0.0135*TempF)*(0.00276*Salinity-0.000344* (pow(Salinity,1.5))));
   return viscosity;
}


/**************************************************************************/
/* ROCKFLOW - Function: MATCalcFluidHeatConductivityMethod2
 */
/* Task:
   Calculate heat conductivity of all fluids
                                                                          */
/* Parameter: (I: Input; R: Return; X: Both)
   I double temperature
                                                                          */
/* Return:
   Value of heat conductivity of fluids as a function of (p,C)
                                                                          */
/* Programming:
   09/2003   CMCD ARL   Implementation
   08/2004   CMCD inclusion in GeoSys v. 4.
                                                                             */
/**************************************************************************/
double CFluidProperties::MATCalcHeatConductivityMethod2(double Press, double TempK, double Conc)
{
   Conc = Conc;
   int i, j;
   double TauTC, PiiTC, Nabla, Delta, Nabla0, Nabla1, Nabla2;
   double heat_conductivity, Rho, temperature_average, pressure_average, viscosity;
   double Rhostar, TstarTC, Lambdastar, Pstar1;
   double nZero[4];
   double n[5][6];
   double A1,A2,A3,A4,A5,A6,A7,A8;                /*constants*/
   double TempC, TempF, Pbar, Salinity;           /* Temperature [K], Temperature [F], Pressure [bar]*/
   double my_Zero,PsatBar, PsatKPa;               /*my pure water, my saline water, Saturation pressure [bar], Saturation pressure [KPa]*/
                                                  /*intermediate values*/
   double sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8, exponent;

   /*************************************************************************************************/
   /*************************************************************************************************/
   /*************************Partial derivatives calculation*****************************************/
   double GammaPi, GammaPiTau, GammaPiPi, Pii, Tau,  GazConst;
   double LGamma[35];
   double JGamma[35];
   double nGamma[35];
   double Tstar, Pstar;
   double TstarTilda, PstarTilda, RhostarTilda;
   double First_derivative, Second_derivative;

   pressure_average = Press;
   temperature_average = TempK;
   Salinity = C_0;
   Tstar = 1386;
   Pstar = 16.53e6;                               // MPa
   GazConst = 0.461526e3;                         //!!!!Given by equation (2.1)
   TstarTilda = 647.226;
   PstarTilda = 22.115e6;
   RhostarTilda = 317.763;
   Lambdastar = 0.4945;
   Rhostar = 317.763;
   TstarTC = 647.226;
   Pstar1 = 22.115e-6;

   /*BEGIN: reduced dimensions*/
   TauTC = TstarTC / temperature_average;
   //Delta = Rho / Rhostar;
   PiiTC = pressure_average / Pstar1;
   /*END: reduced dimensions*/

   nGamma[1] = 0.14632971213167;
   nGamma[2] = -0.84548187169114;
   nGamma[3] = -0.37563603672040e1;
   nGamma[4] = 0.33855169168385e1;
   nGamma[5] = -0.95791963387872;
   nGamma[6] = 0.15772038513228;
   nGamma[7] = -0.16616417199501e-1;
   nGamma[8] = 0.81214629983568e-3;
   nGamma[9] = 0.28319080123804e-3;
   nGamma[10] = -0.60706301565874e-3;
   nGamma[11] = -0.18990068218419e-1;
   nGamma[12] = -0.32529748770505e-1;
   nGamma[13] = -0.21841717175414e-1;
   nGamma[14] = -0.52838357969930e-4;
   nGamma[15] = -0.47184321073267e-3;
   nGamma[16] = -0.30001780793026e-3;
   nGamma[17] = 0.47661393906987e-4;
   nGamma[18] = -0.44141845330846e-5;
   nGamma[19] = -0.72694996297594e-15;
   nGamma[20] = -0.31679644845054e-4;
   nGamma[21] = -0.28270797985312e-5;
   nGamma[22] = -0.85205128120103e-9;
   nGamma[23] = -0.22425281908000e-5;
   nGamma[24] = -0.65171222895601e-6;
   nGamma[25] = -0.14341729937924e-12;
   nGamma[26] = -0.40516996860117e-6;
   nGamma[27] = -0.12734301741641e-8;
   nGamma[28] = -0.17424871230634e-9;
   nGamma[29] = -0.68762131295531e-18;
   nGamma[30] = 0.14478307828521e-19;
   nGamma[31] = 0.26335781662795e-22;
   nGamma[32] = -0.11947622640071e-22;
   nGamma[33] = 0.18228094581404e-23;
   nGamma[34] = -0.93537087292458e-25;

   LGamma[1] =0.;
   LGamma[2] =0.;
   LGamma[3] =0.;
   LGamma[4] =0.;
   LGamma[5] =0.;
   LGamma[6] =0.;
   LGamma[7] =0.;
   LGamma[8] =0.;
   LGamma[9] =1.;
   LGamma[10] =1.;
   LGamma[11] =1.;
   LGamma[12] =1.;
   LGamma[13] =1.;
   LGamma[14] =1.;
   LGamma[15] =2.;
   LGamma[16] =2.;
   LGamma[17] =2.;
   LGamma[18] =2.;
   LGamma[19] =2.;
   LGamma[20] =3.;
   LGamma[21] =3.;
   LGamma[22] =3.;
   LGamma[23] =4.;
   LGamma[24] =4.;
   LGamma[25] =4.;
   LGamma[26] =5.;
   LGamma[27] =8.;
   LGamma[28] =8.;
   LGamma[29] =21.;
   LGamma[30] =23.;
   LGamma[31] =29.;
   LGamma[32] =30.;
   LGamma[33] =31.;
   LGamma[34] =32.;

   JGamma[1] =-2.;
   JGamma[2] =-1.;
   JGamma[3] =0.;
   JGamma[4] =1.;
   JGamma[5] =2.;
   JGamma[6] =3.;
   JGamma[7] =4.;
   JGamma[8] =5.;
   JGamma[9] =-9.;
   JGamma[10] =-7.;
   JGamma[11] =-1.;
   JGamma[12] =0.;
   JGamma[13] =1.;
   JGamma[14] =3.;
   JGamma[15] =-3.;
   JGamma[16] =0.;
   JGamma[17] =1.;
   JGamma[18] =3.;
   JGamma[19] =17.;
   JGamma[20] =-4.;
   JGamma[21] =0.;
   JGamma[22] =6.;
   JGamma[23] =-5.;
   JGamma[24] =-2.;
   JGamma[25] =10.;
   JGamma[26] =-8.;
   JGamma[27] =-11.;
   JGamma[28] =-6.;
   JGamma[29] =-29.;
   JGamma[30] =-31.;
   JGamma[31] =-38.;
   JGamma[32] = -39.;
   JGamma[33] =-40.;
   JGamma[34] =-41.;

   Pii = pressure_average / Pstar;
   Tau = Tstar / temperature_average;

   /*BEGIN:Calculation of GammaPi*/
   GammaPi = 0;

   for (i=1; i<35; i++)
   {
      GammaPi = GammaPi - (nGamma[i]) * (LGamma[i]) * (pow((7.1-Pii),(LGamma[i] -1))) * (pow((Tau-1.222),JGamma[i]));
   }
   /*END: Calculation of GammaPi*/

   /*BEGIN:Calculation of GammaPiTau*/
   GammaPiTau = 0;
   for (i=1; i<35; i++)
   {
      GammaPiTau = GammaPiTau - (nGamma[i]) * (LGamma[i]) * (pow((7.1-Pii),(LGamma[i] -1))) * (JGamma[i]) * (pow((Tau-1.222),(JGamma[i]-1)));
   }
   /*END: Calculation of GammaPiTau*/

   /*BEGIN:Calculation of GammaPiPi*/
   GammaPiPi = 0;
   for (i=1; i<=34; i++)
   {
      GammaPiPi = GammaPiPi + (nGamma[i]) * (LGamma[i]) * (LGamma[i] -1) * (pow((7.1-Pii),(LGamma[i] -2))) * (pow((Tau-1.222),(JGamma[i])));
   }
   /*END: Calculation of GammaPiPi*/

   /*BEGIN:Calculation of derivative*/
   First_derivative = ((TstarTilda) * (Pstar) * ((GammaPiTau)*(Tstar) - (GammaPi) * (temperature_average))) / (PstarTilda * pow(temperature_average,2) * GammaPiPi),
      Second_derivative = ((-1) * (PstarTilda) * (GammaPiPi) ) / ( (RhostarTilda) * (temperature_average) * (GazConst) * (pow(GammaPi,2)));
   /*End:Calculation of derivative*/

   /*BEGIN: Calculation of density*/
   Rho = pressure_average / (GazConst * (temperature_average) * (Pii) * (GammaPi));
   /*END: Calculation of density*/

   /*************************Partial derivatives calculation*****************************************/
   /*************************************************************************************************/
   /*************************************************************************************************/

   /*BEGIN: Constant values*/

   Lambdastar = 0.4945;
   Rhostar = 317.763;
   TstarTC = 647.226;
   Pstar1 = 22.115e6;

   nZero[0] = 0.1e1;
   nZero[1] = 0.6978267e1;
   nZero[2] = 0.2599096e1;
   nZero[3] = -0.998254;

   n[0][0] = 0.13293046e1;
   n[0][1] = -0.40452437;
   n[0][2] = 0.24409490;
   n[0][3] = 0.18660751e-1;
   n[0][4] = -0.12961068;
   n[0][5] = 0.44809953e-1;

   n[1][0] = 0.17018363e1;
   n[1][1] = -0.22156845e1;
   n[1][2] = 0.16511057e1;
   n[1][3] = -0.76736002;
   n[1][4] = 0.37283344;
   n[1][5] = -0.11203160;

   n[2][0] = 0.52246158e1;
   n[2][1] = -0.10124111e2;
   n[2][2] = 0.49874687e1;
   n[2][3] = -0.27297694;
   n[2][4] = -0.43083393;
   n[2][5] = 0.13333849;

   n[3][0] = 0.87127675e1;
   n[3][1] = -0.95000611e1;
   n[3][2] = 0.43786606e1;
   n[3][3] = -0.91783782;
   n[3][4] = 0;
   n[3][5] = 0;

   n[4][0] = -0.18525999e1;
   n[4][1] = 0.93404690;
   n[4][2] = 0;
   n[4][3] = 0;
   n[4][4] = 0;
   n[4][5] = 0;
   /*END: Constant values*/

   /*BEGIN: reduced dimensions*/
   TauTC = TstarTC / temperature_average;
   Delta = Rho / Rhostar;
   PiiTC = pressure_average / Pstar1;
   /*END: reduced dimensions*/

   /*BEGIN: Nabla0*/
   Nabla0 = 0;
   for (i=0; i<=3;i++)
   {
      Nabla0 = Nabla0 + (nZero[i]) * (pow(TauTC,i));
   }

   Nabla0 = Nabla0 * (pow(TauTC, 0.5));
   Nabla0 = 1 / Nabla0;
   /*END: Nabla0*/

   /*BEGIN: Nabla1*/
   Nabla1 = 0;
   for (i=0; i<=4;i++)
   {
      for (j=0; j<=5;j++)
      {
         Nabla1 = Nabla1 + (n[i][j]) * (pow((TauTC-1),i)) * (pow((Delta-1),j));
      }

   }
   Nabla1 = Delta * (Nabla1);
   Nabla1 = exp(Nabla1);
   /*END: Nabla1*/

   /*Calculate Viscosity*/
   /*Link to function from ALR*/

   TempK=temperature_average;
   Press=pressure_average;
   Salinity=C_0/1000.;

   /***constant values*******/
   A1 = -7.419242;
   A2 = -0.29721;
   A3 = -0.1155286;
   A4 = -0.008685635;
   A5 = 0.001094098;
   A6 = 0.00439993;
   A7 = 0.002520658;
   A8 = 0.0005218684;

   /*Unit conversions*/
   TempC = TempK-273.15;
   TempF = TempC*1.8 + 32.0;
   if (TempF < 0.0)
   {
      TempF = 0.0;
   }
   Pbar = Press/100000.0;
   /*end of units conversion*/

   /*Calculation of the saturation pressure*/
   sum1=pow((0.65-0.01*TempK),0)*A1;
   sum2=pow((0.65-0.01*TempK),1)*A2;
   sum3=pow((0.65-0.01*TempK),2)*A3;
   sum4=pow((0.65-0.01*TempK),3)*A4;
   sum5=pow((0.65-0.01*TempK),4)*A5;
   sum6=pow((0.65-0.01*TempK),5)*A6;
   sum7=pow((0.65-0.01*TempK),6)*A7;
   sum8=pow((0.65-0.01*TempK),7)*A8;

                                                  /*intermediate value*/
   exponent = sum1+sum2+sum3+sum4+sum5+sum6+sum7+sum8;
   exponent= exponent*(374.136-TempC)/TempK;      /*intermediate value*/

   PsatKPa = exp(exponent)*22088;                 /*saturation pressure in kPa*/
   PsatBar = PsatKPa/(1000*100000);               /*Saturation pressure in bar*/

   /*Viscosity of pure water in Pa-S*/
   my_Zero = 243.18e-7 * (pow(10.,(247.8/(TempK-140)))) * (1+(Pbar-PsatBar)*1.0467e-6 * (TempK-305));

   /*Viscosity of saline water in Pa-S*/
   viscosity = my_Zero * (1-0.00187* (pow(Salinity,0.5)) + 0.000218* (pow(Salinity,2.5))+(pow(TempF,0.5)-0.0135*TempF)*(0.00276*Salinity-0.000344* (pow(Salinity,1.5))));

   /* End of viscosity function*/

   /*BEGIN: Nabla2*/
   Nabla2 = 0.0013848 / ((viscosity)/55.071e-6) * (pow(((TauTC)*(Delta)),(-2))) * (pow(First_derivative,2)) * (pow((Delta * (Second_derivative)),0.4678)) * (pow(Delta,0.5)) * exp(-18.66 * (pow((1/TauTC-1),2)) - (pow(Delta-1,4)));
   /*END: Nabla2*/

   /*BEGIN: Nabla => heat_conductivity*/
   Nabla = Nabla0 * (Nabla1) + Nabla2;
   heat_conductivity = Nabla * (Lambdastar);
   /*END: Nabla => Lambda*/

   return heat_conductivity;
}


/**************************************************************************
ROCKFLOW - Funktion:
Task:
Programming:
 08/2005 WW Implementation
**************************************************************************/
double CFluidProperties::vaporDensity(const double T_abs)
{
   return 1.0e-3*exp(19.81-4975.9/T_abs);
}


double CFluidProperties::vaporDensity_derivative(const double T_abs)
{
   return 4.9759*exp(19.81-4975.9/T_abs)/(T_abs*T_abs);
}


/**************************************************************************/
/* ROCKFLOW - Function: MATCalcFluidHeatCapacity
 */
/* Task:
   Calculate heat capacity of all fluids
                                                                          */
/* Parameter: (I: Input; R: Return; X: Both)
   I double temperature
                                                                          */
/* Return:
   Value of heat capacity of all fluids
                                                                          */
/* Programming:
   09/2003   CMCD ARL   Implementation
   09/2004   CMCD included in Geosys ver. 4

                                                                          */
/**************************************************************************/
double CFluidProperties::MATCalcFluidHeatCapacityMethod2(double Press, double TempK, double Conc)
{
   Conc = Conc;
   double Pressurevar, Tau, pressure_average, temperature_average, Tstar, Pstar,GazConst;
   double GammaPi, GammaPiTau, GammaPiPi, GammaTauTau;
   double L[35],J[35],n[35];
   int i;
   double salinity;
   double Cp, Cv;

   pressure_average = Press;
   temperature_average = TempK;

   salinity=C_0;

   Tstar = 1386;

   Pstar = 16.53e6;                               // MPa
   GazConst = 0.461526e3;                         //

   n[0] = 0.0;
   n[1] = 0.14632971213167;
   n[2] = -0.84548187169114;
   n[3] = -0.37563603672040e1;
   n[4] = 0.33855169168385e1;
   n[5] = -0.95791963387872;
   n[6] = 0.15772038513228;
   n[7] = -0.16616417199501e-1;
   n[8] = 0.81214629983568e-3;
   n[9] = 0.28319080123804e-3;
   n[10] = -0.60706301565874e-3;
   n[11] = -0.18990068218419e-1;
   n[12] = -0.32529748770505e-1;
   n[13] = -0.21841717175414e-1;
   n[14] = -0.52838357969930e-4;
   n[15] = -0.47184321073267e-3;
   n[16] = -0.30001780793026e-3;
   n[17] = 0.47661393906987e-4;
   n[18] = -0.44141845330846e-5;
   n[19] = -0.72694996297594e-15;
   n[20] = -0.31679644845054e-4;
   n[21] = -0.28270797985312e-5;
   n[22] = -0.85205128120103e-9;
   n[23] = -0.22425281908000e-5;
   n[24] = -0.65171222895601e-6;
   n[25] = -0.14341729937924e-12;
   n[26] = -0.40516996860117e-6;
   n[27] = -0.12734301741641e-8;
   n[28] = -0.17424871230634e-9;
   n[29] = -0.68762131295531e-18;
   n[30] = 0.14478307828521e-19;
   n[31] = 0.26335781662795e-22;
   n[32] = -0.11947622640071e-22;
   n[33] = 0.18228094581404e-23;
   n[34] = -0.93537087292458e-25;

   L[0] =0.;
   L[1] =0.;
   L[2] =0.;
   L[3] =0.;
   L[4] =0.;
   L[5] =0.;
   L[6] =0.;
   L[7] =0.;
   L[8] =0.;
   L[9] =1.;
   L[10] =1.;
   L[11] =1.;
   L[12] =1.;
   L[13] =1.;
   L[14] =1.;
   L[15] =2.;
   L[16] =2.;
   L[17] =2.;
   L[18] =2.;
   L[19] =2.;
   L[20] =3.;
   L[21] =3.;
   L[22] =3.;
   L[23] =4.;
   L[24] =4.;
   L[25] =4.;
   L[26] =5.;
   L[27] =8.;
   L[28] =8.;
   L[29] =21.;
   L[30] =23.;
   L[31] =29.;
   L[32] =30.;
   L[33] =31.;
   L[34] =32.;

   J[0] =-2.;
   J[1] =-2.;
   J[2] =-1.;
   J[3] =0.;
   J[4] =1.;
   J[5] =2.;
   J[6] =3.;
   J[7] =4.;
   J[8] =5.;
   J[9] =-9.;
   J[10] =-7.;
   J[11] =-1.;
   J[12] =0.;
   J[13] =1.;
   J[14] =3.;
   J[15] =-3.;
   J[16] =0.;
   J[17] =1.;
   J[18] =3.;
   J[19] =17.;
   J[20] =-4.;
   J[21] =0.;
   J[22] =6.;
   J[23] =-5.;
   J[24] =-2.;
   J[25] =10.;
   J[26] =-8.;
   J[27] =-11.;
   J[28] =-6.;
   J[29] =-29.;
   J[30] =-31.;
   J[31] =-38.;
   J[32] = -39.;
   J[33] =-40.;
   J[34] =-41.;

   Pressurevar = pressure_average / Pstar;
   Tau = Tstar / temperature_average;

   /*BEGIN:Calculation of GammaPi*/
   GammaPi = 0;

   for (i=1; i<35; i++)
   {
      GammaPi = GammaPi - (n[i]) * (L[i]) * (pow((7.1-Pressurevar),(L[i] -1.))) * (pow((Tau-1.222),J[i]));
   }
   /*END: Calculation of GammaPi*/

   /*BEGIN:Calculation of GammaPiTau*/
   GammaPiTau = 0;
   for (i=1; i<35; i++)
   {
      GammaPiTau = GammaPiTau - (n[i]) * (L[i]) * (pow((7.1-Pressurevar),(L[i] -1.))) * (J[i]) * (pow((Tau-1.222),(J[i]-1.)));
   }
   /*END: Calculation of GammaPiTau*/

   /*BEGIN:Calculation of GammaTauTau*/
   GammaTauTau = 0;
   for (i=1; i<35; i++)
   {
      GammaTauTau = GammaTauTau + (n[i]) * (pow((7.1-Pressurevar),(L[i]))) * (J[i]) * (J[i] -1.) * (pow((Tau - 1.222),(J[i]-2)));
   }
   /*END: Calculation of GammaTauTau*/

   /*BEGIN:Calculation of GammaPiPi*/
   GammaPiPi = 0;
   for (i=1; i<35; i++)
   {
      GammaPiPi = GammaPiPi + (n[i]) * (L[i]) * (L[i] -1) * (pow((7.1-Pressurevar),(L[i] -2.))) * (pow((Tau-1.222),(J[i])));
   }
   /*END: Calculation of GammaPiPi*/

   /*************************Partial derivatives calculation*****************************************/
   /*************************************************************************************************/
   /*************************************************************************************************/

   /*BEGIN: Fluid isobaric heat capacity*/
   Cp = - (pow(Tau,2))* (GammaTauTau) * (GazConst);
   /*END: Fluid isobaric heat capacity*/

   /*BEGIN: Fluid isochoric heat capacity*/
                                                  /* Cv is not used currently 9.2003*/
   Cv = (- (pow(Tau,2))* (GammaTauTau) + pow(GammaPi - Tau * (GammaPiTau),2) / GammaPiPi) * GazConst;
   /*BEGIN: Fluid isochoric heat capacity*/

   return Cp;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
last modified:
**************************************************************************/
void MFPDelete()
{
   long i;
   int no_mfp =(int)mfp_vector.size();
   for(i=0;i<no_mfp;i++)
   {
      delete mfp_vector[i];
   }
   mfp_vector.clear();
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
10/2005 OK/YD Implementation
**************************************************************************/
CFluidProperties* MFPGet(const string &name)
{
   CFluidProperties* m_mfp = NULL;
   for(int i=0;i<(int)mfp_vector.size();i++)
   {
      m_mfp = mfp_vector[i];
      if(m_mfp->name.compare(name)==0)
         return m_mfp;
   }
   return NULL;
}


CFluidProperties* MFPGet(int fluid)               //NB
{
   CFluidProperties* m_mfp = NULL;
   for(int i=0;i<(int)mfp_vector.size();i++)
   {
      m_mfp = mfp_vector[i];
      if(m_mfp->fluid_id==fluid)
         return m_mfp;
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
10/2005 YD Calculate Enthalpy
**************************************************************************/
double CFluidProperties::CalcEnthalpy(double temperature)
{
   //----------------------------------------------------------------------
   CalPrimaryVariable(enthalpy_pcs_name_vector);
   //----------------------------------------------------------------------
   //double temperature = primary_variable[0];
   double val = 0.0;
   double T0_integrate = 0.0;
   double T1_integrate = 0.0;
   double heat_capacity_all;

   switch(GetHeatCapacityModel())
   {
      case 5:
         MFPGet("LIQUID");
         //------------PART 1--------------------------------
         T0_integrate = 273.15;
         T1_integrate = T_Latent1+273.15;
         int npoint = 100;                        //Gauss point
         double DT = (T1_integrate-T0_integrate)/npoint;
         for(int i=0;i<npoint;i++)
         {
            temperature_buffer = T0_integrate+i*DT-273.15;
            heat_capacity_all = Fem_Ele_Std->FluidProp->SpecificHeatCapacity();
            temperature_buffer = T0_integrate+(i+1)*DT-273.15;
            heat_capacity_all += Fem_Ele_Std->FluidProp->SpecificHeatCapacity();
            val += 0.5*DT*heat_capacity_all;
         }

         //------------PART 2--------------------------------
         npoint = 500;
         T0_integrate = T_Latent1+273.15;
         T1_integrate = temperature+273.15;
         DT = (T1_integrate-T0_integrate)/npoint;
         for(int i=0;i<npoint;i++)
         {
            temperature_buffer = T0_integrate+i*DT-273.15;
            heat_capacity_all = Fem_Ele_Std->FluidProp->SpecificHeatCapacity();
            temperature_buffer = T0_integrate+(i+1)*DT-273.15;
            heat_capacity_all += Fem_Ele_Std->FluidProp->SpecificHeatCapacity();
            val += 0.5*DT*heat_capacity_all;
         }
         break;
         //  case 5:
   }
   return val;
}


/**************************************************************************
PCSLib-Method:
08/2008 OK
last change: 11/2008 NB
**************************************************************************/
double MFPGetNodeValue(long node,const string &mfp_name, int phase_number)
{
   double mfp_value = 0.0;                        //OK411
   //  char c;
   double arguments[2];
   string pcs_name1;
   string pcs_name2;
   CRFProcess *tp;
                                                  //NB
   CFluidProperties *m_mfp = mfp_vector[max(phase_number,0)];

   int mfp_id = -1;
   int val_idx=0;                                 // for later use, NB case 'V': mfp_id = 0; //VISCOSITY
   switch (mfp_name[0])
   {
      case 'V': mfp_id = 0;                       //VISCOSITY
      if(m_mfp->viscosity_pcs_name_vector.size()<1)
         {pcs_name1 = "PRESSURE1";}
         else
      {
         pcs_name1 = m_mfp->viscosity_pcs_name_vector[0];
      }
      if(m_mfp->viscosity_pcs_name_vector.size()<2)
         {pcs_name2 = "TEMPERATURE1";}
         else
      {
         pcs_name2 = m_mfp->viscosity_pcs_name_vector[1];
      }
      break;
      case 'D': mfp_id = 1;                       //DENSITY
      if(m_mfp->density_pcs_name_vector.size()<1)
         {pcs_name1 = "PRESSURE1";}
         else
      {
         pcs_name1 = m_mfp->density_pcs_name_vector[0];
      }
      if(m_mfp->density_pcs_name_vector.size()<2)
         {pcs_name2 = "TEMPERATURE1";}
         else
      {
         pcs_name2 = m_mfp->density_pcs_name_vector[1];
      }
      break;
      case 'H': mfp_id = 2;                       //HEAT_CONDUCTIVITY
      if(m_mfp->heat_conductivity_pcs_name_vector.size()<1)
         {pcs_name1 = "PRESSURE1";}
         else
      {
         pcs_name1 = m_mfp->heat_conductivity_pcs_name_vector[0];
      }
      if(m_mfp->heat_conductivity_pcs_name_vector.size()<2)
         {pcs_name2 = "TEMPERATURE1";}
         else
      {
         pcs_name2 = m_mfp->heat_conductivity_pcs_name_vector[1];
      }
      break;
      case 'S': mfp_id = 3;                       //SPECIFIC HEAT CAPACITY
      break;
      default:  mfp_id = -1;
      pcs_name1 = "PRESSURE1";
      pcs_name2 = "TEMPERATURE1";
   }
   //......................................................................

   int restore_mode=m_mfp->mode;
   m_mfp->mode = 0;
   m_mfp->node = node;

   /*  if (phase_number==0)
     {
        tp = PCSGet("PRESSURE1",true);           //NB 4.8.01
        val_idx=tp->GetNodeValueIndex("PRESSURE1"); // NB
        arguments[0] = tp->GetNodeValue(node,val_idx);
     } else if (phase_number==1)
     {
        tp = PCSGet("PRESSURE2",true);           //NB 4.8.01
        val_idx=tp->GetNodeValueIndex("PRESSURE2"); // NB
        arguments[0] = tp->GetNodeValue(node,val_idx);
     }

   tp = PCSGet("TEMPERATURE1",true);
   arguments[1] = tp->GetNodeValue(node,0);
   */
   tp = PCSGet(pcs_name1,true);                   //NB 4.8.01
   val_idx=tp->GetNodeValueIndex(pcs_name1);      // NB
   arguments[0] = tp->GetNodeValue(node,val_idx);

   tp = PCSGet(pcs_name2,true);                   //NB 4.8.01
   val_idx=tp->GetNodeValueIndex(pcs_name2);      // NB
   arguments[1] = tp->GetNodeValue(node,val_idx);

   //......................................................................
   switch(mfp_id)
   {
      case 0: mfp_value = m_mfp->Viscosity(arguments); break;
                                                  //NB 4.8.01
      case 1: mfp_value = m_mfp->Density(arguments); break;
      case 2: mfp_value = m_mfp->HeatConductivity(arguments); break;
                                                  //NB AUG 2009
      case 3: mfp_value = m_mfp->SpecificHeatCapacity(arguments); break;
      default: cout << "MFPGetNodeValue: no MFP data" << endl;
   }
   //......................................................................
   m_mfp->mode = restore_mode;                    //NB changeback
   return mfp_value;
}


/**************************************************************************
Task: derivative of density with respect to pressure
Programing: 09/2009 NB
 **************************************************************************/
double CFluidProperties::drhodP(double P, double T)
{
   double arguments[2];
   double rho1,rho2,drhodP;

   if (P<0) return 0;

   switch(compressibility_model_pressure)
   {
      case 0 :                                    // incompressible
         drhodP = 0;
         break;

      case 1 :                                    // constant slope compressibility
         drhodP = compressibility_pressure;
         break;

      case 2 :                                    // use of fct file
         drhodP = 0;                              // to be done
         break;
      case 3 :                                    // use of difference quotient
         arguments[1] = T;
                                                  // in case 3, compressibility_pressure acts as delta P
         arguments[0] = P+(compressibility_pressure/2.);
         rho1 = Density(arguments);

         arguments[0] = P-(compressibility_pressure/2.);
         rho2 = Density(arguments);

         //drhodP = (rho(p+dP/2)-rho(P-dP/2))/dP
                                                  // in case 3, compressibility_pressure acts as delta P
         drhodP = (rho1-rho2)/compressibility_pressure;

         break;                                   // use of difference quotient

      case 4 :                                    // use of analytical derivation
         drhodP = 0;                              // to be done
         break;
      default :
         drhodP = 0;

   }
   return drhodP;
}


/**************************************************************************
Task: derivative of density with respect to temperature
Programing: 09/2009 NB
 **************************************************************************/
double CFluidProperties::drhodT(double P, double T)
{
   double arguments[2];
   double rho1,rho2,drhodT;

   if (P<0) return 0;

   switch(compressibility_model_temperature)
   {
      case 0 :                                    // fluid is incompressible
         drhodT = 0;
         break;
      case 1 :                                    // constant slope compressibility, (for test cases)
         drhodT = compressibility_temperature;
         break;

      case 2 :                                    // use of fct file
         drhodT = 0;                              // to be done
         break;

      case 3 :                                    // use of difference quotient
         // in case 3, compressibility_temperature acts as delta T
         arguments[0] = P;

         arguments[1] = T+(compressibility_temperature/2.);
         rho1 = Density(arguments);

         arguments[1] = T-(compressibility_temperature/2.);
         rho2 = Density(arguments);

         drhodT = (rho1-rho2)/compressibility_temperature;
         break;

      case 4 :                                    // use of analytical derivation
         drhodT = 0;                              // to be done
         break;

      default :
         drhodT = 0;
         break;
   }
   return drhodT;
}


/**************************************************************************
Task: return super compressibility factor of the mixture
Programing:
05/2010 AKS
 **************************************************************************/
double CFluidProperties::CalCopressibility(long idx_elem, double p,double T)
{
   std::vector<double> roots;
   double a,b,A,B, R=8314.41;                     //KR w,Pc,dff,TG,Tc,PG,zn,z,ff
   double z1,z2,z3,h;                             //KR d
   a=MixtureSubProperity(0, idx_elem, p, T);
   b=MixtureSubProperity(1, idx_elem, p, T);
   A=a*p/(R*R*T*T);
   B=b*p/(R*T);
   z1=-(1-B);
   z2=(A-3*pow(B,2)-2*B);
   z3=-(A*B-pow(B,2)-pow(B,3));
   NsPol3(z1,z2,z3,&roots);                       //derives the roots of the polynomial
   //if(p < Pc)
   h=FindMax(roots);
   //if(p > Pc)
   //h=FindMin(roots);
   return h;
}


/**************************************************************************
Task: returns various parameters account molecular interation of the different components of mixture
Programing: 05/2010 AKS
 **************************************************************************/
double CFluidProperties::MixtureSubProperity(int properties, long idx_elem, double p, double T)
{
   CRFProcess* m_pcs;
   double ax, dens_arg[3];
   double mass_fraction[10], components_properties[10], w[10],pc[10],tc[10],fact[10];
   double variables = 0.0, R=8314.41;
   int component_number = (int) this->component_vector.size();
   dens_arg[0] = p;
   dens_arg[1] = T;
   dens_arg[2] = idx_elem;
   switch(properties)
   {
      case 0 :                                    // attraction parameter 'a'

         for (int i=0 ; i < component_number ; i++)
         {
            m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[i]->compname);
            mass_fraction[i] = this->component_vector[i]->CalcElementMeanConcNew( idx_elem, m_pcs );
            w[i]=this->component_vector[i]->acentric_factor;
            pc[i]=this->component_vector[i]->critical_pressure;
            tc[i]=this->component_vector[i]->critical_teperature;
            fact[i]=(1+(0.37464+1.54226*w[i]-0.2699*w[i]*w[i])*(1-pow(T/tc[i],0.5)));
            components_properties[i] = 0.45724*R*R*tc[i]*tc[i]*fact[i]*fact[i]/pc[i];
            for (int j=0 ; j < component_number ; j++)
            {
               m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[j]->compname);
               mass_fraction[j] = this->component_vector[j]->CalcElementMeanConcNew( idx_elem, m_pcs );
               w[j]=this->component_vector[j]->acentric_factor;
               pc[j]=this->component_vector[j]->critical_pressure;
               tc[j]=this->component_vector[j]->critical_teperature;
               fact[j]=(1+(0.37464+1.54226*w[j]-0.2699*w[j]*w[j])*(1-pow(T/tc[j],0.5)));
               components_properties[j] = 0.45724*R*R*tc[j]*tc[j]*fact[j]*fact[j]/pc[j];

               variables += mass_fraction[i]*mass_fraction[j]*pow(components_properties[i]*components_properties[j],0.5) ;
            }
         }
         break;

      case 1 :                                    // repulsion parameter 'b'
         for (int i=0 ; i < component_number ; i++)
         {
            m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[i]->compname);
            mass_fraction[i] = this->component_vector[i]->CalcElementMeanConcNew( idx_elem, m_pcs );
            pc[i]=this->component_vector[i]->critical_pressure;
            tc[i]=this->component_vector[i]->critical_teperature;
            components_properties[i] =0.07780*R*tc[i]/pc[i];

            variables += mass_fraction[i]*components_properties[i];
         }
         break;

      case 2 :                                    // potential parameter 'sigma'
         for (int i=0 ; i < component_number ; i++)
         {
            m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[i]->compname);
            if (m_pcs)
               mass_fraction[i] = this->component_vector[i]->CalcElementMeanConcNew( idx_elem, m_pcs );
            components_properties[i] =pow(pow(0.809, 3.0)*this->component_vector[i]->critical_volume,0.5);

            for (int j=0 ; j < component_number ; j++)
            {
               m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[j]->compname);
               mass_fraction[j] = this->component_vector[j]->CalcElementMeanConcNew( idx_elem, m_pcs );
               components_properties[j] = pow(pow(0.809, 3.0)*this->component_vector[j]->critical_volume,0.5);

               variables += mass_fraction[i]*mass_fraction[j]*components_properties[i]*components_properties[j];
            }
         }
         variables = pow(variables, 0.3333);
         break;

      case 3 :                                    // energy parameter 'epsilon'
         for (int i=0 ; i < component_number ; i++)
         {
            m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[i]->compname);
            if (m_pcs)
               mass_fraction[i] = this->component_vector[i]->CalcElementMeanConcNew( idx_elem, m_pcs );
            components_properties[i] = pow((this->component_vector[i]->critical_teperature/1.2593)*pow(0.809, 3.0)*this->component_vector[i]->critical_volume, 0.5);

            for (int j=0 ; j < component_number ; j++)
            {
               m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[j]->compname);
               mass_fraction[j] = this->component_vector[j]->CalcElementMeanConcNew( idx_elem, m_pcs );
               components_properties[j] =pow((this->component_vector[j]->critical_teperature/1.2593)*pow(0.809, 3.0)*this->component_vector[j]->critical_volume, 0.5);

               variables += mass_fraction[i]*mass_fraction[j]*components_properties[i]*components_properties[j];
            }
         }
         ax=pow(MixtureSubProperity(2, idx_elem, p, T), 3.0);
         variables /= ax;
         break;

      case 4 :                                    // acentric factor 'w'
         for (int i=0 ; i < component_number ; i++)
         {
            m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[i]->compname);
            mass_fraction[i] = this->component_vector[i]->CalcElementMeanConcNew( idx_elem, m_pcs );
            components_properties[i] =this->component_vector[i]->acentric_factor;
            variables += mass_fraction[i]*components_properties[i];
         }
         break;

      case 5 :                                    // molecular weight 'M'
         for (int i=0 ; i < component_number ; i++)
         {
            m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[i]->compname);
            mass_fraction[i] = this->component_vector[i]->CalcElementMeanConcNew( idx_elem, m_pcs );
            components_properties[i] =this->component_vector[i]->mol_mass;
            variables += mass_fraction[i]*components_properties[i];
         }
         break;

      case 6 :                                    // thermal conductivity 'k'
         for (int i=0 ; i < component_number ; i++)
         {
            m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[i]->compname);
            mass_fraction[i] = this->component_vector[i]->CalcElementMeanConcNew( idx_elem, m_pcs );
            components_properties[i] =  Fluid_Heat_Conductivity(Density(dens_arg)*mass_fraction[i], T, 2*i);
            variables += mass_fraction[i]*components_properties[i];
         }
         break;

      case 7 :                                    // heat capacity 'cp'
         for (int i=0 ; i < component_number ; i++)
         {
            m_pcs = PCSGetNew("MASS_TRANSPORT", this->component_vector[i]->compname);
            mass_fraction[i] = this->component_vector[i]->CalcElementMeanConcNew( idx_elem, m_pcs );
            components_properties[i] =  this->component_vector[i]->comp_capacity;
            variables += mass_fraction[i]*components_properties[i];
         }
         break;

      default :
         variables = 0;
   }
   return variables ;
}

