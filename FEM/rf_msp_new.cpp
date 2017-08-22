/**************************************************************************
FEMLib - Object: MSP Solid Properties
Task:
Programing:
08/2004 WW Implementation
last modified:
**************************************************************************/

// C++ STL
//#include <math.h>
//#include <string>
//#include <fstream>
//#include <iostream>
//#include <sstream>
#include <cfloat>

// FEM-Makros
#include "makros.h"
#include "rf_pcs.h"

// Time
#include "rf_tim_new.h"
#include "rf_msp_new.h"
#include "fem_ele_vec.h"
#include "fem_ele_std.h"
//#include "rf_mmp_new.h"
#include "pcs_dm.h"

#include "files0.h"                               // GetLineFromFile1
#include "tools.h"                                // GetLineFromFile

using namespace std;

vector<SolidProp::CSolidProperties*> msp_vector;
vector<string> msp_key_word_vector;               //OK

using FiniteElement::ElementValue_DM;
namespace SolidProp
{

   /**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   08/2004 OK Implementation for fuild properties
   08/2004 WW Modification for solid properties
   12/2005 WW Creep properties
   **************************************************************************/
   std::ios::pos_type CSolidProperties::Read(std::ifstream *msp_file)
   {
      char buffer[MAX_ZEILE];
      std::string sub_line;
      std::string line_string;
      std::string delimiter(" ");
      bool new_keyword = false;
      std::string hash("#");
      std::ios::pos_type position;
      //  ios::pos_type position0;
      std::string sub_string;
      bool new_subkeyword = false;
      std::string dollar("$");
      std::string delimiter_type(":");
      std::stringstream in_sd;

      int i=0, Size = 0;

      //========================================================================
      // Schleife ueber alle Phasen bzw. Komponenten
      while (!new_keyword)
      {
         new_subkeyword = false;

         position = msp_file->tellg();
         if(!GetLineFromFile(buffer,msp_file)) break;
         line_string = buffer;
         if(line_string.find(hash)!=string::npos)
         {
            new_keyword = true;
            break;
         }
         //....................................................................
         //NAME //OK
                                                  //subkeyword found
         if(line_string.find("$NAME")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd >> name;                        //sub_line
            in_sd.clear();
            continue;
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$SWELLING_PRESSURE_TYPE")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>>SwellingPressureType;
            if(SwellingPressureType==1||SwellingPressureType==2)
            {
               in_sd>>Max_SwellingPressure;
               in_sd.clear();
            }
                                                  //10.03.2008 WW
            else if (SwellingPressureType==3||SwellingPressureType==4)
            {
               if(!PCSGet("MULTI_PHASE_FLOW")||!PCSGet("RICHARDS_FLOW"))
               {
                  data_Youngs = new Matrix(9);
                  //! 0: \f$ \kappa_i0     \f$
                  //! 1: \f$ \alpha_i     \f$
                  //! 2: \f$ \kappa_{s0}  \f$
                  //! 3: \f$ \alpha_{sp}  \f$
                  //! 4: \f$ \alpha_{ss}  \f$
                  //! 5: \f$ p_ref        \f$
                  //! 6: \f$ buffer       \f$
                  if (SwellingPressureType==3)
                     in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >> (*data_Youngs)(2)
                        >> (*data_Youngs)(3) >> (*data_Youngs)(4)>> (*data_Youngs)(5);
                  else if (SwellingPressureType==4)
                     in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >> (*data_Youngs)(2);
                  in_sd.clear();
               }
               else
               {
                  std::cout<<"No multi-phase flow coupled. The thermal elatic model can only be used in H2 coupled proccess."<< std::endl;
                  std::cout<<"Quit the simulation now!"<< std::endl;
                  abort();
               }
            }
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$DENSITY")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>> Density_mode;
            if(Density_mode==0)                   // rho = f(x)
            {
               in_sd >> Size;
               in_sd.clear();
               data_Density = new Matrix(Size, 2);
               for(i=0; i<Size; i++)
               {
                  in_sd.str(GetLineFromFile1(msp_file));
                  in_sd >>(*data_Density)(i,0)>>(*data_Density)(i,1);
                  in_sd.clear();
               }
            }
            else if(Density_mode==1)              // rho = const
            {
               data_Density = new Matrix(1);
               in_sd >> (*data_Density)(0);
               in_sd.clear();
            }
         }
         //....................................................................
         if(line_string.find("$THERMAL")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>>sub_line>>ThermalExpansion;
            in_sd.clear();
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("CAPACITY")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd >> Capacity_mode;
            switch(Capacity_mode)
            {
               case 0:                            //  = f(x)
                  in_sd >> Size;
                  in_sd.clear();
                  data_Capacity = new Matrix(Size, 2);
                  for(i=0; i<Size; i++)
                  {
                     in_sd.str(GetLineFromFile1(msp_file));
                     in_sd >>(*data_Capacity)(i,0)>>(*data_Capacity)(i,1);
                     in_sd.clear();
                  }
                  break;
               case 1:                            //  = const
                  data_Capacity = new Matrix(1);
                  in_sd>> (*data_Capacity)(0);
                  in_sd.clear();
                  break;
               case 2:                            // boiling model for rock. WW
                  // 0. Wet capacity
                  // 1. Dry capacity
                  // 2. Boiling temperature
                  // 3. Boiling temperature range
                  // 4. Latent of vaporization
                  data_Capacity = new Matrix(5);
                  for(i=0; i<5; i++)
                     in_sd>> (*data_Capacity)(i);
                  in_sd.clear();
                  break;
               case 3:                            // DECOVALEX THM1, Bentonite
                  in_sd.clear();
                  break;
            }
         }

         //....................................................................
                                                  // subkeyword found
         if(line_string.find("CONDUCTIVITY")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>> Conductivity_mode;
            switch(Conductivity_mode)
            {
               case 0:                            //  = f(T) //21.12.2009 WW
                  in_sd>> Size;
                  in_sd.clear();
                  data_Conductivity = new Matrix(Size, 2);
                  for(i=0; i<Size; i++)
                  {
                     in_sd.str(GetLineFromFile1(msp_file));
                     in_sd >>(*data_Conductivity)(i,0)>>(*data_Conductivity)(i,1);
                     in_sd.clear();
                  }
                                                  //WW
                  conductivity_pcs_name_vector.push_back("TEMPERATURE1");
                  break;
               case 1:                            //  = const
                  data_Conductivity = new Matrix(1);
                  in_sd >> (*data_Conductivity)(0);
                  in_sd.clear();
                  break;
               case 2:                            // boiling model for rock. WW
                  // 0. Wet conductivity
                  // 1. Dry conductivity
                  // 2. Boiling temperature
                  // 3. Boiling temperature range
                  data_Conductivity = new Matrix(4);
                  for(i=0; i<4; i++)
                     in_sd>> (*data_Conductivity)(i);
                  in_sd.clear();
                  capacity_pcs_name_vector.push_back("TEMPERATURE1");
                  capacity_pcs_name_vector.push_back("SATURATION1");
                  break;
               case 3:                            // DECOVALEX THM1, Bentonite
                  in_sd.clear();
                  capacity_pcs_name_vector.push_back("TEMPERATURE1");
                  capacity_pcs_name_vector.push_back("SATURATION1");
                  break;
               case 4:                            //  = f(S) //21.12.2009 WW
                  in_sd>> Size;
                  in_sd.clear();
                  data_Conductivity = new Matrix(Size, 2);
                  for(i=0; i<Size; i++)
                  {
                     in_sd.str(GetLineFromFile1(msp_file));
                     in_sd >>(*data_Conductivity)(i,0)>>(*data_Conductivity)(i,1);
                     in_sd.clear();
                  }
                  break;
            }
         }
         /*
             //....................................................................
             if(line_string.find("$THERMAL_CONDUCTIVITY_TENSOR")!=string::npos) { //subkeyword found
              *msp_file >> thermal_conductivity_tensor_type_name;
               switch(thermal_conductivity_tensor_type_name[0]) {
                 case 'I': // isotropic
                   thermal_conductivity_tensor_type = 0;
                  *msp_file >> thermal_conductivity_tensor[0];
                   break;
                 case 'O': // orthotropic
                   break;
         case 'A': // anisotropic
         break;
         default:
         cout << "Error in CSolidProperties::Read: no valid thermal conductivity tensor type" << endl;
         break;
         }
         msp_file->ignore(MAX_ZEILE,'\n');
         continue;
         }
         */
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$ELASTICITY")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>>sub_line>>PoissonRatio;
            in_sd.clear();
         }
         //....................................................................
                                                  //12.2009. WW
         if(line_string.find("$EXCAVATION")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>>excavation;
            in_sd.clear();
         }
         //....................................................................
         if(line_string.find("YOUNGS_MODULUS")!=string::npos)
         {                                        // subkeyword found
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd >> Youngs_mode;
            int type = Youngs_mode;
            if(type>9&&type<14)  type = 1000;
            switch(type)                          // 15.03.2008 WW
            {
               case 0:                            //  = f(x)
                  in_sd >> Size;
                  in_sd.clear();
                  data_Youngs = new Matrix(Size, 2);
                  for(i=0; i<Size; i++)
                  {
                     in_sd.str(GetLineFromFile1(msp_file));
                     in_sd>>(*data_Youngs)(i,0)>>(*data_Youngs)(i,1);
                     in_sd.clear();
                  }
                  break;
               case 1:                            //  = const
                  data_Youngs = new Matrix(1);
                  in_sd >> (*data_Youngs)(0);
                  in_sd.clear();
                  break;                          // UJG 24.11.2009
               case 2:                            //  = const
                  // data_Youngs Lubby1 model
                  //  0: E_0
                  //  1: a (factor)
                  //  2: n (exponent)
                  data_Youngs = new Matrix(3);
                  in_sd >> (*data_Youngs)(0)>>(*data_Youngs)(1)>>(*data_Youngs)(2);
                  in_sd.clear();
                  break;
               case 1000:                         // case 10-13: transverse isotropic linear elasticity (UJG 24.11.2009)
                  // data_Youngs transverse isotropic linear elasticity
                  //  0: E_i     (Young's modulus of the plane of isotropy)
                  //  1: E_a     (Young's modulus w.r.t. the anisotropy direction)
                  //  2: nu_{ia} (Poisson's ratio w.r.t. the anisotropy direction)
                  //  3: G_a     (shear modulus w.r.t. the anisotropy direction)
                  //  4: n_x     (x-coefficient of the local axis of anisotropy (2D case: -\sin\phi))
                  //  5: n_y     (y-coefficient of the local axis of anisotropy (2D case: \cos\phi))
                  //  6: n_z     (z-coefficient of the local axis of anisotropy (2D case: 0))
                  data_Youngs = new Matrix(7);
                  in_sd >> (*data_Youngs)(0)>>(*data_Youngs)(1)>>(*data_Youngs)(2)>>(*data_Youngs)(3)
                     >>(*data_Youngs)(4)>>(*data_Youngs)(5)>>(*data_Youngs)(6);
                  in_sd.clear();
                  break;
#ifdef RFW_FRACTURE
               case 2:                            // = f(aperture), for fracture material groups, RFW 04/2005,  RFW 09/12/2005
                  data_Youngs = new Matrix(3);
                  in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >> (*data_Youngs)(2);
                  in_sd.clear();
                  break;
#endif
            }
         }
         //....................................................................
         if(line_string.find("$CREEP")!=string::npos)
         {
            if(line_string.find("NORTON")!=string::npos)
            {
               Creep_mode=1;
               /*! \subsection Norton creep model */
               /*! \f$\dot\epsilon_s=A \left(\frac{\sigma_v}{\sigma^\ast}\right)^n\f$ */
               // data_Creep:
               //  0: A,   coefficient
               //  1: n,   exponential
               data_Creep = new Matrix(2);
               in_sd.str(GetLineFromFile1(msp_file));
               in_sd>>(*data_Creep)(0);
               in_sd>>(*data_Creep)(1);
               in_sd.clear();
            }
            if(line_string.find("BGRA")!=string::npos)
            {
               Creep_mode=2;
               /*! \subsection Temperature dependent creep model by BGR */
               /*! \f$\dot\epsilon_s=A\exp^{-Q/RT}\left(\frac{\sigma_v}{\sigma^\ast}\right)^n\f$ */
               // data_Creep:
               //  0: A,   coefficient
               //  1: n,   exponential
               //  2: Q,   activation energy
               data_Creep = new Matrix(3);
               in_sd.str(GetLineFromFile1(msp_file));
               in_sd>>(*data_Creep)(0);
               in_sd>>(*data_Creep)(1);
               in_sd>>(*data_Creep)(2);
               in_sd.clear();
            }
            //....................................................................
            if(line_string.find("LUBBY2")!=string::npos)
            {
               Creep_mode=1000;
               // data_Creep:
               //  0: eta_m
               //  1: m
               //  2: l
               //  3: eta_k
               //  4: k1
               //  5: k2
               //  6: G_k
               data_Creep = new Matrix(7,2);
               in_sd.str(GetLineFromFile1(msp_file));
               for(i=0; i<7; i++)
                  in_sd>>(*data_Creep)(i,0);
               in_sd.clear();
            }
         }
         //....................................................................
         if(line_string.find("$BIOT_CONSTANT")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>>biot_const;
            in_sd.clear();
         }
         //....................................................................
         if(line_string.find("$STRESS_INTEGRATION_TOLERANCE")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>>f_tol>>s_tol;
            in_sd.clear();
         }
         if(line_string.find("$STRESS_UNIT")!=string::npos)
         {

            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>>sub_line;
            if(sub_line.compare("MegaPascal")==0)
               grav_const = 9.81e-6;
            else if(sub_line.find("KiloPascal")==0)
               grav_const = 9.81e-3;
            in_sd.clear();
         }
         //....................................................................
                                                  // subkeyword found
         if(line_string.find("$PLASTICITY")!=string::npos)
         {
            in_sd.str(GetLineFromFile1(msp_file));
            in_sd>>sub_line;
            in_sd.clear();
            if(sub_line.find("DRUCKER-PRAGER")!=string::npos)
            {
               devS = new double[6];
               Plasticity_type=1;
                                                  // No return mapping
               if(sub_line.find("NORETURNMAPPING")!=string::npos)
               {
                  Plasticity_type=10;
                  dFds = new double[6];
                  dGds = new double[6];
                  D_dFds = new double[6];
                  D_dGds = new double[6];
               }
               Size = 5;
               /*
               Material parameters for Cam-Clay model
               i : parameter
               0 : The initial yield stress
               1 : Plastic hardening parameter
               2 : Internal frection angle
               3 : Dilatancy angle
               4 : Localized softening modulus
               */
            }
            else if(sub_line.find("SINGLE_YIELD_SURFACE")!=string::npos)
            {
               Plasticity_type=2;
               Size = 23;
               AllocateMemoryforSYS();
               /*
               Material parameters for Single yield surface model
                 i: parameter
                 0: alpha0
                 1: beta0
                 2: delta0
                 3: epsilon0
                 4: kappa0
                 5: gamma0
                 6: m0

               7: alpha1
               8: beta1
               9: delta1
               10: epsilon1
               11: kappa1
               12: gamma1
               13: m1

               14: Psi1
               15: Psi2

               16: Ch
               17: Cd
               18: br
               19: mr

               20: Initial stress_xx
               21: Initial stress_yy
               22: Initial stress_zz
               */
            }
            else if(sub_line.find("CAM-CLAY")!=string::npos)
            {
               Plasticity_type=3;
               Size = 10;
               /*
               Material parameters for Cam-Clay model
                i: parameter
                0 : M: slope of the critical line
                1 : Lamda, the virgin compression index
                2 : Kappa, swelling index
                3 : p0, preconsolidation pressure
                4 : e0, initial void ratio
                5 : OCR
                6 : Initial stress_xx
                7 : Initial stress_yy
               8 : Initial stress_zz
               9 : Mimimum p: ( stress_xx + stress_yy + stress_zz )/3
               */

            }
            data_Plasticity = new Matrix(Size);
            for(i=0; i<Size; i++)
            {
               in_sd.str(GetLineFromFile1(msp_file));
               in_sd >>(*data_Plasticity)(i);
               in_sd.clear();
            }
         }
         in_sd.clear();
      }
      return position;
   }

   //==========================================================================

   /**************************************************************************
   FEMLib-Method:
   Task: Constructor and destructor
   Programing:
   08/2004 WW Implementation
   04/2005 WW Curve dependency
   **************************************************************************/
   CSolidProperties::CSolidProperties()
      : data_Youngs(NULL), data_Density(NULL), data_Capacity(NULL),
      data_Conductivity(NULL), data_Plasticity(NULL), data_Creep(NULL)
   {
      PoissonRatio = 0.2;
      ThermalExpansion = 0.0;
      biot_const = 1.0;
      // Default, all data are constant
      Density_mode = -1;
      Youngs_mode = -1;
      Capacity_mode = -1;
      Conductivity_mode = -1;
      Creep_mode = -1;
      grav_const = 9.81;                          //WW
      excavation = -1;                            //12.2009. WW
      excavated = false;                          //To be .....  12.2009. WW

      SwellingPressureType = -1;
      Max_SwellingPressure = 0.0;
      // Default, elasticity
      Plasticity_type = -1;

      E = Lambda = G = K = 0.0;
      devS = NULL;
      axisymmetry = false;
      dl2 = 0.0;
      // SYS
      d2G_dSdS=NULL;
      d2G_dSdM=NULL;
      LocalJacobi=NULL;
      inv_Jac=NULL;
      sumA_Matrix=NULL;
      rhs_l=NULL;
      x_l=NULL;
      Li=NULL;
      // Drucker-Prager
      dFds = NULL;
      dGds = NULL;
      D_dFds = NULL;
      D_dGds = NULL;
      // Curve variable type
      // 0: Time
      // 1: ...
      CurveVariableType_Conductivity=-1;
      mode = 0;                                   // Gauss point values //OK
      //
      s_tol = 1e-9;
      f_tol = 1e-9;

      Crotm = NULL;                               // rotation matrix for matrices: UJG 25.11.2009
      D_tran = NULL;                              //UJG/WW
   }
   CSolidProperties::~CSolidProperties()
   {
      if(data_Density) delete data_Density;
      if(data_Density) delete data_Youngs;
      if(data_Plasticity) delete data_Plasticity;
      if(data_Capacity) delete data_Capacity;
      if(data_Conductivity) delete data_Conductivity;
      if(data_Creep) delete data_Creep;
      if(devS) delete [] devS;

      if(Crotm) delete Crotm;                     // rotation matrix for matrices: UJG 25.11.2009
      if(D_tran) delete D_tran;                   // rotation matrix for matrices: UJG 25.11.2009
      data_Density=NULL;
      data_Youngs=NULL;
      data_Plasticity=NULL;
      data_Capacity=NULL;
      data_Conductivity=NULL;
      data_Creep = NULL;
      devS=NULL;
      Crotm = NULL;
      D_tran = NULL;

      if(d2G_dSdS) delete d2G_dSdS;
      if(d2G_dSdM) delete d2G_dSdM;
      if(LocalJacobi) delete LocalJacobi;         // To store local Jacobi matrix
      if(inv_Jac) delete inv_Jac;                 // To store the inverse of the  Jacobi matrix
      if(sumA_Matrix) delete sumA_Matrix;
      if(rhs_l) delete [] rhs_l;                  // To store local unknowns of 15
      if(x_l) delete [] x_l;                      // To store local unknowns of 15
      if(Li) delete [] Li;

      if(dFds) delete [] dFds;
      if(dGds) delete [] dGds;
      if(D_dFds) delete [] D_dFds;
      if(D_dGds) delete [] D_dGds;
      dFds = NULL;
      dGds = NULL;
      D_dFds = NULL;
      D_dGds = NULL;
      d2G_dSdS=NULL;
      d2G_dSdM=NULL;
      LocalJacobi=NULL;
      inv_Jac=NULL;
      sumA_Matrix=NULL;
      rhs_l=NULL;
      x_l=NULL;
      Li=NULL;
   }
   //----------------------------------------------------------------------------

   /**************************************************************************
   FEMLib-Method: CSolidProperties::CalulateValue
   Task: Linear interpolation
   Programing:
   08/2004 WW Implementation
   **************************************************************************/
   double CSolidProperties::CalulateValue
      (const Matrix *data,  const double x) const
   {
      int i;
      double Val = 0.0;

      if(!data) return 0.0;

      const int Size = data->Rows();
      // If given varial is not in the range of data
      if(x<(*data)(0,0)) Val = (*data)(0,1);
      else if(x>=(*data)(Size-1,0)) Val = (*data)(Size-1,1);
      else
      {
         for(i=0; i<Size-1; i++)
         {
            if((x>=(*data)(i,0))&&(x<(*data)(i+1,0)))
               Val = (x-(*data)(i,0))*((*data)(i+1,1)-(*data)(i,1))
                  /((*data)(i+1,0)-(*data)(i,0))+(*data)(i,1);
         }
      }
      return Val;
   }

   /**************************************************************************
   FEMLib-Method: CSolidProperties::Density(double refence = 0.0)
   Task: Get density
   Programing:
   08/2004 WW Implementation
   **************************************************************************/
   double CSolidProperties::Density(double refence )
   {
      double val = 0.0;
      switch(Density_mode)
      {
         case 0:
            val = CalulateValue(data_Density, refence);
            break;
         case 1:
            val = (*data_Density)(0);
            break;
      }
      return val;
   }
   // Initialize density
   void CSolidProperties::NullDensity()
   {
      (*data_Density) = 0.0;
   }

   /**************************************************************************
   FEMLib-Method: CSolidProperties::Heat_Capacity(const double refence = 0.0) const
   Task: Get heat capacity
   Programing:
   08/2004 WW Implementation
   **************************************************************************/
   double CSolidProperties::Heat_Capacity(double refence)
   {
      double val = 0.0;
      switch(Capacity_mode)
      {
         case 0:
            val = CalulateValue(data_Capacity, refence);
            break;
         case 1:
            val = (*data_Capacity)(0);
            break;
         case 3:
            //WW        val=1.38*(273.15+refence)+732.5;
            val=1.38*refence+732.5;
            break;
         default:
            val = (*data_Capacity)(0);
            break;
      }
      return val;
   }

   /*************************************************************************
   FEMLib-Method:
   Task: Get heat phase change temperature
   Programing:
   09/2005 WW Implementation
   **************************************************************************/
   bool CSolidProperties:: CheckTemperature_in_PhaseChange
      (const double T0, const double T1)
   {
      bool stat = false;
      double T_a=(*data_Capacity)(2);
      double T_b=(*data_Capacity)(2)+(*data_Capacity)(3);
      switch(Capacity_mode)
      {
         case 0:
            break;
         case 1:
            break;
         case 2:
            if((T1>T0)&&(*data_Capacity)(4)>0.0)
            {
               if((T0<T_a)&&(T1>T_a))
                  stat = true;
               else if((T0<T_a)&&(T1>T_b))
                  stat = true;
               else if((T0<T_b)&&(T1>T_b))
                  stat = true;
               else if((T1>=T_a)&&(T1<=T_b))
                  stat = true;
            }
            break;
      }
      return stat;
   }

   /*************************************************************************
   FEMLib-Method:
   Task: Get heat capacity with boiling model
         latent_factor=density_w*porosity*saturation
   Programing:
   09/2005 WW Implementation
   **************************************************************************/
   double CSolidProperties::Heat_Capacity(double temperature, double porosity, double Sat)
   {
      double val = 0.0;
      double sign = 1;
      double dens = fabs(Density());
      // If sign =1, temperature increases in evolution. Otherwise, decrease.
      if(fabs(temperature)>1.e-9) sign = fabs(temperature)/temperature;
      temperature *= sign;
      CFluidProperties *m_mfp = NULL;
      m_mfp = mfp_vector[0];

      // 0. Wet capacity
      // 1. Dry capacity
      // 2. Boiling temperature
      // 3. Boiling temperature range
      // 4. Latent of vaporization
      if(temperature< (*data_Capacity)(2))        // Wet
         val =  dens*(*data_Capacity)(0);
      else if((temperature>=(*data_Capacity)(2))&&
         (temperature<((*data_Capacity)(2)+(*data_Capacity)(3))))
      {
         if(sign>0.0)                             // Temperature increase
            val =  dens*(*data_Capacity)(0)+dens*((*data_Capacity)(1)-(*data_Capacity)(0))
               *(temperature-(*data_Capacity)(2))/(*data_Capacity)(3)
               + porosity*Sat*m_mfp->Density()*(*data_Capacity)(4)/(*data_Capacity)(3) ;
         else
            val =  dens*(*data_Capacity)(0)+dens*((*data_Capacity)(1)-(*data_Capacity)(0))
               *(temperature-(*data_Capacity)(2))/(*data_Capacity)(3);
      }
      else
         val =  (*data_Capacity)(1);
      return val;
   }

   /**************************************************************************
   FEMLib-Method:
   Task: Get heat capacity with boiling model
         latent_factor=density_w*porosity*saturation
        temperature increases in evolution
   Programing:
   09/2005 WW Implementation
   **************************************************************************/
   double CSolidProperties::Enthalpy(double temperature, const double latent_factor )
   {
      double val = 0.0;
      double T0=0.0;                              //�C, a reference temperature
      double dT=0.0;

      // 0. Wet capacity
      // 1. Dry capacity
      // 2. Boiling temperature
      // 3. Boiling temperature range
      // 4. Latent of vaporization
      if(temperature< (*data_Capacity)(2))        // Wet
         val =  (*data_Capacity)(0)*(temperature-T0);
      else if((temperature>=(*data_Capacity)(2))&&
         (temperature<((*data_Capacity)(2)+(*data_Capacity)(3))))
      {
         dT = temperature-(*data_Capacity)(2)-T0;
         val =  (*data_Capacity)(0)*((*data_Capacity)(2)-T0);
         val += dT*latent_factor*(*data_Capacity)(4)/(*data_Capacity)(3)/Density();
         val +=  (*data_Capacity)(0)*dT
            +0.5*dT*dT*((*data_Capacity)(1)-(*data_Capacity)(0))/(*data_Capacity)(3);
      }
      else
      {
         dT =   (*data_Capacity)(3);
         val =  (*data_Capacity)(0)*((*data_Capacity)(2)-T0);
         val += dT*latent_factor*(*data_Capacity)(4)/(*data_Capacity)(3)/Density();
         val +=  (*data_Capacity)(0)*dT
            +0.5*dT*dT*((*data_Capacity)(1)-(*data_Capacity)(0))/(*data_Capacity)(3);
         val +=  (*data_Capacity)(1)*(temperature-(*data_Capacity)(3)-(*data_Capacity)(2)-T0);
      }
      return val*fabs(Density());
   }

   /**************************************************************************
   FEMLib-Method: CSolidProperties::Heat_Capacity(const double refence = 0.0) const
   Task: Get density
   Programing:
   08/2004 WW Implementation
   **************************************************************************/
   double CSolidProperties::Heat_Conductivity(double refence)
   {
      double val = 0.0;
      switch(Conductivity_mode)
      {
         case 0:
            val = CalulateValue(data_Conductivity, refence);
            break;
         case 1:
            val = (*data_Conductivity)(0);
            break;
         case 2:
            // 0. Wet conductivity
            // 1. Dry conductivity
            // 2. Boiling temperature
            // 3. Boiling temperature range
            if(refence< (*data_Conductivity)(2))  // Wet
               val =  (*data_Conductivity)(0);
            else if((refence>=(*data_Conductivity)(2))&&
               (refence<((*data_Conductivity)(2)+(*data_Conductivity)(3))))
               val =  (*data_Conductivity)(0)+((*data_Conductivity)(0)-(*data_Conductivity)(1))
                     *(refence-(*data_Conductivity)(2))/(*data_Conductivity)(3);
            else
               val =  (*data_Conductivity)(1);
            break;
         case 3:                                  // refence: saturation
            //val = 1.28-0.71/(1+10.0*exp(refence-0.65));  //MX
            val = 1.28-0.71/(1+exp(10.0*(refence-0.65)));
            break;
         case 4:                                  //21.12.2009. WW
            val = CalulateValue(data_Conductivity, refence);
            break;

      }
      return val;
   }

   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2004 OK MSP implementation
   03/2005 WW Conductivity from input file
   last modification:
   ToDo: geo_dimension
   **************************************************************************/
   void CSolidProperties::HeatConductivityTensor(const int dim, double* tensor, int group)
   {
      group = group;                              //OK411
      //static double tensor[9];
      double temperature = 0.0;
      double saturation = 0.0;
      int i = 0;
                                                  //WW
      CalPrimaryVariable(conductivity_pcs_name_vector);
      //--------------------------------------------------------------------
      // MMP medium properties
      //WW CMediumProperties *m_mmp = NULL;
      //WW m_mmp = mmp_vector[group];  //MX
      // Test for DECOVALEX
      // thermal_conductivity_tensor[0] =0.5+0.8*PCSGetELEValue(number,NULL,theta,"SATURATION1");

      //There are a number of cases where the heat conductivity tensor is defined by the capacity model;
      switch (Conductivity_mode)
      {
         case 0:
                                                  //WW
            thermal_conductivity_tensor[0] = Heat_Conductivity(primary_variable[0]);
            for(i=0; i<dim; i++)
               tensor[i*dim+i] =  thermal_conductivity_tensor[0];
            break;
         case 1:
                                                  //WW
            thermal_conductivity_tensor[0] = Heat_Conductivity(0);
            for(i=0; i<dim; i++)
               tensor[i*dim+i] =  thermal_conductivity_tensor[0];
            break;
         case 2:                                  // Boiling model. DECOVALEX THM2
            temperature = primary_variable[0];
            //for(i=0; i<dim*dim; i++) mat[i] = 0.0;
            for(i=0; i<dim; i++)
               tensor[i*dim+i] = Heat_Conductivity(temperature);
            break;
         case 3:                                  // DECOVALEX THM1
            saturation = primary_variable[1];
            //for(i=0; i<dim*dim; i++) mat[i] = 0.0;
            for(i=0; i<dim; i++)
               tensor[i*dim+i] = Heat_Conductivity(saturation);
            break;
         default:                                 //Normal case
                                                  //WW
            thermal_conductivity_tensor[0] = Heat_Conductivity();
            thermal_conductivity_tensor_type = 0;
            switch(dim)
            {
               case 1:                            // 1-D
                  tensor[0] = thermal_conductivity_tensor[0];
                  break;
               case 2:                            // 2-D
                  if(thermal_conductivity_tensor_type==0)
                  {
                     tensor[0] = thermal_conductivity_tensor[0];
                     tensor[1] = 0.0;
                     tensor[2] = 0.0;
                     tensor[3] = thermal_conductivity_tensor[0];
                  }
                  else if(thermal_conductivity_tensor_type==1)
                  {
                     tensor[0] = thermal_conductivity_tensor[0];
                     tensor[1] = 0.0;
                     tensor[2] = 0.0;
                     tensor[3] = thermal_conductivity_tensor[1];
                  }
                  else if(thermal_conductivity_tensor_type==2)
                  {
                     tensor[0] = thermal_conductivity_tensor[0];
                     tensor[1] = thermal_conductivity_tensor[1];
                     tensor[2] = thermal_conductivity_tensor[2];
                     tensor[3] = thermal_conductivity_tensor[3];
                  }
                  break;
               case 3:                            // 3-D
                  if(thermal_conductivity_tensor_type==0)
                  {
                     tensor[0] = thermal_conductivity_tensor[0];
                     tensor[1] = 0.0;
                     tensor[2] = 0.0;
                     tensor[3] = 0.0;
                     tensor[4] = thermal_conductivity_tensor[0];
                     tensor[5] = 0.0;
                     tensor[6] = 0.0;
                     tensor[7] = 0.0;
                     tensor[8] = thermal_conductivity_tensor[0];
                  }
                  else if(thermal_conductivity_tensor_type==1)
                  {
                     tensor[0] = thermal_conductivity_tensor[0];
                     tensor[1] = 0.0;
                     tensor[2] = 0.0;
                     tensor[3] = 0.0;
                     tensor[4] = thermal_conductivity_tensor[1];
                     tensor[5] = 0.0;
                     tensor[6] = 0.0;
                     tensor[7] = 0.0;
                     tensor[8] = thermal_conductivity_tensor[2];
                  }
                  else if(thermal_conductivity_tensor_type==2)
                  {
                     tensor[0] = thermal_conductivity_tensor[0];
                     tensor[1] = thermal_conductivity_tensor[1];
                     tensor[2] = thermal_conductivity_tensor[2];
                     tensor[3] = thermal_conductivity_tensor[3];
                     tensor[4] = thermal_conductivity_tensor[4];
                     tensor[5] = thermal_conductivity_tensor[5];
                     tensor[6] = thermal_conductivity_tensor[6];
                     tensor[7] = thermal_conductivity_tensor[7];
                     tensor[8] = thermal_conductivity_tensor[8];
                  }
                  break;
            }
      }
      //  return tensor;
   }

   /**************************************************************************
   FEMLib-Method: CSolidProperties::Youngs_Modulus(const double refence = 0.0) const
   Task: Get density
   Programing:
   08/2004 WW Implementation
   **************************************************************************/
#ifndef RFW_FRACTURE
   double CSolidProperties::Youngs_Modulus(double refence)
   {
      double val = 0.0;
      switch(Youngs_mode)
      {
         case 0:
            val = CalulateValue(data_Youngs, refence);
            break;
         case 1:
            val = (*data_Youngs)(0);
            break;
      }
      return val;
   }
#endif

#ifdef RFW_FRACTURE
                                                  //RFW changed
   double CSolidProperties::Youngs_Modulus(CElem* elem, double refence)
   {
      double val = 0.0;

      switch(Youngs_mode)
      {
         case 0:
            val = CalulateValue(data_Youngs, refence);
            break;
         case 1:
            val = (*data_Youngs)(0);
            break;
         case 2:
            double E1 = (*data_Youngs)(0);
            double E2 = (*data_Youngs)(1);
            double min_aperture = (*data_Youngs)(2);
            double step_size = min_aperture/10;

            CMediumProperties *m_mmp = NULL;
            m_mmp = mmp_vector[elem->GetPatchIndex()];
            double aperture = m_mmp->CalculateFracAperture(elem, step_size);

            if (  (aperture > (min_aperture+0.00000001))  ||  (aperture == -10) )
            {
               val=E1;
            }
            else
            {
               val = E2;
            }

            break;

      }
      return val;
   }
#endif

#ifdef RFW_FRACTURE
   /**************************************************************************
   FEMLib-Method: CSolidProperties::Youngs_Min_Aper(const double refence = 0.0) const
   Task: Get minimum aperture value for fracture deformation
   Programing:
   07/2005 RFW Implementation
   **************************************************************************/
                                                  //RFW 07/2005
   double CSolidProperties::Get_Youngs_Min_Aperture(CElem *elem)
   {

      int group = elem->GetPatchIndex();
      CSolidProperties* mat_pointer;
      Matrix* data;
      double min_aperture;

      mat_pointer = msp_vector[group];
      data = mat_pointer->data_Youngs;
      min_aperture = (*data)(2);
      return min_aperture;
   }
#endif

   //-------------------------------------------------------------------------
   // Kronecker function
   double CSolidProperties::Kronecker(const int ii, const int jj)
   {
      double delta=0.0;
      if(ii==jj) delta=1.0;
      return delta;
   }

   /**************************************************************************
   FEMLib-Method: CSolidProperties::Calculate_Lame_Constant()
   Task: Get density
   Programing:
   08/2004 WW Implementation
   **************************************************************************/
#ifndef RFW_FRACTURE
   void  CSolidProperties::Calculate_Lame_Constant()
   {
      double nv = Poisson_Ratio();
      E = Youngs_Modulus();                       // Constant at present
      Lambda = E * nv / ((1. + nv) * (1. - 2. * nv));
      G = 0.5 * E / (1. + nv);
      K=(3.0*Lambda+2.0*G)/3.0;
   }
#endif
#ifdef RFW_FRACTURE
   //with a CElem pointer, for the fracture deformation Youngs Modulus is element specific, not material specific
                                                  //RFW changed
   void  CSolidProperties::Calculate_Lame_Constant( CElem* elem)
   {
      E = Youngs_Modulus(elem);                   // RFW 06/2006
      Lambda = E * Poisson_Ratio() /
         ((1. + Poisson_Ratio()) * (1. - 2. * Poisson_Ratio()));
      G = 0.5 * E / (1. + Poisson_Ratio());
      K=(3.0*Lambda+2.0*G)/3.0;
   }
#endif
   /*************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::ElasticConsitutive

    Aufgabe:
      Fill the lastic constitutive

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
          const int Dimension  :  Dimension of the real space
          double *D_e          :  Elastic constitutive
    Ergebnis:
      - void -

   Programmaenderungen:
   11/2003   WW   Set plastic parameter

   *************************************************************************/
   void CSolidProperties::ElasticConsitutive(const int Dimension, Matrix *D_e) const
   {
      (*D_e) = 0.0;
      (*D_e)(0,0) = Lambda + 2 * G;
      (*D_e)(0,1) = Lambda;
      (*D_e)(0,2) = Lambda;

      (*D_e)(1,0) = Lambda;
      (*D_e)(1,1) = Lambda + 2 * G;
      (*D_e)(1,2) = Lambda;

      (*D_e)(2,0) = Lambda;
      (*D_e)(2,1) = Lambda;
      (*D_e)(2,2) = Lambda + 2 * G;

      (*D_e)(3,3) = G;
      //Plane stress
      // plane stress, only for test
      //(*D_e)(0,0) = (1.0-Mu)*Lambda + 2 * G;
      //(*D_e)(0,1) = Lambda;

      //(*D_e)(1,0) = Lambda;
      // (*D_e)(1,1) = (1.0-Mu)*Lambda + 2 * G;
      // (*D_e)(3,3) = G;

      if(Dimension==3)
      {
         (*D_e)(4,4) = G;
         (*D_e)(5,5) = G;
      }
   }
   /*************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::ElasticConstitutiveTransverseIsotropic

    Aufgabe:
      Generate material matrix in case of transverse isotropic linear elasticity

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
          const int Dimension  :  Dimension of the real space
          double *D_e          :  Elastic constitutive
    Ergebnis:
      - void -

   Programmaenderungen:
   24.11.2009  UJG  first implementation

   *************************************************************************/
   void CSolidProperties::ElasticConstitutiveTransverseIsotropic(const int Dimension)
   {
      double aii,aai,bii,bai,cii,cai;

      double ni  = Poisson_Ratio();
      double Ei  = (*data_Youngs)(0);
      double Ea  = (*data_Youngs)(1);
      double nia = (*data_Youngs)(2);
      double Ga  = (*data_Youngs)(3);

      double nai  = nia*(Ea/Ei);
      double Disk = 1.0+ni;

      Matrix *D_e = NULL;

      Disk *= (1.0-ni-2.0*(nia*nai));
      Disk /= (Ei*Ei*Ea);

      aii = (1.0-(nai*nia))/(Ei*Ea*Disk);
      aai = (1.0-(ni*ni))/(Ei*Ei*Disk);
      bii = (ni+(nia*nai))/(Ei*Ea*Disk);
      bai = (nia*(1.0+ni))/(Ei*Ei*Disk);
      cii = Ei/(2.0*(1+ni));
      cai = Ga;

      int size;
      size = Dimension*2;

      switch(Youngs_mode)
      {
         case 10:
            int i,j,l,m;
            D_e = new Matrix(size,size);
            (*D_e) = 0.0;

            if(Dimension==3)
            {
               (*D_e)(0,0) = aii;
               (*D_e)(0,1) = bii;
               (*D_e)(0,2) = bai;

               (*D_e)(1,0) = bii;
               (*D_e)(1,1) = aii;
               (*D_e)(1,2) = bai;

               (*D_e)(2,0) = bai;
               (*D_e)(2,1) = bai;
               (*D_e)(2,2) = aai;

               (*D_e)(3,3) = cii;
               (*D_e)(4,4) = cai;
               (*D_e)(5,5) = cai;
            }
            else
            {
               (*D_e)(0,0) = aii;
               (*D_e)(0,1) = bai;
               (*D_e)(0,2) = bii;

               (*D_e)(1,0) = bai;
               (*D_e)(1,1) = aai;
               (*D_e)(1,2) = bai;

               (*D_e)(2,0) = bii;
               (*D_e)(2,1) = bai;
               (*D_e)(2,2) = aii;

               (*D_e)(3,3) = cai;
            }

            D_tran = new Matrix(size, size);
            (*D_tran) = 0.;
            for(i=0; i<size; i++)
            {
               for(j=0; j<size; j++)
               {
                  for(l=0; l<size; l++)
                  {
                     for(m=0; m<size; m++)
                        (*D_tran)(i,j) += (*Crotm)(l,i)*(*D_e)(l,m)*(*Crotm)(m,j);
                  }
               }
            }
            delete D_e;
            D_e = NULL;
            break;
         case 11:
            D_tran = new Matrix(size, size);
            (*D_tran)(0,0) = aai;
            (*D_tran)(0,1) = bai;
            (*D_tran)(0,2) = bai;

            (*D_tran)(1,0) = bai;
            (*D_tran)(1,1) = aii;
            (*D_tran)(1,2) = bii;

            (*D_tran)(2,0) = bai;
            (*D_tran)(2,1) = bii;
            (*D_tran)(2,2) = aii;

            (*D_tran)(3,3) = cai;

            if(Dimension==3)
            {
               (*D_tran)(4,4) = cai;
               (*D_tran)(5,5) = cii;
            }
            break;
         case 12:
            D_tran = new Matrix(size, size);
            (*D_tran)(0,0) = aii;
            (*D_tran)(0,1) = bai;
            (*D_tran)(0,2) = bii;

            (*D_tran)(1,0) = bai;
            (*D_tran)(1,1) = aai;
            (*D_tran)(1,2) = bai;

            (*D_tran)(2,0) = bii;
            (*D_tran)(2,1) = bai;
            (*D_tran)(2,2) = aii;

            (*D_tran)(3,3) = cai;

            if(Dimension==3)
            {
               (*D_tran)(4,4) = cii;
               (*D_tran)(5,5) = cai;
            }
            break;
         case 13:
            D_tran = new Matrix(size, size);
            (*D_tran)(0,0) = aii;
            (*D_tran)(0,1) = bii;
            (*D_tran)(0,2) = bai;

            (*D_tran)(1,0) = bii;
            (*D_tran)(1,1) = aii;
            (*D_tran)(1,2) = bai;

            (*D_tran)(2,0) = bai;
            (*D_tran)(2,1) = bai;
            (*D_tran)(2,2) = aai;

            (*D_tran)(3,3) = cii;

            if(Dimension==3)
            {
               (*D_tran)(4,4) = cai;
               (*D_tran)(5,5) = cai;
            }
            break;
      }
   }
   /*************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::CalculateTransformMatrixFromNormalVector

    Aufgabe:
      Generate rotation matrices for vectors and matrices based on a given normal vector

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
          const int Dimension  :  Dimension of the real space
          double *Crotv       :  Rotation matrix for vectors
          double *Crotm       :  Rotation matrix for matrices
    Ergebnis:
   - void -

   Programmaenderungen:
   25.11.2009  UJG  first implementation

   *************************************************************************/
   void CSolidProperties::CalculateTransformMatrixFromNormalVector(const int Dimension)
   {
      if(!(Youngs_mode>9&&Youngs_mode<14))        //WW
         return;

      if(Youngs_mode!=10)
      {
         ElasticConstitutiveTransverseIsotropic(Dimension);
         return;
      }

      double e1[3]=
      {
         0.0
      }
      ,e2[3]=
      {
         0.0
      }
      ,e3[3]=
      {
         0.0
      };

      double nx = (*data_Youngs)(4);
      double ny = (*data_Youngs)(5);
      double nz = (*data_Youngs)(6);
      double Disk = 1.0, t1, t2, t3, ax(1.0), ay(0.0), az(0.0);

      t1 = t2 = t3 = 0.;

                                                  // rotation matrix for vectors: UJG 25.11.2009
      Matrix *Crotv = new Matrix(Dimension,Dimension);
      Crotm = new Matrix(Dimension*2,Dimension*2);// rotation matrix for matrices: UJG 25.11.2009

      if(Dimension==3)
      {
         if(nz<-1.0)
         {
            t1    = 0.01745329252*nx;
            t2    = 0.01745329252*ny;
            e3[0] = sin(t1)*cos(t2);
            e3[1] = -sin(t2);
            e3[2] = cos(t1)*cos(t2);
         }
         else if(nz>1.0)
         {
            t1    = 0.01745329252*nx;
            t2    = 0.01745329252*ny;
            e3[0] = sin(t2)*sin(t1);
            e3[1] = sin(t2)*cos(t1);
            e3[2] = cos(t2);
         }
         else
         {
            e3[0] = nx;
            e3[1] = ny;
            e3[2] = nz;
         }

         ax = e3[0];
         ay = e3[1];
         az = e3[2];

         if(ax==0.0)
         {
            e1[0] = 1.0;
            e1[1] = 0.0;
            e1[2] = 0.0;

            e2[0] = 0.0;
            e2[1] = az;
            e2[2] = -ay;

            e3[0] = 0.0;
            e3[1] = ay;
            e3[2] = az;
         }
         else
         {
            Disk = 1.0/sqrt(1.0+((az*az)/(ax*ax)));

            e1[0] = (az/ax)*Disk;
            e1[1] = 0.0;
            e1[2] = -Disk;

            t1    = ay*e1[2];
            t2    = az*e1[0]-ax*e1[2];
            t3    = -ay*e1[0];

            Disk  = t1*t1;
            Disk += t2*t2;
            Disk += t3*t3;

            Disk = 1.0/sqrt(Disk);

            e2[0] = Disk*t1;
            e2[1] = Disk*t2;
            e2[2] = Disk*t3;
         }

         (*Crotv)(0,0) = e1[0];
         (*Crotv)(0,1) = e1[1];
         (*Crotv)(0,2) = e1[2];

         (*Crotv)(1,0) = e2[0];
         (*Crotv)(1,1) = e2[1];
         (*Crotv)(1,2) = e2[2];

         (*Crotv)(2,0) = e3[0];
         (*Crotv)(2,1) = e3[1];
         (*Crotv)(2,2) = e3[2];

         (*Crotm)(0,0) = (*Crotv)(0,0)*(*Crotv)(0,0);
         (*Crotm)(0,1) = (*Crotv)(0,1)*(*Crotv)(0,1);
         (*Crotm)(0,2) = (*Crotv)(0,2)*(*Crotv)(0,2);
         (*Crotm)(0,3) = (*Crotv)(0,0)*(*Crotv)(0,1);
         (*Crotm)(0,4) = (*Crotv)(0,0)*(*Crotv)(0,2);
         (*Crotm)(0,5) = (*Crotv)(0,1)*(*Crotv)(0,2);

         (*Crotm)(1,0) = (*Crotv)(1,0)*(*Crotv)(1,0);
         (*Crotm)(1,1) = (*Crotv)(1,1)*(*Crotv)(1,1);
         (*Crotm)(1,2) = (*Crotv)(1,2)*(*Crotv)(1,2);
         (*Crotm)(1,3) = (*Crotv)(1,0)*(*Crotv)(1,1);
         (*Crotm)(1,4) = (*Crotv)(1,0)*(*Crotv)(1,2);
         (*Crotm)(1,5) = (*Crotv)(1,1)*(*Crotv)(1,2);

         (*Crotm)(2,0) = (*Crotv)(2,0)*(*Crotv)(2,0);
         (*Crotm)(2,1) = (*Crotv)(2,1)*(*Crotv)(2,1);
         (*Crotm)(2,2) = (*Crotv)(2,2)*(*Crotv)(2,2);
         (*Crotm)(2,3) = (*Crotv)(2,0)*(*Crotv)(2,1);
         (*Crotm)(2,4) = (*Crotv)(2,0)*(*Crotv)(2,2);
         (*Crotm)(2,5) = (*Crotv)(2,1)*(*Crotv)(2,2);

         (*Crotm)(3,0) = 2.0*(*Crotv)(0,0)*(*Crotv)(1,0);
         (*Crotm)(3,1) = 2.0*(*Crotv)(0,1)*(*Crotv)(1,1);
         (*Crotm)(3,2) = 2.0*(*Crotv)(0,2)*(*Crotv)(1,2);
         (*Crotm)(3,3) = (*Crotv)(0,0)*(*Crotv)(1,1)+(*Crotv)(1,0)*(*Crotv)(0,1);
         (*Crotm)(3,4) = (*Crotv)(0,0)*(*Crotv)(1,2)+(*Crotv)(0,2)*(*Crotv)(1,0);
         (*Crotm)(3,5) = (*Crotv)(0,1)*(*Crotv)(1,2)+(*Crotv)(0,2)*(*Crotv)(1,1);

         (*Crotm)(4,0) = 2.0*(*Crotv)(0,0)*(*Crotv)(2,0);
         (*Crotm)(4,1) = 2.0*(*Crotv)(0,1)*(*Crotv)(2,1);
         (*Crotm)(4,2) = 2.0*(*Crotv)(0,2)*(*Crotv)(2,2);
         (*Crotm)(4,3) = (*Crotv)(0,0)*(*Crotv)(2,1)+(*Crotv)(0,1)*(*Crotv)(2,0);
         (*Crotm)(4,4) = (*Crotv)(0,0)*(*Crotv)(2,2)+(*Crotv)(2,0)*(*Crotv)(0,2);
         (*Crotm)(4,5) = (*Crotv)(0,1)*(*Crotv)(2,2)+(*Crotv)(0,2)*(*Crotv)(2,1);

         (*Crotm)(5,0) = 2.0*(*Crotv)(1,0)*(*Crotv)(2,0);
         (*Crotm)(5,1) = 2.0*(*Crotv)(1,1)*(*Crotv)(2,1);
         (*Crotm)(5,2) = 2.0*(*Crotv)(1,2)*(*Crotv)(2,2);
         (*Crotm)(5,3) = (*Crotv)(1,0)*(*Crotv)(2,1)+(*Crotv)(1,1)*(*Crotv)(2,0);
         (*Crotm)(5,4) = (*Crotv)(1,0)*(*Crotv)(2,2)+(*Crotv)(1,2)*(*Crotv)(2,0);
         (*Crotm)(5,5) = (*Crotv)(1,1)*(*Crotv)(2,2)+(*Crotv)(2,1)*(*Crotv)(1,2);
      }
      else
      {
         if(ny>1.0)
         {
            e1[0] = cos(0.01745329252*nx);
            e1[1] = sin(0.01745329252*nx);

            e2[0] = -e1[1];
            e2[1] = e1[0];
         }
         else
         {
            e2[0] = nx;
            e2[1] = ny;

            e1[0] = ny;
            e1[1] = -nx;
         }

         (*Crotv)(0,0) = e1[0];
         (*Crotv)(0,1) = e1[1];

         (*Crotv)(1,0) = e2[0];
         (*Crotv)(1,1) = e2[1];

         (*Crotm)(0,0) = (*Crotv)(0,0)*(*Crotv)(0,0);
         (*Crotm)(0,1) = (*Crotv)(1,0)*(*Crotv)(1,0);
         (*Crotm)(0,2) = 0.0;
         (*Crotm)(0,3) = -(*Crotv)(0,0)*(*Crotv)(1,0);

         (*Crotm)(1,0) = (*Crotv)(0,1)*(*Crotv)(0,1);
         (*Crotm)(1,1) = (*Crotv)(1,1)*(*Crotv)(1,1);
         (*Crotm)(1,2) = 0.0;
         (*Crotm)(1,3) = -(*Crotv)(0,1)*(*Crotv)(1,1);

         (*Crotm)(2,0) = 0.0;
         (*Crotm)(2,1) = 0.0;
         (*Crotm)(2,2) = 1.0;
         (*Crotm)(2,3) = 0.0;

         (*Crotm)(3,0) = -2.0*(*Crotv)(0,0)*(*Crotv)(0,1);
         (*Crotm)(3,1) = -2.0*(*Crotv)(1,0)*(*Crotv)(1,1);
         (*Crotm)(3,2) = 0.0;
         (*Crotm)(3,3) = (*Crotv)(0,0)*(*Crotv)(1,1)+(*Crotv)(1,0)*(*Crotv)(0,1);
      }
      delete Crotv;
      Crotv = NULL;

      ElasticConstitutiveTransverseIsotropic(Dimension);
   }
   //
   // WW. 09/02. Compute dilatancy for yield function
   // For plane strain
   double CSolidProperties::GetAngleCoefficent_DP(const double Angle)
   {
      double val = 0.0;
      // Input as a coefficent
      if(Angle<0.0) val = fabs(Angle);
      else                                        // Input as a dialatant angle
      {
         double D_Angle = Angle*PI/180.0;
         double sinA = sin(D_Angle);
         //     val = sinA/sqrt(9.0+4.0*sinA*sinA);
         //     val = 2.0*MSqrt2Over3*sinA/(3.0+sinA);
         val = 2.0*MSqrt2Over3*sinA/(3.0+sinA);   //(3.0-sinA)
         //     val = 2.0*MSqrt2Over3*sinA/(3.0-sinA);
      }
      return val;
   }
   // WW. 09/02. Cumpute yield coefficient, beta
   // For plane strain
   // al = 6.0*c*cos(a)/sqrt(3)/(3-sin(a))
   double CSolidProperties::GetYieldCoefficent_DP(const double Angle)
   {
      double val = 0.0;
      // Input as a coefficent
      if(Angle<0.0||Angle<MKleinsteZahl) val = 1.0;
      else                                        // Input as a dialatant angle
      {
         double D_Angle = Angle*PI/180.0;
         double sinA = sin(D_Angle);
         val = 6.0*cos(D_Angle)/(3.0+sinA);       // (3.0-sinA)
         //     val = 2.0*sqrt(3.0)*cos(D_Angle)/sqrt(9.0+4.0*sinA*sinA);
      }
      return val;
   }

   void CSolidProperties::CalulateCoefficent_DP()
   {
      Y0 = (*data_Plasticity)(0);
      Hard = (*data_Plasticity)(1);
      Al = GetAngleCoefficent_DP((*data_Plasticity)(2));
      Xi = GetAngleCoefficent_DP((*data_Plasticity)(3));
      BetaN = GetYieldCoefficent_DP((*data_Plasticity)(2));
      if(fabs(Al)<MKleinsteZahl&&fabs(Xi)< MKleinsteZahl) BetaN = 1.0;
      BetaN *= sqrt(2.0/3.0);
      Hard_Loc = (*data_Plasticity)(4);
   }

   /**************************************************************************
     ROCKFLOW - Funktion: StressIntegrationDP

      Aufgabe:
      Computing the stresses at a point and return the plastical status of this
      point (Return mapping method)

      Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
        long index: Elementnummer
   double *TryStress                   :   Try incremental stresses as input
   New stresses as output
   const int GPi, const int GPj        :   Local indeces of Gauss Points

   const double G, double K            :   Shear modulus, (2*lambda+2*G)/3
   const double Al, const double Xi,
   const double Y0, const double BetaN :   Coefficient for Drucker-Prager model
   double* dPhi                        :   Plastic multiplier.
   double &jt                          :   enhanced parameters

   const int Update                    :   Indicator to store the stress or not.

   Ergebnis:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   09/2002   WW  Erste Version
   02/2004   WW  Modify for the 3D case
   08/2004   WW  Set As a member of material class
   03/2007   WW  Multi-yield surface
   **************************************************************************/
   bool CSolidProperties::StressIntegrationDP(const int GPiGPj,
      const ElementValue_DM *ele_val, double *TryStress,
      double& dPhi, const int Update)
   {
      int i = 0;
      double I1 = 0.0;
      double ep, ep0;
      double F = 0.0, F0 = 0.0;                   //, yl;
      //  double RF0 = 0.0;
      double sqrtJ2 = 0.0;
      double Beta = 0.0;
      double p3, normXi, Jac, err_corner=0.0, fac=0.0;
      int max_ite = 10;
      int ite=0;
      // double dstrs[6];

      bool isLoop = true;                         // Used only to avoid warnings with .net
      bool ploading = false;

      const int Size = ele_val->Stress->Rows();
      //static double DevStress[6];

      int Dim = 2;
      if(Size>4) Dim = 3;

      // Get the total effective plastic strain
      ep = (*ele_val->pStrain)(GPiGPj);

      for(i=0; i<Size; i++)
      {
         //     dstrs[i] = TryStress[i];  // d_stress
         TryStress[i] += (*ele_val->Stress)(i, GPiGPj);
         devS[i] = TryStress[i];
      }
      // I_tr
      I1 = DeviatoricStress(devS);
      // s_tr
      sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));
      //
      normXi = sqrtJ2;
      p3 = I1;
      /* If yield, compute plastic multiplier dPhi */
      dPhi = 0.0;
      F0 = sqrtJ2 + Al*I1;
      F = F0 - BetaN*(Y0+Hard*ep);
      //if(pcs_deformation==1) F=-1.0;
      // yl = sqrt(TensorMutiplication2(devS, dstrs, Dim))/sqrtJ2;
      // yl += (dstrs[0]+dstrs[1]+dstrs[2])*Al;
      // if(yl<=0.0)
      //   F = -1.0;
      //
      if(F0<=(*ele_val->y_surface)(GPiGPj))       // unloading
         F = -1.0;
      //
      if(F>0.0&&(!PreLoad))                       // in yield status
      {
         ploading = true;
         //err = 1.0e+5;
         ep0 = ep;

         //TEST
         // Check the corner region
         err_corner = 4.5*Xi*K*normXi/G+BetaN*
            (0.5*Hard*sqrt(1.0+3.0*Xi*Xi)*normXi/G +Y0+Hard*ep )/Al;
         // Multi-surface
         if(p3>err_corner)
         {
            // RF0 = F;
            dl2 = 0.0;
            while(isLoop)
            {
               ite++;
               // dl1
               dPhi = 0.5*normXi/G;
               fac = sqrt(dPhi*dPhi+3.0*Xi*Xi*(dPhi+dl2)*(dPhi+dl2));
               Jac = 9.0*Xi*K+3.0*Xi*Xi*BetaN*Hard*(dPhi+dl2)/(fac*Al);
               F = 9.0*Xi*K*(dPhi+dl2)+BetaN*(Y0+Hard*(ep0+fac))/Al-p3;
               dl2 -= F/Jac;
               if(fabs(F)<1000.0*Tolerance_Local_Newton) break;
               if(ite>max_ite) break;
            }
            ep = ep0 + fac;
            for(i=0; i<3; i++)
               TryStress[i] =  I1/3.0-3.0*(dPhi+dl2)*K*Xi;
            for(i=3; i<Size; i++)
               TryStress[4] = 0.0;
         }
         else
         {
            // Local Newton-Raphson procedure
            // If non-perfect plasticity, the below line has to be change
            Jac = -2.0*G-9.0*K*Al*Xi
               -BetaN*Hard*sqrt(1.0+3.0*Xi*Xi);
            // RF0 = F;
            while(isLoop)
            {
               ite++;
               if(ite>max_ite) break;
               if(F<0.0||fabs(F)<10.0*Tolerance_Local_Newton) break;
               //if(err<TolLocalNewT) break;
               dPhi -= F/Jac;
               //
               p3 = I1 - 9.0*dPhi*Xi*K;
               normXi = sqrtJ2 - 2.0*G*dPhi;
               ep = ep0 + dPhi*sqrt(1.0+3.0*Xi*Xi);
               F0 =  normXi + Al*p3;
               F = F0 - BetaN*(Y0+Hard*ep);
               /* Jac = fun(); if non-linear hardening is involved */
               //err = fabs(F)/RF0;
            }
            // update stress
            Beta = 1.0-2.0*dPhi*G/sqrtJ2;
            for(i=0; i<Size; i++)
               TryStress[i] = Beta*devS[i];
            for(i=0; i<3; i++)
               TryStress[i] += I1/3.0-3.0*dPhi*K*Xi;
         }
      }
      else
      {
         //
         for(i=0; i<Size; i++)
            TryStress[i] = devS[i];
         //
         for(i=0; i<3; i++)
            TryStress[i] += I1/3.0;
      }
      // Save the current stresses
      if(Update>0)
      {
         if(dPhi>0.0)
            (*ele_val->pStrain)(GPiGPj) = ep;
         (*ele_val->y_surface)(GPiGPj) = F0;
      }
      return ploading;
   }

   /**************************************************************************
     ROCKFLOW - Funktion: DirectStressIntegrationDP
      Computing the stresses at a point and return the plastical status of this
      point by direct integration.

      Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
        long index: Elementnummer
        double *TryStress                   :   Try incremental stresses as input
                                                New stresses as output
        const int GPi, const int GPj        :   Local indeces of Gauss Points

   const double G, double K            :   Shear modulus, (2*lambda+2*G)/3
   const double Al, const double Xi,
   const double Y0, const double BetaN :   Coefficient for Drucker-Prager model
   double* dPhi                        :   Plastic multiplier.
   double &jt                          :   enhanced parameters

   const int Update                    :   Indicator to store the stress or not.

   return:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   02/2006   WW  Erste Version
   **************************************************************************/
   bool CSolidProperties::DirectStressIntegrationDP(const int GPiGPj,
      const ElementValue_DM *ele_val, double *TryStress, const int Update)
   {
      int i = 0, m=0, m_max=100;
      double I1 = 0.0;
      double sy0, sy, ep, dlambda=0.0;
      double F = 0.0, yy;                         //, yl;
      double R=1.0;
      double sqrtJ2 = 0.0;
      double A_H = 0.0, domA=0.0;
      double dstrs[6];
      bool ploading = false;
      const int Size = ele_val->Stress->Rows();
      //static double DevStress[6];

      int Dim = 2;
      if(Size>4) Dim = 3;

      // Get the total effective plastic strain
      ep = (*ele_val->pStrain)(GPiGPj);

      for(i=0; i<Size; i++)
      {
         dstrs[i] = TryStress[i];                 // d_stress
                                                  // stress_0
         TryStress[i] = (*ele_val->Stress)(i, GPiGPj);
         devS[i] = TryStress[i] + dstrs[i];
      }

      I1 = DeviatoricStress(devS);
      sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));
      //
      sy = sqrtJ2 + Al*I1;
      yy = BetaN*(Y0+Hard*ep);
      F = sy - yy;
      sy0 = (*ele_val->y_surface)(GPiGPj);
      //yl = sqrt(TensorMutiplication2(devS, dstrs, Dim))/sqrtJ2;
      //yl += (dstrs[0]+dstrs[1]+dstrs[2])*Al;
      //if(yl<=0.0)
      //  F = -1.0;
      if(sy<=sy0)                                 // unloading
         F = -1.0;
      if(F>0.0&&(!PreLoad))                       // in yield status
      {
         if(ep<MKleinsteZahl)                     // Elastic in previous load step
            R = F/(sy-sy0);
         m=(int)(8.0*F/yy)+1;
         for(i=0; i<Size; i++)
         {
            TryStress[i] += (1.0-R)*dstrs[i];
            dstrs[i] *= R/(double)m;
         }
         if(m>m_max)
            m=m_max;
         // sub-inrement
         while(m>0)
         {
            // Compute dlamda
            A_H = BetaN*Hard*sqrt(1+3.0*Xi*Xi);   // Hard: if it is not constant....
            for(i=0; i<Size; i++)
               devS[i] = TryStress[i];
            I1 = DeviatoricStress(devS);
            sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));
            for(i=0; i<Size; i++)
            {
               devS[i] /= sqrtJ2;
               dFds[i] = devS[i];
            }
            for(i=0; i<3; i++)
               dFds[i] += Al;
            // dlambda
            dlambda = 0.0;
            domA = A_H+2.0*G+9.0*Al*Xi*K;
            for(i=0; i<Size; i++)
               dlambda += dFds[i]*dstrs[i];
            dlambda /= domA;
            if(dlambda<0.0) dlambda = 0.0;
            ep += dlambda*sqrt(1.0+3.0*Xi*Xi);
            // Update stress
            for(i=0; i<Size; i++)
               TryStress[i] += dstrs[i]-2.0*dlambda*G*devS[i];
            dlambda *= 3.0*Xi*K;
            for(i=0; i<3; i++)
               TryStress[i] -= dlambda;
            m--;
         }
         for(i=0; i<Size; i++)
            devS[i] = TryStress[i];
         I1 = DeviatoricStress(devS);
         sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));
         sy = sqrtJ2 + Al*I1;
         yy = BetaN*(Y0+Hard*ep);
         R=1.0;
         if(sy>yy)
            R = yy/sy;
         for(i=0; i<Size; i++)
         {
            TryStress[i] *= R;
            devS[i] *= R/sqrtJ2;
         }
         sy *= R;
         ploading = true;
      }
      else
      {
         for(i=0; i<Size; i++)
            TryStress[i] += dstrs[i];
      }
      // Save the current stresses
      if(Update>0)
      {
         (*ele_val->pStrain)(GPiGPj) = ep;
         (*ele_val->y_surface)(GPiGPj) = sy;
      }
      return ploading;
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::ConsistentTangentialDP

      Local assembly of elasto-plastic tangential matrix C^ep
      (Drucker-Prager model)

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E
      double *Dep            : Consistent tangential matrix
      const Dim              : Space dimension

   Programmaenderungen:
   10/2006   WW  Erste Version
   02/2006   WW  programmed
   **************************************************************************/
   void CSolidProperties::TangentialDP(Matrix *Dep)
   {
      int i, j, Size;
      double domA;
      //
      Size=Dep->Rows();
      domA = BetaN*Hard*sqrt(1+3.0*Xi*Xi);        // Hard: if it is not constant....
      //
      for(i=0; i<Size; i++)
      {
         D_dFds[i] = 2.0*G*devS[i];
         D_dGds[i] = 2.0*G*devS[i];
      }
      for(i=0; i<3; i++)
      {
         D_dFds[i] += 3.0*Al*K;
         D_dGds[i] += 3.0*Xi*K;
      }
      //
      domA += 2.0*G+9.0*Al*Xi*K;
      //
      for(i=0; i<Size; i++)
      {
         for(j=0; j<Size; j++)
            (*Dep)(i,j) -= D_dGds[i]*D_dFds[j]/domA;
      }

      //Dep->Write();
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::ConsistentTangentialDP

    Aufgabe:
      Local assembly of elasto-plastic consistent tangential matrix C^ep
      (Drucker-Prager model)

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E
      double *Dep            : Consistent tangential matrix
      const double dPhi      : Plastic multiplier
   const Dim              : Space dimension
   Ergebnis:
   - double* - Matrix of 4x4, stored in an array

   Programmaenderungen:
   10/2002   WW  Erste Version
   02/2004   WW  Modification for the 3D case
   03/2007   WW  Multi-phase yield surface
   **************************************************************************/
   void CSolidProperties::ConsistentTangentialDP(Matrix *Dep, const double dPhi, const int Dim)
   {
      double s11, s22, s12, s33, s13, s23;
      double NormX = 0.0;
      double d = 0.0;
      double c1,c2,c3,c4;
      //
      s11 = devS[0];
      s22 = devS[1];
      s33 = devS[2];
      s12 = devS[3];
      s13 = 0.0;
      s23 = 0.0;
      if(Dim==3)
      {
         s13 =  devS[4];
         s23 =  devS[5];
         NormX = sqrt(s11*s11+s22*s22+s33*s33
            +2.0*s12*s12+2.0*s13*s13+2.0*s23*s23);
      }
      else
         NormX = sqrt(s11*s11+s22*s22+s33*s33+2.0*s12*s12);
      //
      if(dl2>0.0)                                 // Multi-surface
      {
         c1 = Xi*BetaN*Hard*(dPhi+dl2);
         c3 = K*c1/(c1+3.0*Al*K*sqrt(dPhi*dPhi+3.0*Xi*Xi*(dPhi+dl2)*(dPhi+dl2)));
         c2 = 0.5*c3/(Xi*G*(dPhi+dl2));
         (*Dep) = 0.0;
         //
         // Row 1
         (*Dep)(0,0) = c3+c2*s11;
         (*Dep)(0,1) = c3+c2*s22;
         (*Dep)(0,3) = c2*s12;
         // Row 2
         (*Dep)(1,0) = c3+c2*s11;
         (*Dep)(1,1) = c3+c2*s22;
         (*Dep)(1,3) = c2*s12;
         // Row 3
         (*Dep)(2,0) = c3+c2*s11;
         (*Dep)(2,1) = c3+c2*s22;
         (*Dep)(2,3) = c2*s12;
         // Row 4
         if(axisymmetry||Dim==3)
         {
            (*Dep)(0,2) = c3+c2*s33;
            (*Dep)(1,2) = c3+c2*s33;
            (*Dep)(2,2) = c3+c2*s33;
         }
         if(Dim==3)
         {
            (*Dep)(0,4) = c2*s13;
            (*Dep)(0,5) = c2*s23;
            (*Dep)(1,4) = c2*s13;
            (*Dep)(1,5) = c2*s23;
            (*Dep)(2,4) = c2*s13;
            (*Dep)(2,5) = c2*s23;
         }
      }
      else
      {
         //
         d=9.0*K*Al*Xi+2.0*G+BetaN*Hard*
            sqrt(1.0+3.0*Xi*Xi);

         c1 = 2.0*G*(1.0-2.0*G*dPhi/NormX);
         c2 = (1.0-9.0*Al*Xi*K/d)*K;
         c3 = -6.0*G*K/(NormX*d);
         c4 = -4.0*G*G*(1.0/d-dPhi/NormX)/(NormX*NormX);

         switch(Dim)
         {
            case 2:                               // 2D
               // Row 1
               (*Dep)(0,0)
                  =  2.0*c1/3.0+c2+c3*(Al+Xi)*s11+c4*s11*s11;
               (*Dep)(0,1)
                  = -c1/3.0+c2+c3*(Xi*s22+Al*s11)+c4*s11*s22;
               (*Dep)(0,2)
                  = 0.0;                          //-c1/3.0+c2+c3*(Xi*s33+Al*s11)+c4*s11*s33;
               (*Dep)(0,3)
                  = c3*Xi*s12+c4*s11*s12;

               // Row 2
               (*Dep)(1,0)
                  = -c1/3.0+c2+c3*(Xi*s11+Al*s22)+c4*s11*s22;
               (*Dep)(1,1)
                  =  2.0*c1/3.0+c2+c3*(Al+Xi)*s22+c4*s22*s22;
               (*Dep)(1,2)
                  = 0.0;                          //-c1/3.0+c2+c3*(Xi*s33+Al*s22)+c4*s33*s22;
               (*Dep)(1,3)
                  = c3*Xi*s12+c4*s22*s12;

               // Row 3
               (*Dep)(2,0)
                  = 0.0;                          //-c1/3.0+c2+c3*(Xi*s11+Al*s33)+c4*s11*s33;
               (*Dep)(2,1)
                  = 0.0;                          // -c1/3.0+c2+c3*(Xi*s22+Al*s33)+c4*s22*s33;
               (*Dep)(2,2)
                  = 0.0;                          // 2.0*c1/3.0+c2+c3*(Al+Xi)*s33+c4*s33*s33;
               (*Dep)(2,3)
                  = 0.0;                          // c3*Xi*s12+c4*s33*s12;

               // Row 4
               (*Dep)(3,0)
                  = c3*Al*s12+c4*s12*s11;
               (*Dep)(3,1)
                  = c3*Al*s12+c4*s12*s22;
               (*Dep)(3,2)
                  = 0.0;                          //c3*Al*s12+c4*s12*s33;
               (*Dep)(3,3)
                  =  c1/2.0+c4*s12*s12;

               if(axisymmetry)
               {
                  (*Dep)(0,2)
                     = -c1/3.0+c2+c3*(Xi*s33+Al*s11)+c4*s11*s33;
                  // Row 2
                  (*Dep)(1,2)
                     = -c1/3.0+c2+c3*(Xi*s33+Al*s22)+c4*s33*s22;
                  // Row 3
                  (*Dep)(2,0)
                     = -c1/3.0+c2+c3*(Xi*s11+Al*s33)+c4*s11*s33;
                  (*Dep)(2,1)
                     = -c1/3.0+c2+c3*(Xi*s22+Al*s33)+c4*s22*s33;
                  (*Dep)(2,2)
                     =  2.0*c1/3.0+c2+c3*(Al+Xi)*s33+c4*s33*s33;
                  (*Dep)(2,3)
                     = c3*Xi*s12+c4*s33*s12;
                  // Row 4
                  (*Dep)(3,2)
                     = c3*Al*s12+c4*s12*s33;
               }
               break;
            case 3:                               // 3D
               // Row 1
               (*Dep)(0,0)
                  =  2.0*c1/3.0+c2+c3*(Al+Xi)*s11+c4*s11*s11;
               (*Dep)(0,1)
                  = -c1/3.0+c2+c3*(Xi*s22+Al*s11)+c4*s11*s22;
               (*Dep)(0,2)
                  = -c1/3.0+c2+c3*(Xi*s33+Al*s11)+c4*s11*s33;
               (*Dep)(0,3)
                  = c3*Xi*s12+c4*s11*s12;
               (*Dep)(0,4)
                  = c3*Xi*s13+c4*s11*s13;
               (*Dep)(0,5)
                  = c3*Xi*s23+c4*s11*s23;
               // Row 2
               (*Dep)(1,0)
                  = -c1/3.0+c2+c3*(Xi*s11+Al*s22)+c4*s11*s22;
               (*Dep)(1,1)
                  =  2.0*c1/3.0+c2+c3*(Al+Xi)*s22+c4*s22*s22;
               (*Dep)(1,2)
                  = -c1/3.0+c2+c3*(Xi*s33+Al*s22)+c4*s33*s22;
               (*Dep)(1,3)
                  = c3*Xi*s12+c4*s22*s12;
               (*Dep)(1,4)
                  = c3*Xi*s13+c4*s22*s13;
               (*Dep)(1,5)
                  = c3*Xi*s23+c4*s22*s23;
               // Row 3
               (*Dep)(2,0)
                  = -c1/3.0+c2+c3*(Xi*s11+Al*s33)+c4*s11*s33;
               (*Dep)(2,1)
                  =  -c1/3.0+c2+c3*(Xi*s22+Al*s33)+c4*s22*s33;
               (*Dep)(2,2)
                  =  2.0*c1/3.0+c2+c3*(Al+Xi)*s33+c4*s33*s33;
               (*Dep)(2,3)
                  =  c3*Xi*s12+c4*s33*s12;
               (*Dep)(2,4)
                  =  c3*Xi*s13+c4*s33*s23;
               (*Dep)(2,5)
                  =  c3*Xi*s23+c4*s33*s23;
               // Row 4
               (*Dep)(3,0)
                  = c3*Al*s12+c4*s12*s11;
               (*Dep)(3,1)
                  = c3*Al*s12+c4*s12*s22;
               (*Dep)(3,2)
                  = c3*Al*s12+c4*s12*s33;
               (*Dep)(3,3)
                  =  c1/2.0+c4*s12*s12;
               (*Dep)(3,4)
                  =  c4*s12*s13;
               (*Dep)(3,5)
                  =  c4*s12*s23;
               // Row 5
               (*Dep)(4,0)
                  = c3*Al*s13+c4*s13*s11;
               (*Dep)(4,1)
                  = c3*Al*s13+c4*s13*s22;
               (*Dep)(4,2)
                  = c3*Al*s13+c4*s13*s33;
               (*Dep)(4,3)
                  =  c4*s13*s12;
               (*Dep)(4,4)
                  =  c1/2.0+c4*s13*s13;
               (*Dep)(4,5)
                  =  c4*s13*s23;
               // Row 6
               (*Dep)(5,0)
                  = c3*Al*s23+c4*s23*s11;
               (*Dep)(5,1)
                  = c3*Al*s23+c4*s23*s22;
               (*Dep)(5,2)
                  = c3*Al*s23+c4*s23*s33;
               (*Dep)(5,3)
                  =  c4*s23*s12;
               (*Dep)(5,4)
                  =  c4*s23*s13;
               (*Dep)(5,5)
                  =  c1/2.0+c4*s23*s23;
               break;
         }
      }
   }
   //-------------------------------------------------------------------------
   void CSolidProperties::AllocateMemoryforSYS()
   {
      d2G_dSdS = new Matrix(6,6);
      d2G_dSdM = new Matrix(6,4);
      // Stresses, 0-5, xi, 6-10, mat, 11-17, plastci multiplier, 18.
      LocalJacobi = new Matrix(19,19);
      inv_Jac = new Matrix(18,18);
      sumA_Matrix = new Matrix(18,6);
      rhs_l = new double[18];
      x_l = new double[18];
      Li = new int[18];
   }

   void CSolidProperties::ResizeMatricesSYS(const int Dim)
   {
      if(Dim==2)
      {
         d2G_dSdS->LimitSize(4,4);
         d2G_dSdM->LimitSize(4,4);
         LocalJacobi->LimitSize(15,15);
         inv_Jac->LimitSize(14,14);
         sumA_Matrix->LimitSize(14,4);
      }
      else
      {
         d2G_dSdS->LimitSize(6,6);
         d2G_dSdM->LimitSize(6,4);
         LocalJacobi->LimitSize(19,19);
         inv_Jac->LimitSize(18,18);
         sumA_Matrix->LimitSize(18,6);
      }
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::CalStress_and_TangentialMatrix_SYS

    Aufgabe:
        Using the substteping techniques to compute the integral of stresses and
      the consistent tangential matrix.
      (Weimar's model)

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
   long index: Elementnummer
   const int GPi, const int GPj        :   Local indeces of Gauss Points
   double *D_ep                        :   the consistent tangential matrix
   const double *De                    :   Elastic constitutive matrix
   double *dStress                     :   Incremental stresses as input
   New stress as output
   const int Update                    :   0: Do not save Gauss points values
   1: Save Gauss points values

   Ergebnis:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   04/2003   WW  Erste Version
   02/2004   WW  Modify for the 3D case
   08/2004   WW  Set as a member of CSolidProperties
   **************************************************************************/
   int CSolidProperties::CalStress_and_TangentialMatrix_SYS
      (const int GPiGPj, const ElementValue_DM *ele_val,
      const Matrix *De,  Matrix *D_ep, double *dStress,
      const int Update)
   {
      const int LengthMat=7;
      const int LengthStrs=De->Cols();
      int dim = 2;
      const int minSub=4;
      const int maxSub=1000;
      const double den=0.01;                      //0.01
      const int MaxIter = 1000;

      const double TolF = Tolerance_Local_Newton;

      double ErrLoc, ErrLoc0;
      int PLASTIC, subPLASTIC;

      int NPStep;
      int preSub = 1;

      int i,j,k,l,iSub, elasticSub;
      double factor, dSub;

      double I_n = 0.0;
      double I_n1 = 0.0;

      double I_p2 = 0.0;
      double I_p4 = 0.0;

      double ep = 0.0;
      double PHI = 0.0;
      double PSI = 0.0;
      double ang = 0.0;
      double dlmd = 0.0;

      const double psi1= (*data_Plasticity)(14);
      const double Ch= (*data_Plasticity)(16);
      const double Cd= (*data_Plasticity)(17);
      const double br= (*data_Plasticity)(18);
      const double mr= (*data_Plasticity)(19);

      double Co=0.0;
      double t1,t2,c1,c2, dW;
      double damping;
      double NormR0, NormR1;                      //, NormU0, NormU1;

      double F = 0.0;

      static double  Stress_Inv[7];

      static double  Stress_n[6];
      static double  Stress_n1[6];
      static double  nStress[6];
      static double  xi_n[6];
      static double  xi_n1[6];
      static double  Mat_n[7];
      static double  Mat_n1[7];
      static double  supMat[7];

      static double  dF_dS[6];
      static double  dF_dM[7];
      static double  dG_dS[6];
      const int LocDim=LocalJacobi->Cols();

      // Limits of material parameters
      // Get material parameters of the previous step
      for(i=0; i<LengthMat; i++)
      {
         supMat[i] = (*data_Plasticity)(i+LengthMat);
         Mat_n[i] =  (*ele_val->MatP)(i, GPiGPj);
      }

      // Get the converged stresses of the previous step
      for(i=0; i<LengthStrs; i++)
         Stress_n[i] = (*ele_val->Stress)(i,GPiGPj);

      // Get the total effective plastic strain
      ep = (*ele_val->pStrain)(GPiGPj);

      // Get the rotational hardening varriables of the previous step
      for(i=0; i<LengthStrs-1; i++)
         xi_n[i] = (*ele_val->xi)(i,GPiGPj);
      xi_n[2] = -xi_n[0]-xi_n[1];                 //x_ii*delta_ii = 0

      l = 0;

      // Total elastic try
      for(i=0; i<LengthStrs; i++)
      {
         Stress_n1[i] = Stress_n[i]+dStress[i];
         nStress[i] = Stress_n1[i];
      }

      I_n1 = DeviatoricStress(nStress);
      for(i=0; i<LengthStrs; i++)
      {
         nStress[i]  -= I_n1*xi_n[i]/3.0;
         xi_n1[i] = xi_n[i];
      }

      for(i=0; i<LengthMat; i++)
         Mat_n1[i] = Mat_n[i];

      Stress_Inv[1] = 0.5*TensorMutiplication2(nStress, nStress, dim);
      Stress_Inv[2] = TensorMutiplication3(nStress, nStress, nStress, dim);

      I_p2 = I_n1*I_n1;
      I_p4 = I_p2*I_p2;

      if(Stress_Inv[1]>0.0)
      {
         ang = Stress_Inv[2]/pow(Stress_Inv[1],1.5);

         PHI = Stress_Inv[1]*pow(1+Mat_n[5]*ang, Mat_n[6])+0.5*Mat_n[0]*I_p2
            +Mat_n[2]*Mat_n[2]*I_p4;

         /* Compute yield surface */
         F = sqrt(PHI) + Mat_n[1]*I_n1 +Mat_n[3]*I_p2- Mat_n[4];
      }
      else
         F=-1.0;

      if(pcs_deformation==1) F=-1.0;

      PLASTIC = 0;
      if(F>TolF&&!PreLoad)                        /* In Yield Status */
      {
         PLASTIC=1;
         subPLASTIC=0;

         /* Compute the ratio of elastic increment */
         I_n = DeviatoricStress(Stress_n);
         for(i=0; i<LengthStrs; i++)
            Stress_n[i] -= I_n*xi_n[i]/3.0;

         /* Compute the number of sub step */
         Co = sqrt(2.0*Stress_Inv[1]*TensorMutiplication2(Stress_n, Stress_n, dim));
         if(Co>0.0)
            preSub = 1+(int)(acos(TensorMutiplication2(Stress_n, nStress, dim)/Co)/den);
         else
            preSub = minSub;
         if(preSub>maxSub) preSub=maxSub;
         if(preSub<minSub) preSub=minSub;

         //Recover stress at n_th step
         for(i=0; i<LengthStrs; i++)
            Stress_n[i] += I_n*xi_n[i]/3.0;
         for(i=0; i<3; i++)                       //principle stress
            Stress_n[i] += I_n/3.0;

         dSub = 1.0/(float)preSub;

         // Approximate the yield surface
         elasticSub = 0;
         // The previous solution as the initial approximation
         for(i=0; i<LengthMat; i++)
            Mat_n1[i] = Mat_n[i];
         for(iSub=0; iSub<preSub; iSub++)
         {
            // Elastic try in each sub step
            for(i=0; i<LengthStrs; i++)
            {
               // Elastic try in each sub step as the initial approximation
               Stress_n1[i] = Stress_n[i] + (iSub+1.0)*dSub*dStress[i];
               // The previous solution as the initial approximation
               xi_n1[i] = xi_n[i];
               nStress[i] = Stress_n1[i];
            }
            I_n1 = DeviatoricStress(nStress);
            for(i=0; i<LengthStrs; i++)
               nStress[i] -= I_n1*xi_n1[i]/3.0;

            Stress_Inv[0] = I_n1;
            Stress_Inv[1] = 0.5*TensorMutiplication2(nStress, nStress,dim);
            Stress_Inv[2] = TensorMutiplication3(nStress, nStress, nStress, dim);

            I_p2 = I_n1*I_n1;
            I_p4 = I_p2*I_p2;
            ang = Stress_Inv[2]/pow(Stress_Inv[1],1.5);
            PHI = sqrt(Stress_Inv[1]*pow(1+Mat_n1[5]*ang, Mat_n1[6])+0.5*Mat_n1[0]*I_p2
               +Mat_n1[2]*Mat_n1[2]*I_p4);
            // The yield function
            F =  PHI + Mat_n1[1]*I_n1 +Mat_n1[3]*I_p2- Mat_n1[4];
            if(F>TolF)
            {
               elasticSub = iSub;
               Co = dSub*(float)elasticSub;

               /*
               for(i=0; i<LengthStrs; i++)
               {
                  Stress_n[i] += Co*dStress[i];
                  dStress[i] *= 1-Co;
               }
               */
               break;
            }

         }

         subPLASTIC=0;
         /* Compute stresses by substepping  */
         //for(iSub=0; iSub<preSub; iSub++)
         for(iSub=elasticSub; iSub<preSub; iSub++)
         {
            factor = dSub;
            if(iSub==elasticSub) factor += Co;
            /* Elastic try in each sub step */
            for(i=0; i<LengthStrs; i++)
            {
               /* Elastic try in each sub step as the initial approximation*/
               Stress_n1[i] = Stress_n[i] + factor*dStress[i];
               /* The previous solution as the initial approximation */
               xi_n1[i] = xi_n[i];

               nStress[i] = Stress_n1[i];
            }

            I_n1 = DeviatoricStress(nStress);
            for(i=0; i<LengthStrs; i++)
               nStress[i] -= I_n1*xi_n1[i]/3.0;

            /* The previous solution as the initial approximation */
            for(i=0; i<LengthMat; i++)
               Mat_n1[i] = Mat_n[i];

            /* Local Newton-Raphson method */
            ErrLoc0=1.0e8;
            ErrLoc=1.0e8;
            //if(iSub==0)
            dlmd=0.0;
            NPStep=0;
            //NormU0 = 1.0;
            NormR0 = 1.0;
            NormR1 = 1.0;
            while(ErrLoc>Tolerance_Local_Newton)
            {
               NPStep++;

               Stress_Inv[0] = I_n1;
               Stress_Inv[1] = 0.5*TensorMutiplication2(nStress, nStress,dim);
               Stress_Inv[2] = TensorMutiplication3(nStress, nStress, nStress, dim);

               I_p2 = I_n1*I_n1;
               I_p4 = I_p2*I_p2;
               ang = Stress_Inv[2]/pow(Stress_Inv[1],1.5);
               PHI = sqrt(Stress_Inv[1]*pow(1+Mat_n1[5]*ang, Mat_n1[6])+0.5*Mat_n1[0]*I_p2
                  +Mat_n1[2]*Mat_n1[2]*I_p4);
               PSI = sqrt(Stress_Inv[1]/psi1+0.5*Mat_n1[0]*I_p2+Mat_n1[2]*Mat_n1[2]*I_p4);

               /* The yield function */
               F =  PHI + Mat_n1[1]*I_n1 +Mat_n1[3]*I_p2- Mat_n1[4];

               if(F<=TolF&&NPStep==1&&subPLASTIC==0)
                  break;
               else
               {
                  if(NPStep==1) subPLASTIC++;
                  Stress_Inv[3] = ang;
                  Stress_Inv[4] = PHI;
                  Stress_Inv[5] = PSI;
                  Stress_Inv[6] = PSI*PSI*PSI;

                  dG_dNStress(dG_dS, nStress, Stress_Inv, Mat_n1, LengthStrs);
                  dG__dNStress_dNStress(nStress, Stress_Inv, Mat_n1, LengthStrs);
                  dG_dSTress_dMat(nStress, Stress_Inv, Mat_n1, LengthStrs);

                  //--------------- Local Jacibin ------------------
                  (*LocalJacobi) = 0.0;

                  // dr_stress/d... --------------------------
                  for(i=0; i<LengthStrs; i++)
                  {
                     l = i;
                     // dr_stress/dxi and dr_stress/dM  -------
                     for(j=0; j<LengthStrs; j++)
                     {
                        c1=1.0;
                        if(j>2) c1=2.0;
                        t1=0.0;
                        t2=0.0;
                        for(k=0; k<LengthStrs; k++)
                        {
                           c2=1.0;

                           if(k>2) c2=2.0;
                           t1 += c1*c2*(*De)(i,k)*(*d2G_dSdS)(k,j);
                           if(j>=4) continue;
                           t2 += c2*(*De)(i,k)*(*d2G_dSdM)(k,j);
                        }

                        // dG_dStress_dXi
                        if(br>0.0)
                        {
                           if(j==2)               // xi_33=-(xi_11+xi_22+x_33)
                           {
                              (*LocalJacobi)(l,LengthStrs) += dlmd*t1*I_n1/3.0;
                              (*LocalJacobi)(l,LengthStrs+1) += dlmd*t1*I_n1/3.0;
                           }
                           else if(j<2)
                              (*LocalJacobi)(l,LengthStrs+j) += -dlmd*t1*I_n1/3.0;
                           else if(j>2)
                              (*LocalJacobi)(l, LengthStrs+j-1) += -dlmd*t1*I_n1/3.0;
                        }

                        // dG_dStress_dMat
                        if((Ch>0.0)&&(Cd>0.0))
                           (*LocalJacobi)(l,2*LengthStrs-1+j) += dlmd*t2;
                     }

                     // dG_dStress_dLambda
                     t1 = 0.0;
                     for(k=0; k<LengthStrs; k++)
                     {
                        c2=1.0;
                        if(k>2) c2=2.0;
                        t1 += c2*(*De)(i,k)*dG_dS[k];
                     }
                     (*LocalJacobi)(l, LocDim-1) += t1;
                  }

                  dW=0.0;
                  for(k=0; k<LengthStrs; k++)
                  {
                     c2=1.0;
                     if(k>2) c2=2.0;
                     dW += c2*Stress_n1[k]*dG_dS[k];
                  }

                  // dr_M/d... --------------------------
                  for(i=0; i<LengthMat; i++)
                  {
                     l = i+2*LengthStrs-1;
                     Co = Ch;
                     if(i>4) Co = Cd;
                     for(j=0; j<LengthMat; j++)
                     {
                        if(i==j) (*LocalJacobi)(l, 2*LengthStrs-1+j) += 1.0+dlmd*Co*dW;
                     }

                     // dr_m/dnStress --------------------------
                     for(j=0; j<LengthStrs; j++)
                     {
                        c1=1.0;
                        if(j>2) c1=2.0;
                        t1=0.0;
                        for(k=0; k<LengthStrs; k++)
                        {
                           c2=1.0;
                           if(k>2) c2=2.0;
                           t1 += c1*c2*Stress_n1[k]*(*d2G_dSdS)(k,j);
                        }
                        t2 = -dlmd*Co*(supMat[i]-Mat_n1[i])*t1;

                        // d(.)/dxi
                        if(br>0.0)
                        {
                           if(j==2)               //xi_33 = -(xi_11+xi_22)
                           {
                              (*LocalJacobi)(l,LengthStrs) += I_n1*t2/3.0;
                              (*LocalJacobi)(l,LengthStrs+1) += I_n1*t2/3.0;
                           }
                           else if(j<2)
                              (*LocalJacobi)(l,LengthStrs+j) += -I_n1*t2/3.0;
                           else if(j>2)
                              (*LocalJacobi)(l,LengthStrs+j-1) += -I_n1*t2/3.0;
                        }
                     }

                     // dr_m/dMat --------------------------
                     for(j=0; j<LengthStrs; j++)
                     {
                        t1=0.0;                   // Stress times d2G_dSdM
                        for(k=0; k<4; k++)
                        {
                           c1=1.0;
                           if(k>2) c1=2.0;
                           t1 += c1*Stress_n1[k]*(*d2G_dSdM)(k,j);
                        }
                        t2 = (supMat[i]-Mat_n1[i])*t1;
                        (*LocalJacobi)(l, 2*LengthStrs-1+j) += -dlmd*Co*t2;
                     }

                     // dr_M/dLambda --------------------------
                     (*LocalJacobi)(l, LocDim-1)
                        += -Co*(supMat[i]-Mat_n1[i])*dW;
                  }

                  // d(r_stress)__dstress --------------------------
                  dG__dStress_dStress(nStress, xi_n1, Stress_Inv, Mat_n1, LengthStrs);
                  for(i=0; i<LengthStrs; i++)
                  {
                     l = i;
                     for(j=0; j<LengthStrs; j++)
                     {
                        if(i==j) (*LocalJacobi)(l,j) +=1.0;
                        t1=0.0;
                        t2=0.0;
                        for(k=0; k<LengthStrs; k++)
                        {
                           c2=1.0;
                           if(k>2)
                           {
                              if(j>2) c2=4.0;
                              else c2=2.0;
                           }
                           t1 += c2*(*De)(i,k)*(*d2G_dSdS)(k,j);
                        }
                        // dG_dStress_dStress
                        (*LocalJacobi)(i,j) += dlmd*t1;

                     }
                  }

                  // d(r_M)__dstress --------------------------
                  for(i=0; i<LengthMat; i++)
                  {
                     Co = Ch;
                     if(i>4) Co = Cd;
                     l = i+2*LengthStrs-1;
                     // df3_dStress --------------------------
                     for(j=0; j<LengthStrs; j++)
                     {
                        c1=1.0;
                        if(j>2) c1=2.0;
                        t1=0.0;
                        for(k=0; k<LengthStrs; k++)
                        {
                           c2=1.0;
                           if(k>2) c2=2.0;
                           t1 += c1*c2*Stress_n1[k]*(*d2G_dSdS)(k,j);
                        }
                        t2 = -dlmd*Co*(supMat[i]-Mat_n1[i])*(c1*dG_dS[j]+t1);
                        (*LocalJacobi)(l,j) += t2;
                     }
                  }

                  // d(r_xi)/d... --------------------------
                  // d2G_dSdS: df2_dS. d2G_dSdM: df2_dXi.
                  dfun2(nStress, xi_n1, Stress_Inv, Mat_n1, LengthStrs);
                  for(i=0; i<LengthStrs; i++)
                  {
                     if(i<2)  l = i+LengthStrs;
                     if(i==2) continue;
                     if(i>2) l = i-1+LengthStrs;

                     for(j=0; j<LengthStrs; j++)
                     {
                        c1=1.0;
                        if(j>2) c1=2.0;
                        (*LocalJacobi)(l,j)
                           += c1*dlmd*br*(*d2G_dSdS)(i,j)/psi1;
                        if(j==2)                  //x_ii*delta_ii = 0;
                        {
                           (*LocalJacobi)(l, LengthStrs)
                              -= c1*dlmd*br*(*d2G_dSdM)(i,j)/psi1;
                           (*LocalJacobi)(l, LengthStrs+1)
                              -= c1*dlmd*br*(*d2G_dSdM)(i,j)/psi1;
                        }
                        if(j<2)
                        {
                           (*LocalJacobi)(l, LengthStrs+j)
                              += c1*dlmd*br*(*d2G_dSdM)(i,j)/psi1;
                        }
                        if(j>2)
                        {
                           (*LocalJacobi)(l, LengthStrs+j)
                              += c1*dlmd*br*(*d2G_dSdM)(i,j-1)/psi1;
                        }

                        if(i==j)                  //dxi/dxi
                        {
                           if(j<2) (*LocalJacobi)(l, LengthStrs+j) += 1.0;
                           if(j==2) continue;
                           if(j>2) (*LocalJacobi)(l, LengthStrs+j-1) += 1.0;
                        }

                     }
                     // d(r_xi)/dMat
                     t1 = sqrt(Stress_Inv[1]/3.0);
                     (*LocalJacobi)(l, 2*LengthStrs-1)
                        -= 0.25*dlmd*br*t1*I_n1
                        *(I_n1*xi_n1[i]-mr*nStress[i])/(Stress_Inv[6]*psi1*I_n1);
                     (*LocalJacobi)(l, 2*LengthStrs+1)
                        -=      dlmd*br*t1*Mat_n1[2]*I_p4
                        *(I_n1*xi_n1[i]-mr*nStress[i])/(Stress_Inv[6]*psi1*I_n1);
                     // d(r_xi)/dLambda
                     (*LocalJacobi)(l, LocDim-1)
                        -=  br*t1*(I_n1*xi_n1[i]-mr*nStress[i])/(PSI*psi1*I_n1);
                  }

                  // d(r_F)/d...  --------------------------
                  dF_dNStress(dF_dS, nStress, Stress_Inv, Mat_n1, LengthStrs);
                  dF_dMat(dF_dM, Stress_Inv, Mat_n1);
                  l= LocDim-1;
                  for(j=0; j<LengthStrs; j++)
                  {
                     c1=1.0;
                     if(j>2) c1=2.0;
                     //dr_F/dxi
                     if(br>0.0)
                     {
                        if(j==2)
                        {
                           (*LocalJacobi)(l,LengthStrs) += c1*I_n1*dF_dS[j]/3.0;
                           (*LocalJacobi)(l, LengthStrs+1) += c1*I_n1*dF_dS[j]/3.0;
                        }
                        else if (j<2)
                           (*LocalJacobi)(l, LengthStrs+j) += -c1*I_n1*dF_dS[j]/3.0;
                        else if (j>2)
                           (*LocalJacobi)(l, LengthStrs+j-1) += -c1*I_n1*dF_dS[j]/3.0;
                     }
                  }

                  // dr_F/dStress
                  dF_dStress(dF_dS, xi_n1, Stress_Inv, Mat_n1, LengthStrs);
                  for(j=0; j<LengthStrs; j++)
                  {
                     c1=1.0;
                     if(j>2) c1=2.0;
                     (*LocalJacobi)(l,j) += c1*dF_dS[j];
                  }

                  // dr_F/dM
                  if((Ch>0.0)&&(Cd>0.0))
                     for(j=0; j<LengthMat; j++)
                  {
                     (*LocalJacobi)(l, 2*LengthStrs-1+j) += dF_dM[j];
                  }
                  //-------- End Local Jacibin ----------

                  //-------- RHS ------------
                  for(i=0; i<LengthStrs; i++)
                  {
                     rhs_l[i] = Stress_n1[i];
                     for(k=0; k<LengthStrs; k++)
                     {
                        c2=1.0;
                        if(k>2) c2=2.0;
                        rhs_l[i] += c2*dlmd*(*De)(i,k)*dG_dS[k];
                     }
                     rhs_l[i] -= Stress_n[i] + factor*dStress[i];

                     //For Xi
                     if(i==2) continue;
                     if(i<2) l = LengthStrs+i;
                     else if(i>2) l = LengthStrs+i-1;
                     rhs_l[l]
                        = xi_n1[i]-xi_n[i]
                        -dlmd*br*sqrt(Stress_Inv[1]/3.0)
                        *(I_n1*xi_n1[i]-mr*nStress[i])/(PSI*psi1*I_n1);
                  }
                  for(i=0; i<LengthMat; i++)
                  {
                     Co = Ch;
                     if(i>4) Co = Cd;
                     rhs_l[2*LengthStrs-1+i]
                        = Mat_n1[i]-Mat_n[i]
                        -Co*dlmd*(supMat[i]-Mat_n1[i])*dW;
                  }
                  // The yield function
                  rhs_l[LocDim-1] = F;
                  //-------- End RHS -----------------

                  //--------- Compute the error of the residual  --------
                  if(NPStep==1)
                  {
                     NormR0 = 0.0;
                     for(i=0; i<LocDim; i++)
                     {
                        NormR0 += rhs_l[i]*rhs_l[i];
                     }
                  }

                  if(sqrt(NormR0)<TolF) ErrLoc = 0.01*TolF;
                  else
                  {
                     NormR1 = 0.0;
                     for(i=0; i<LocDim; i++)
                     {
                        NormR1 += rhs_l[i]*rhs_l[i];
                     }
                     ErrLoc = sqrt(NormR1/NormR0);
                  }
                  if(fabs(F)<TolF||F<0.0) ErrLoc = 0.01*TolF;
                  if(NormR1<TolF) ErrLoc = 0.01*TolF;
                  //------ End: Compute the error of the residual

                  damping=1.0;
                  if(NPStep>1)
                  {
                     if((ErrLoc/ErrLoc0)>0.01) damping=0.5;
                  }
                  if(NPStep>MaxIter-1) damping=0.2;

                  //------  Solve the linear equation
                  Gauss_Elimination(LocDim, *LocalJacobi, Li, x_l);
                  Gauss_Back(LocDim, *LocalJacobi, rhs_l, Li, x_l);
                  //------  End Solve the linear equation

                  //------ Compute the error of the solution
                  /*
                  if(NPStep==1)
                  {
                     NormU0 = 0.0;
                     for(i=0; i<LocDim; i++)
                     {
                        NormU0 += x_l[i]*x_l[i];
                     }
                  }

                  if(sqrt(NormU0)<TolLocalNewT) ErrLoc = 0.1*TolLocalNewT;
                  else
                  {
                  NormU1 = 0.0;
                  for(i=0; i<LocDim; i++)
                  {
                  NormU1 += x_l[i]*x_l[i];
                  }
                  ErrLoc = min(ErrLoc, sqrt(NormU1/NormU0));
                  }
                  if(fabs(F)<TolLocalNewT) ErrLoc = 0.1*TolLocalNewT;
                  */
                  //----- End: Compute the error of the solution
                  ErrLoc0 = ErrLoc;

                  //----- Update the Newton-Raphson step
                  for(i=0; i<LocDim; i++)
                     x_l[i] *= damping;

                  for(i=0; i<LengthStrs; i++)
                  {
                     Stress_n1[i] -= x_l[i];
                     nStress[i] = Stress_n1[i];
                     if(i==2) continue;
                     else if(i<2)  xi_n1[i] -= x_l[i+LengthStrs];
                     else if(i>2)  xi_n1[i] -= x_l[i-1+LengthStrs];
                  }
                  I_n1 = DeviatoricStress(nStress);
                  xi_n1[2] = -xi_n1[0]-xi_n1[1];

                  for(i=0; i<LengthStrs; i++)
                     nStress[i] -= I_n1*xi_n1[i]/3.0;

                  for(i=0; i<LengthMat; i++)
                     Mat_n1[i] -= x_l[i+2*LengthStrs-1];

                  dlmd -= x_l[LocDim-1];

               }                                  // End if (F>0.0)

               if(NPStep>MaxIter)
               {
                  printf("\n ~O~ ~O~  Local Newton-Raphson has problem in convergence in MatCalcStressWM().. \n");
                  printf("\n F = %g\n", F);
                  break;                          //abort();
               }

            }                                     // End of local Newton-Raphson

            // Accumulated plastic strain
            if(subPLASTIC>0) ep += dlmd*sqrt(2.0*TensorMutiplication2(dG_dS, dG_dS, dim)/3.0);

            //-------  Compute the consistent tangential matrix
            if(Update<=0&&subPLASTIC>0)
            {

               //--- 1.  Compute the inverse of the Jacobian -
               for(j=0; j<LocDim-1; j++)
               {
                  for(i=0; i<LocDim; i++)
                  {
                     rhs_l[i] = 0.0;
                     if(i==j) rhs_l[i] = 1.0;
                  }
                  // the i_th column of the invJac matrix
                  Gauss_Back(LocDim, *LocalJacobi, rhs_l, Li, x_l);
                  for(i=0; i<LocDim-1; i++)
                  {
                     (*inv_Jac)(i, j) = x_l[i];
                  }
               }

               //- 2.  A*A*A*... -
               if(subPLASTIC==1)                  //First substep, which has plasticity
               {
                  for(i=0; i<LocDim-1; i++)
                  {
                     for(j=0; j<LengthStrs; j++)
                        (*sumA_Matrix)(i, j) = (*inv_Jac)(i, j)*factor;
                  }
               }
               else
               {
                  for(i=0; i<LocDim-1; i++)
                  {
                     for(j=0; j<LengthStrs; j++)
                     {
                        if(i==j)
                           (*sumA_Matrix)(i,j) += factor;

                     }
                  }

                  LocalJacobi->LimitSize(LocDim-1,LengthStrs);
                  (*LocalJacobi) = 0.0;
                  inv_Jac->multi(*sumA_Matrix, *LocalJacobi);
                  for(i=0; i<LocDim-1; i++)
                  {
                     for(j=0; j<LengthStrs; j++)
                     {
                        (*sumA_Matrix)(i,j) = (*LocalJacobi)(i,j);
                     }
                  }
                  LocalJacobi->LimitSize(LocDim, LocDim);

               }
               //- 3.  D_ep -
               if(iSub==preSub-1)
               {
                  for(i=0; i<LengthStrs; i++)
                  {
                     for(j=0; j<LengthStrs; j++)
                     {
                        (*D_ep)(i,j) = 0.0;
                        for(k=0; k<LengthStrs; k++)
                        {
                           (*D_ep)(i,j)
                              += (*sumA_Matrix)(i,k)*(*De)(k,j);
                        }
                     }
                  }
               }
            }
            ///  End Compute the consistent tangential matrix

            // Update the substep
            for(i=0; i<LengthStrs; i++)
            {
               Stress_n[i] = Stress_n1[i];
               xi_n[i] = xi_n1[i];
            }

            for(i=0; i<LengthStrs; i++)
               Mat_n[i] = Mat_n1[i];

         }                                        // End of Compute stresses by substepping

      }                                           // If F>0.0

      // Save the current stresses
      if(Update>0)
      {

         (*ele_val->pStrain)(GPiGPj) = ep;
         //for(i=0; i<LengthStrs; i++)
         //    (*ele_val->Stress)(i, GPiGPj) = Stress_n1[i];
         for(i=0; i<LengthStrs-1; i++)
            (*ele_val->xi)(i, GPiGPj) = xi_n1[i];
         for(i=0; i<LengthMat; i++)
            (*ele_val->MatP)(i, GPiGPj) = Mat_n1[i];

      }
      //   else
      //   {  // New stresses passed through dStress for the residual computation
      for(i=0; i<LengthStrs; i++)
         dStress[i] = Stress_n1[i];
      //   }

      /*
      // Already considered in element_dm, ele_val_dm
      // For the case of the initial stress being given, the contribution of
      // the initial stress to the right hand side should be taken into account
      dStress[0] -= (*data_Plasticity)(20); // Initial stress_xx
      dStress[1] -= (*data_Plasticity)(21); // Initial stress_yy
      dStress[2] -= (*data_Plasticity)(22); // Initial stress_zz
      */

      return PLASTIC;
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::dF_dNStress

     Aufgabe:
      Computing the derivatives of the yield function with respect to
      the normalized stresses (Weimar's model)

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
          double *dFdS               : The derivatives of the yield function
   with respect to stresses
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Modification for the 3D case
   08/2004   WW  Set as a member of CSolidProperties
   ***************************************************************************/
   void CSolidProperties::dF_dNStress(double *dFdS, const double *DevS,
      const double *S_Invariants, const double *MatN1, const int LengthStrs)
   {
      int i;
      double In1 = S_Invariants[0];
      double I_p3 = In1*In1*In1;
      double ang = S_Invariants[3];
      double PHI = S_Invariants[4];

      /* 1. Compute the derivatives of the yield function with respect to stresses */
      /* 1.1 s_{jk}s_{ki} */
      if(LengthStrs==4)                           // 2D
      {
                                                  // s11*s11+s12*s21
         dFdS[0] =  DevS[0]*DevS[0]+DevS[3]*DevS[3];
                                                  // s21*s12+s22*s22
         dFdS[1] =  DevS[3]*DevS[3]+DevS[1]*DevS[1];
         dFdS[2] =  DevS[2]*DevS[2];              // s33*s33
                                                  // s11*s12+s12*s22
         dFdS[3] =  DevS[0]*DevS[3]+DevS[3]*DevS[1];
      }
      else
      {
                                                  // s11*s11+s12*s21+s13*s31
         dFdS[0] =  DevS[0]*DevS[0]+DevS[3]*DevS[3]+DevS[4]*DevS[4];
                                                  // s21*s12+s22*s22+s23*s32
         dFdS[1] =  DevS[3]*DevS[3]+DevS[1]*DevS[1]+DevS[5]*DevS[5];
                                                  // s31*s13+s32*s23+s33*s33
         dFdS[2] =  DevS[4]*DevS[4]+DevS[5]*DevS[5]+DevS[2]*DevS[2];
                                                  // s11*s12+s12*s22+s13*s32
         dFdS[3] =  DevS[0]*DevS[3]+DevS[3]*DevS[1]+DevS[4]*DevS[5];
                                                  // s11*s13+s12*s23+s13*s33
         dFdS[4] =  DevS[0]*DevS[4]+DevS[3]*DevS[5]+DevS[4]*DevS[2];
                                                  // s21*s13+s22*s23+s23*s33
         dFdS[5] =  DevS[3]*DevS[4]+DevS[1]*DevS[5]+DevS[5]*DevS[2];
      }
      /* 1.2 dtheta/ds*/
      for(i=0; i<LengthStrs; i++)
      {
         dFdS[i] -= 1.5*S_Invariants[2]*DevS[i]/S_Invariants[1];
         if(i<3)  dFdS[i] -= 2.0*S_Invariants[1]/3.0;
         dFdS[i] /= pow(S_Invariants[1],1.5);
      }
      /* 1.3 dF/ds   ..  */
      for(i=0; i<LengthStrs; i++)
      {
         dFdS[i] *= MatN1[5]*MatN1[6]*S_Invariants[1]*pow(1+MatN1[5]*ang, MatN1[6]-1.0);
         dFdS[i] += pow(1+MatN1[5]*ang, MatN1[6])*DevS[i];
         dFdS[i] /= (2*PHI);
         if(i<3)
         {
            dFdS[i] += (MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)/(2*PHI);
            dFdS[i] += MatN1[1]+2.0*MatN1[3]*In1;
         }
      }
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::dF_dStress

     Aufgabe:
      Computing the derivatives of the yield function with respect to
      stresses (Weimar's model)
      = dF_dNStress- dF_dStress;

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
   double *dFdS               : The derivatives of the yield function
   with respect to stresses
   const double *RotV         : Rotational variables
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Modification for the 3D case
   08/2004   WW  Set as a member of CSolidProperties
   ***************************************************************************/
   void CSolidProperties::dF_dStress(double *dFdS,  const double *RotV,
      const double *S_Invariants, const double *MatN1, const int LengthStrs)
   {
      int i;
      double In1 = S_Invariants[0];
      double I_p3 = In1*In1*In1;
      double PHI = S_Invariants[4];

      for(i=0; i<LengthStrs; i++)
      {
         dFdS[i] -= ((MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)/(2*PHI)
            + MatN1[1]+2.0*MatN1[3]*In1)*RotV[i];
      }
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::dF_dMat

     Aufgabe:
      Computing the derivatives of the yield function with respect to
      material parameters (Weimar's model)

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
          double *dFdM               : The derivatives of the yield function
   with respect to material parameters
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   08/2004   WW  Set as a member of CSolidProperties
   ***************************************************************************/
   void CSolidProperties::dF_dMat(double *dFdM, const double *S_Invariants, const double *MatN1)
   {
      double In1 = S_Invariants[0];
      double I_p2 = In1*In1;
      double I_p4 = I_p2*I_p2;
      double ang = S_Invariants[3];
      double PHI = S_Invariants[4];
      /*dF_dAlpha*/
      dFdM[0] = 0.25*I_p2/PHI;
      /*dF_dBeta*/
      dFdM[1] = In1;
      /*dF_dDelta*/
      dFdM[2] = MatN1[2]*I_p4/PHI;
      /*dF_dEpsilon*/
      dFdM[3] = I_p2;
      /*dF_dKappa*/
      dFdM[4] = -1.0;
      /*dF_dGamma*/
      dFdM[5] = 0.5*S_Invariants[1]*MatN1[6]*ang*pow(1.0+MatN1[5]*ang, MatN1[6]-1.0)/PHI;
      /*dF_dM*/
      dFdM[6] = 0.5*S_Invariants[1]*pow(1.0+MatN1[5]*ang, MatN1[6])
         *log(1.0+MatN1[5]*ang)/PHI;

   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::dG_dNStress

     Aufgabe:
      Computing the derivatives of the plastic potential with respect to
      normal stresses (Weimar's model)

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
          double *dGdS               : The derivatives of the plastic potential
   with respect to stresses
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Modification for the 3D case
   08/2004   WW  Set as a member of CSolidProperties
   ***************************************************************************/
   void CSolidProperties::dG_dNStress(double *dGdS, const double *DevS,
      const double *S_Invariants, const double *MatN1, const int LengthStrs)
   {
      int i;
      const double psi1= (*data_Plasticity)(14);
      const double psi2= (*data_Plasticity)(15);
      double In1 = S_Invariants[0];
      double I_p3 = In1*In1*In1;
      double PSI = S_Invariants[5];

      for(i=0; i<LengthStrs; i++)
      {
         dGdS[i] = 0.5*DevS[i]/PSI/psi1;
         if(i<3)
            dGdS[i] += 0.5*(MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)/PSI
               +(1.0+psi2)*MatN1[1]+2.0*MatN1[3]*In1;
      }
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::dG_dNStress_dNStress

     Aufgabe:
      Computing the second derivatives of the plastic potential with respect to
      normal stresses (Weimar's model)

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
          double *dG_dSdS            : The second derivatives of the plastic potential
   with respect to stresses
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Extended to 3D
   08/2004   WW  Set as a member of CSolidProperties
   ***************************************************************************/
   void CSolidProperties::dG__dNStress_dNStress(const double *DevS,
      const double *S_Invariants, const double *MatN1, const int LengthStrs )
   {
      int i,j;
      int ii, jj, kk, ll;
      double delta_ij_kl;
      const double psi1= (*data_Plasticity)(14);
      double In1 = S_Invariants[0];
      double I_p2 = In1*In1;
      double I_p3 = In1*In1*In1;
      double PSI = S_Invariants[5];
      double PSI_p3 = S_Invariants[6];

      for(i=0; i<LengthStrs; i++)
      {
         ii=i;
         jj=i;
         if(i==3)                                 //s_12
         {
            ii=1;
            jj=2;
         }
         //For 3D
         else if(i==4)                            //s_13
         {
            ii=1;
            jj=3;
         }
         else if(i==5)                            //s_23
         {
            ii=2;
            jj=3;
         }

         for(j=0; j<LengthStrs; j++)
         {
            kk=j;
            ll=j;
            if(j==3)
            {
               kk=1;
               ll=2;
            }
            //For 3D
            else if(j==4)                         //s_13
            {
               kk=1;
               ll=3;
            }
            else if(j==5)                         //s_23
            {
               kk=2;
               ll=3;
            }

            delta_ij_kl = Kronecker(ii,jj)*Kronecker(kk,ll);
            //dG_dSdS[i*LengthStrs+j] =
            (*d2G_dSdS)(i,j) =
               0.5*((Kronecker(ii,kk)*Kronecker(jj,ll)-delta_ij_kl/3.0)/psi1
               +(MatN1[0]+12.0*MatN1[2]*MatN1[2]*I_p2)*delta_ij_kl)/PSI
               -0.25*(DevS[i]/psi1+(MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)*Kronecker(ii,jj))
               *(DevS[j]/psi1+(MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)*Kronecker(kk,ll))/PSI_p3
               +2.0*MatN1[3]*delta_ij_kl;
         }

      }
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::dG_dStress_dStress

     Aufgabe:
      Computing the second derivatives of the plastic potential with respect to
      stresses (Weimar's model)

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
          double *dG_dSdS            : The second derivatives of the plastic potential
   with respect to stresses
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *RotV         : Rotational variables
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Extended to 3D
   08/2004   WW  Set as a member of CSolidProperties
   ***************************************************************************/
   void CSolidProperties::dG__dStress_dStress(const double *DevS,  const double *RotV,
      const double *S_Invariants, const double *MatN1, const int LengthStrs)
   {
      int i,j;
      int ii, jj;
      const double psi1= (*data_Plasticity)(14);
      double In1 = S_Invariants[0];
      double I_p2 = In1*In1;
      double I_p3 = In1*In1*In1;
      double PSI = S_Invariants[5];
      double PSI_p3 = S_Invariants[6];

      for(i=0; i<LengthStrs; i++)
      {
         ii=i;
         jj=i;
         if(i==3)
         {
            ii=1;
            jj=2;
         }
         //For 3D
         else if(i==4)                            //s_13
         {
            ii=1;
            jj=3;
         }
         else if(i==5)                            //s_23
         {
            ii=2;
            jj=3;
         }

         for(j=0; j<LengthStrs; j++)
         {
            (*d2G_dSdS)(i,j) -=
               0.5*(MatN1[0]+12.0*MatN1[2]*MatN1[2]*I_p2)*Kronecker(ii,jj)*RotV[j]/PSI
               -0.25*(DevS[i]/psi1+(MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)*Kronecker(ii,jj))
               *(MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)*RotV[j] /PSI_p3
               +2.0*MatN1[3]*Kronecker(ii,jj)*RotV[j];
         }

      }
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::dG_dSTress_dMat

     Aufgabe:
      Computing the second derivatives of the plastic potential with respect to
      stresses and with recpect to material parameters (Weimar's model)

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
          double *dG_dSdM            : The second derivatives of the plastic potential
   with respect to stresses and material parameters

   2D:   dG_dS11dM1 dG_dS11dM2 dG_dS11dM3 dG_dS11dM4
   dG_dS22dM1 dG_dS22dM2 dG_dS22dM3 dG_dS22dM4
   dG_dS12dM1 dG_dS12dM2 dG_dS12dM3 dG_dS12dM4
   dG_dS33dM1 dG_dS33dM2 dG_dS33dM3 dG_dS33dM4

   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Extended to 3D
   08/2004   WW  Set as a member of CSolidProperties
   ***************************************************************************/
   void CSolidProperties::dG_dSTress_dMat(const double *DevS,
      const double *S_Invariants, const double *MatN1, const int LengthStrs)
   {
      int i;
      const double psi1= (*data_Plasticity)(14);
      const double psi2= (*data_Plasticity)(15);
      double In1 = S_Invariants[0];
      double I_p2 = In1*In1;
      double I_p3 = In1*I_p2;
      double I_p4 = In1*I_p3;
      double PSI = S_Invariants[5];
      double PSI_p3 = S_Invariants[6];
      int ii, jj;
      double delta_ij;

      for(i=0; i<LengthStrs; i++)
      {
         ii=i;
         jj=i;
         if(i==3)
         {
            ii=1;
            jj=2;
         }
         //For 3D
         else if(i==4)                            //s_13
         {
            ii=1;
            jj=3;
         }
         else if(i==5)                            //s_23
         {
            ii=2;
            jj=3;
         }

         delta_ij = Kronecker(ii,jj);
         //dG_dSdM[i*LengthStrs]   //dG_dS_dAlpha
         (*d2G_dSdM)(i,0)                         //dG_dS_dAlpha
            =   0.5*In1*delta_ij/PSI
            -0.125*(DevS[i]/psi1+(MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)*delta_ij)
            *I_p2/PSI_p3;
         //dG_dSdM[i*LengthStrs+1]   //dG_dS_dBeta
         (*d2G_dSdM)(i,1)                         //dG_dS_dBeta
            =  (1+psi2)*delta_ij;
         //dG_dSdM[i*LengthStrs+2]   // dG_dS_dDelta
         (*d2G_dSdM)(i,2)                         // dG_dS_dDelta
            =   4.0*MatN1[2]*I_p3*delta_ij/PSI
            -0.5*(DevS[i]/psi1+(MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)*delta_ij)
            *MatN1[2]*I_p4/PSI_p3;
         //dG_dSdM[i*LengthStrs+3]   // dG_dS_dEpsilon
         (*d2G_dSdM)(i,3)                         // dG_dS_dEpsilon
            =   2.0*In1*delta_ij;
      }
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::dfun2

     Aufgabe:
      Computing (Weimar's model)
       1:  the derivatives of the rotational variable equation with
            respect to stresses
       2:  the derivatives of the rotational variable equation with
            respect to rotational variables

      Note : xi_11+xi_22+xi_33=0.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *RotV         : The hardening rotational variables
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:
   double *d2G_dSdS              : Results 1
   double *d2G_dSdM              : Results 2

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Extended to 3D
   08/2004   WW  Set as a member of CSolidProperties
   ***************************************************************************/
   void CSolidProperties::dfun2(const double *DevS,
      const double *RotV, const double *S_Invariants,
      const double *MatN1, const int LengthStrs)
   {
      int i,j, l;
      int ii, jj, kk, ll;
      const double mr= (*data_Plasticity)(19);
      const double psi1= (*data_Plasticity)(14);

      double In1 = S_Invariants[0];
      double I_p3 = In1*In1*In1;
      double PSI = S_Invariants[5];
      double PSI_p3 = S_Invariants[6];
      double var = sqrt(S_Invariants[1]/3.0);
      l = 0;
      for(i=0; i<LengthStrs; i++)
      {
         // if(i<2) l = i*LengthStrs;
         if(i==2) continue;                       //
         // if(i>2) l = (i-1)*LengthStrs;

         ii=i;
         jj=i;
         if(i==3)
         {
            ii=1;
            jj=2;
         }
         //For 3D
         else if(i==4)                            //s_13
         {
            ii=1;
            jj=3;
         }
         else if(i==5)                            //s_23
         {
            ii=2;
            jj=3;
         }

         for(j=0; j<LengthStrs; j++)
         {
            kk=j;
            ll=j;
            if(j==3)
            {
               kk=1;
               ll=2;
            }
            //For 3D
            else if(j==4)                         //s_13
            {
               kk=1;
               ll=3;
            }
            else if(j==5)                         //s_23
            {
               kk=2;
               ll=3;
            }

            // Derivative with respect to normal stresses
            (*d2G_dSdS)(l,j) = -(
               (In1*RotV[i]-mr*DevS[i])*DevS[j]/(6.0*var*PSI)
               -0.5*var*(In1*RotV[i]-mr*DevS[i])
               *(DevS[j]/psi1+(MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)*Kronecker(kk,ll))/PSI_p3
               +var*(RotV[i]*Kronecker(kk,ll)
               -mr*(Kronecker(ii,kk)*Kronecker(jj,ll)-Kronecker(ii,jj)*Kronecker(kk,ll)/3.0))/PSI);

            (*d2G_dSdS)(l,j) /= In1;
            (*d2G_dSdS)(l,j) -= 0.5*(In1*RotV[i]-mr*DevS[i])*Kronecker(kk,ll)/(PSI*In1*In1);

            if(j<4)
               (*d2G_dSdM)(l,j) = - In1*(*d2G_dSdS)(i,j)/3.0
                  - var*In1*Kronecker(ii,kk)*Kronecker(jj,ll)/(PSI*In1);

            // Derivative with respect to stresses
            (*d2G_dSdS)(l,j) +=
               -0.5*var*(In1*RotV[i]-mr*DevS[i])
               *(MatN1[0]*In1+4.0*MatN1[2]*MatN1[2]*I_p3)*RotV[j]/PSI_p3+var*RotV[i]*RotV[j]/PSI
               +0.5*(In1*RotV[i]-mr*DevS[i])*RotV[j]/(PSI*In1*In1);
         }

      }
   }

   /**************************************************************************
    ROCKFLOW - Funktion: CSolidProperties::Gauss_Elimination and Gauss_Back

     Aufgabe: Mini linear equation solver by Gasssian elemination with
              scaled partial pivoting.

    Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
          const int DimE             : Dimension of the equations
          double *AA                 : Matrix
   double * rhs               : Right hand side
   int *L                     : Teporarily used array
   double *xx                 : Solution
   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version
   08/2004   WW  Set as a member of CSolidProperties
   ***************************************************************************/
   void CSolidProperties::Gauss_Elimination(const int DimE, Matrix& AA, int *L,  double *xx)
   {
      int i,j, k, jj, lk;
      double var, R;

      for(i=0; i<DimE; i++)
      {
         L[i] = i;
         var = 0.0;
         for(j=0; j<DimE; j++)
         {
            if(fabs(AA(i,j))>var) var = fabs(AA(i,j));
            L[i] = i;
         }
         xx[i]=var;
      }

      for(k=0; k<DimE-1; k++)
      {
         var = 0.0;
         jj=0;
         for(i=k; i<DimE; i++)
         {
            R = fabs(AA(L[i], k)/xx[L[i]]);

            if(R>var)
            {
               jj = i;
               var = R;
            }
         }
         lk = L[jj];
         L[jj] = L[k];
         L[k] = lk;

         for(i=k+1; i<DimE; i++)
         {
            var = AA(L[i], k)/AA(lk, k);

            for(j=k+1; j<DimE; j++)
            {
               AA(L[i],j) -= var*AA(lk, j);
            }
            AA(L[i], k) = var;
         }
      }
   }

   void CSolidProperties::Gauss_Back(const int DimE, Matrix& AA, double * rhs, int *L, double *xx)
   {
      int i,j, k;
      double var;

      /* Back substituting */
      for(k=0; k<DimE-1; k++)
      {
         for(i=k+1; i<DimE; i++)
         {
            rhs[L[i]] -= AA(L[i], k)*rhs[L[k]];
         }
      }
      xx[DimE-1] = rhs[L[DimE-1]]/AA(L[DimE-1], DimE-1);
      for(i=DimE-2; i>=0; i--)
      {
         var = rhs[L[i]];
         for(j=i+1; j<DimE; j++)
         {
            var -= AA(L[i], j)*xx[j];
         }
         xx[i] = var/AA(L[i], i);
      }
   }

   /*=========================================================================

                       CAM-CLAY Model
   =========================================================================*/
   /**************************************************************************
      ROCKFLOW - Funktion: CSolidProperties::CalStress_and_TangentialMatrix_CC

     Aufgabe:
      Computing the stresses at Gauss point and return the plastical status of this
      point. Plastic model: Cam-Clay.
      (Explicit integration of elastic model)

     Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
   double *TryStress                   :   Incremental straines as input
   New stresses as output
   double *Dep                         :   The consist tangential matrix
   const int GPiGPj                    :   Local indeces of Gauss Points
   const int Update                    :   Indicator. If true, update the gauss point
   values

   Ergebnis:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   11/2003   WW  Erste Version

   **************************************************************************/
   //#define New
#define associative
   void CSolidProperties::CalStress_and_TangentialMatrix_CC(const int GPiGPj,
      const ElementValue_DM *ele_val, double *dStrain,  Matrix *Dep, const int Update)
   {
      int i, ns;
      double p, q, p_tr, q_tr, p_c, p_cn;
      double ep;
      double e_0, dev;
      double F0, F, vep, vep0, TolP;
      double var1, alpha1, alpha2, beta1, beta2;
      double gamma1, gamma2,gamma3,gamma4,gamma5;
      double dfdp, dfdq, dampFac;
      const double fac = sqrt(2.0/3.0);
#ifdef New
      double J11,J12, J21, J22;
#else
      double  alpha3, Jac;
#endif
      static double DevStress[6], TryStress[6];

      const int MaxI = 400;
      int NPStep;

      const double SwellIndex = (*data_Plasticity)(2);
      const double CompressIndex = (*data_Plasticity)(1);
      const double M_s = (*data_Plasticity)(0);
      const double M2 = M_s*M_s;
      const double pmin = (*data_Plasticity)(9);

      double vartheta=0.0;

      TolP = Tolerance_Local_Newton*1.0e4;

      int dim = 2;
      ns =4;
      if(Dep->Cols()>4)
      {
         dim = 3;
         ns = 6;
      }
      bool isLoop = true;                         // Used only to avoid warnings with .net

      // Get the total effective plastic strain
      ep = (*ele_val->pStrain)(GPiGPj);

      // Get the preconsolidation pressure of the previous time step
      p_cn = (*ele_val->prep0)(GPiGPj);

      // Get the void ratio of the previous time step
      e_0 = (*ele_val->e_i)(GPiGPj);

      p = -(  (*ele_val->Stress)(0,GPiGPj)
         +(*ele_val->Stress)(1,GPiGPj)
         +(*ele_val->Stress)(2,GPiGPj))/3.0;

      if(fabs(p)<pmin)
         p = pmin;

      // Volume strain increment
      dev = -(dStrain[0]+dStrain[1]+dStrain[2]);
                                                  //TEST. Sign?
      vartheta=(1.0+e_0)/(CompressIndex-SwellIndex);

      // if(fabs(dev)<MKleinsteZahl)
      if(SwellingPressureType==3)
      {
         K = (1.0+e_0)*fabs(p)/(*data_Youngs)(6);
         if(K<10.e6) K = 10.e6;
      }
      else
         K = (1.0+e_0)*fabs(p)/SwellIndex;
      // else
      //    K = (1.0+e_0)*fabs(p)*(exp(PoissonRatio*dev/SwellIndex)-1.0)/dev;

      G = 1.5*K*(1-2.0*PoissonRatio)/(1+PoissonRatio);

      Lambda = K-2.0*G/3.0;

      ElasticConsitutive(dim, Dep);

      for(i=0; i<ns; i++)
         TryStress[i] = 0.0;
      Dep->multi(dStrain, TryStress);

      for(i=0; i<ns; i++)
      {
         TryStress[i] += (*ele_val->Stress)(i,GPiGPj);
         DevStress[i] = TryStress[i];
      }

      p_tr = -DeviatoricStress(DevStress)/3.0;

      q_tr = sqrt(1.5*TensorMutiplication2(DevStress, DevStress, dim));

      // If yield, integrate the stress
      F = q_tr*q_tr/M2 + p_tr*(p_tr-p_cn);

      p = p_tr;
      q = q_tr;
      p_c = p_cn;
      vep = 0.0;
      vep0 = 0.0;

      dampFac=1.0;
      NPStep = 0;

      F0 = F;
      if(pcs_deformation==1) F=-1.0;
      if((*data_Plasticity)(3)<MKleinsteZahl)     // p_c0=0
         F=-1.0;
      //TEST CAM-CLAY
      if(p_tr<0)
         F = -1.0;
      if(F>0.0&&!PreLoad)                         // in yield status
      {

         // Local Newton-Raphson procedure to compute the volume plastic strain
         vep = 0.0;

#ifdef associative
#ifdef New
         //Associative flow rule
         while(isLoop)                            // Newton step for the plastic multiplier
         {
            NPStep++;
            if(NPStep>MaxI)
            {
               // printf("\n Too much iteration in Newton step in the integration of Cam-Clay \n");
               break;
               // abort();
            }
            dfdp = 2.0*p-p_c;
            dfdq = 2.0*q/M2;
            gamma1 = 1.0+2.0*vep*K;
            // J11: df/ddphi. J12: df/dpc
            // dp/dphi
            alpha1 = -K*dfdp/gamma1;
            // dq/dphi
            alpha2 = -q/(vep+M2/(G*6.0));
            //
            J11 = alpha1*dfdp+dfdq*alpha2;
            J12 = dfdp*vep*K/gamma1;
            // J21: dg/ddphi. J22: dg/dpc
            gamma2 = p_cn*exp(vartheta*dev*dfdp);
            J21 = gamma2*vartheta*(dfdp+2.0*vep*alpha1);
            J22 = gamma2*vartheta*vep*(2.0*dev*K/gamma1-1.0)-1.0;

            F = q*q/M2 + p*(p-p_c);
            beta1 = gamma2-p_c;                   // g;
            if(fabs(beta1)<TolP)
            {
               //         if(fabs(F)<TolF) break;
               if(fabs(F/F0)<TolP) break;
            }
            beta2 = J11*J22-J21*J12;

            p_c -= (beta1*J11-F*J21)/beta2;
            vep -= (F*J22-beta1*J12)/beta2;

            p = (p_tr+vep*K*p_c)/gamma1;
            K =
               q = q_tr/(1.0+6.0*G*vep/M2);
         }
#else
         //Associative flow rule
         while(isLoop)                            // Newton step for the plastic multiplier
         {
            NPStep++;

            alpha1 = (2.0*p-p_c)/(1.0+(2.0*K+vartheta*p_c)*vep);
            alpha2 = -q/(vep+M2/(G*6.0));
            alpha3 = vartheta*p_c*alpha1;
            alpha1 *= -K;

            dfdp = 2.0*p-p_c;
            dfdq = 2.0*q/M2;

            Jac = alpha1*dfdp+dfdq*alpha2-p*alpha3;

            while(isLoop)                         // Damp
            {

               vep = vep0 - dampFac*F/Jac;
               p_c = 0.0;

               dampFac = 1.0;
               while(isLoop)                      // Newton step for p_c
               {
                  NPStep++;
                  if(NPStep>MaxI)
                  {
                     //    printf("\n Too much iteration in Newton step the integration of Cam-Clay \n");
                     //TEST	abort();
                     //test
                     break;
                  }

                  //
                  alpha1 = vartheta*vep/(1.0+2.0*vep*K);
                  alpha2 = p_cn*exp(alpha1*(2.0*p_tr-p_c));
                  beta1 = -alpha1*alpha2-1.0;     //dG(p_c)
                  alpha2 -= p_c;                  //G(p_c)
                  p_c -= alpha2/beta1;
                  if(fabs(alpha2)<TolP) break;
                  if(p_c<0.0) break;
               }

               p = (p_tr+vep*K*p_c)/(1.0+2.0*vep*K);
               q = q_tr/(1.0+6.0*G*vep/M2);

               F = q*q/M2 + p*(p-p_c);
               if(F>0.0) break;
               if(fabs(F/F0)<TolP) break;
               dampFac = 0.8;
            }
            vep0 = vep;
            if(fabs(F/F0)<TolP) break;
         }
#endif
         // Plastic strain
         //    alpha1 = 6.0*q*q/(M2*M2*(2.0*p-p_c)*(2.0*p-p_c));
         //    ep += fabs(vep)*sqrt(2.0*(1.0/9.0+alpha1)/3.0);
         ep += 3.0*fabs(vep)*q/M2;

#else
         double var2 = 0.0;
         while(1)
         {
            NPStep++;
            if(NPStep>MaxI)
            {
               printf("\n Too much iteration in Newton step in the integration of Cam-Clay \n");
               abort();
            }

            var1 = (p_tr-p)/(p_tr-0.5*p_c);

            dfdp = 2.0*p-p_c;
            dfdq = 2.0*q/M2;
            Jac = -K*dfdp                         // dF/dp *  dp/dv
               -dfdq*q_tr*(K+0.5*p_c*vartheta*var1)
               /(p_tr-0.5*p_c)                    //dF/dq  * dq/dv
               -p*vartheta*p_c;                   // df/dp_c  * dp_c/dv

            //Update
            dampFac=1.0;

            while(1)
            {
               vep = vep0 - dampFac*F/Jac;

               p_c = p_cn*exp(vep*vartheta);
               p = p_tr-K*vep;
               q = q_tr*(p-0.5*p_c)/(p_tr-0.5*p_c);

               F = q*q/M2 + p*(p-p_c);
               if(F>0.0) break;
               //if(fabs(F)/RF0<TolLocalNewT*1.0e-3) break;
               if(fabs(F)<Tolerance_Local_Newton) break;
               dampFac = 0.5;
            }
            vep0 = vep;
            //if(fabs(F)/RF0<TolLocalNewT*1.0e-3) break;
            if(fabs(F)<Tolerance_Local_Newton) break;
         }

         ep += 3.0*fabs(vep)*q/M2;
#endif

         //-------------------------------------------------------------
         // Consistent tangential matrix

         if(Update<1)
         {
#ifdef associative
                                                  //a
            alpha1 = 1.0+2.0*K*vep+p_c*vartheta*vep;
                                                  //a1
            double a1 = (1.0+p_c*vartheta*vep)/alpha1;
            double a2 = -(2.0*p-p_c)/alpha1;      //a2
                                                  //a3
            double a3 = 2.0*p_c*vartheta*vep/alpha1;
                                                  //a4
            double a4 = vartheta*p_c*(2.0*p-p_c)/(K*alpha1);
                                                  //a5
            double a5 = sqrt(1.5)/(1.0+6.0*G*vep/M2);
                                                  //a6
            double a6 =  -3.0*q/(1.0+6.0*G*vep/M2)/M2;

            alpha1 = -4.0*G*q*a6/M2-K*((2.0*a2-a4)*p-a2*p_c);
            double b1 = -K*((a3-2.0*a1)*p+a1*p_c)/alpha1;
            double b2 = 4.0*G*a5*q/(alpha1*M2);

            gamma1 = 2.0*G*q/q_tr;
            gamma2 = K*(a1+a2*b1)-gamma1/3.0;
            gamma3 = -K*a2*b2;
            gamma4 = -2.0*fac*G*a6*b1;
            gamma5 = 2.0*G*fac*(a5+a6*b2)-gamma1;
#else
            dfdp = 2.0*p-p_c;
            dfdq = 2.0*q/M2;
            var1 = (p_tr-p)/(p_tr-0.5*p_c);
            var2 = (p-0.5*p_c)/(p_tr-0.5*p_c);

            alpha1 = q_tr*var1/(p_tr-0.5*p_c);
            alpha2 = sqrt(3.0/2.0)*var2;
            alpha3 = -q_tr*(1.0+0.5*vartheta*p_c*var1/K)/(p_tr-0.5*p_c);

            double beta0 = K*dfdp-K*dfdq*alpha3+vartheta*p*p_c;
            beta1 = K*(dfdq*alpha1+dfdp)/beta0;
            beta2 = 2.0*G*dfdq*alpha2/beta0;

            gamma1 = 2.0*G*q/q_tr;
            gamma2 = K*(1.0-beta1)-gamma1/3.0;
            gamma3 = K*beta2;
            gamma4 = -fac*K*(alpha1+alpha3*beta1);
            gamma5 = 2.0*G*fac*alpha2-gamma1+fac*alpha3*beta2*K;
#endif

            // Normalize the stress
            for(i=0; i<ns; i++)
            {
               TryStress[i] = DevStress[i];
               TryStress[i] /=q_tr*fac ;
            }

            // Row 1
            beta1 = gamma2+gamma4*TryStress[0];
            beta2 = gamma3+gamma5*TryStress[0];
            (*Dep)(0,0) = gamma1 + beta1 +beta2*TryStress[0];
            (*Dep)(0,1) =          beta1 +beta2*TryStress[1];
            (*Dep)(0,2) =          beta1 +beta2*TryStress[2];
            (*Dep)(0,3) =                 beta2*TryStress[3];
            if(dim==3)
            {
               (*Dep)(0,4) =              beta2*TryStress[4];
               (*Dep)(0,5) =              beta2*TryStress[5];
            }

            // Row 2
            beta1 = gamma2+gamma4*TryStress[1];
            beta2 = gamma3+gamma5*TryStress[1];
            (*Dep)(1,0) =         beta1 +beta2*TryStress[0];
            (*Dep)(1,1) = gamma1 +beta1 +beta2*TryStress[1];
            (*Dep)(1,2) =         beta1 +beta2*TryStress[2];
            (*Dep)(1,3) =                beta2*TryStress[3];
            if(dim==3)
            {
               (*Dep)(1,4) =             beta2*TryStress[4];
               (*Dep)(1,5) =             beta2*TryStress[5];
            }

            // Row 3
            beta1 = gamma2+gamma4*TryStress[2];
            beta2 = gamma3+gamma5*TryStress[2];
            (*Dep)(2,0) =          beta1 +beta2*TryStress[0];
            (*Dep)(2,1) =          beta1 +beta2*TryStress[1];
            (*Dep)(2,2) =   gamma1+beta1 +beta2*TryStress[2];
            (*Dep)(2,3) =                 beta2*TryStress[3];
            if(dim==3)
            {
               (*Dep)(2,4) =              beta2*TryStress[4];
               (*Dep)(2,5) =              beta2*TryStress[5];
            }

            // Row 4
            beta1 = gamma4*TryStress[3];
            beta2 = gamma5*TryStress[3];
            (*Dep)(3,0) = beta1 + beta2*TryStress[0];
            (*Dep)(3,1) = beta1 + beta2*TryStress[1];
            (*Dep)(3,2) = beta1 + beta2*TryStress[2];
            (*Dep)(3,3) =  gamma1+beta2*TryStress[3];
            if(dim==3)
            {
               (*Dep)(3,4) =      beta2*TryStress[4];
               (*Dep)(3,5) =      beta2*TryStress[5];
               // End row 4

               // Row 5
               beta1 = gamma4*TryStress[4];
               beta2 = gamma5*TryStress[4];
               (*Dep)(4,0) = beta1 + beta2*TryStress[0];
               (*Dep)(4,1) = beta1 + beta2*TryStress[1];
               (*Dep)(4,2) = beta1 + beta2*TryStress[2];
               (*Dep)(4,3) =         beta2*TryStress[3];
               (*Dep)(4,4) = gamma1+ beta2*TryStress[4];
               (*Dep)(4,5) =         beta2*TryStress[5];

               // Row 6
               beta1 = gamma4*TryStress[5];
               beta2 = gamma5*TryStress[5];
               (*Dep)(5,0) = beta1 + beta2*TryStress[0];
               (*Dep)(5,1) = beta1 + beta2*TryStress[1];
               (*Dep)(5,2) = beta1 + beta2*TryStress[2];
               (*Dep)(5,3) =         beta2*TryStress[3];
               (*Dep)(5,4) =         beta2*TryStress[4];
               (*Dep)(5,5) = gamma1+ beta2*TryStress[5];
            }

         }

         //-------------------------------------------------------------

         // Update stresses
#ifdef associative
         for(i=0; i<ns; i++)
         {
            DevStress[i] /=1.0+6.0*G*vep/M2;
            TryStress[i] = DevStress[i];
         }
#else
         for(i=0; i<ns; i++)
         {
            DevStress[i] /=1.0+2.0*K*vep/(2.0*p+p_c);
            TryStress[i] = DevStress[i];
         }
#endif

         // True stress
         for(i=0; i<3; i++)
            TryStress[i] -= p;
      }
      else
      {
         if(Update<1)
            ElasticConsitutive(dim, Dep);
      }

      for(i=0; i<ns; i++)
         dStrain[i] = TryStress[i];

      // Save the current stresses
      if(Update>0)
      {
         if((*data_Plasticity)(3)<MKleinsteZahl)  // p_c0=0
         {
            p_c = p+q*q/(M2*p);
            (*ele_val->prep0)(GPiGPj) = p_c;
         }
         if(ep>0.0)
         {
            var1 = dev*(1.+e_0);
            e_0 -= var1;
            (*ele_val->pStrain)(GPiGPj) = ep;
            (*ele_val->e_i)(GPiGPj) = e_0;
            (*ele_val->prep0)(GPiGPj) = p_c;
         }
      }

      // For the case of the initial stress being given, the contribution of
      // the initial stress to the right hand side should be taken into account
      // If the initial stress is not accounted into the final stress
      //
      //  for(i=0; i<3; i++)
      //    dStrain[i] -= (*data_Plasticity)(6+i); // Initial stress
      //
   }

   /**************************************************************************
     GeoSys - Function: Integration with substep for CAM-Clay like model
     Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
      E :
        double *TryStress                   :   Incremental straines as input
                                                New stresses as output
        double *Dep                         :   The consist tangential matrix
        const int GPiGPj                    :   Local indeces of Gauss Points
        const int Update                    :   Indicator. If true, update the gauss point
                                                values

   Ergebnis:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   06/2008   WW  Programming
   **************************************************************************/
   void CSolidProperties::CalStress_and_TangentialMatrix_CC_SubStep(const int GPiGPj,
      const ElementValue_DM *ele_val, double *dStrain,  Matrix *Dep, const int Update)
   {
      int i, ns;
      double p, q, p_tr, q_tr, p_c, p_cn;
      double ep;
      double e_0, dev;
      double F0, F, vep, vep0;
      double var1, alpha1, alpha2, beta1, beta2;
      double gamma1, gamma2,gamma3,gamma4,gamma5;
      double dfdp, dfdq, dampFac;
      const double fac = sqrt(2.0/3.0);
      double  alpha3, Jac;
      static double DevStress[6], TryStress[6];
      static double dStress0[6], dStress1[6];

      const int MaxI = 400;
      int NPStep;
      double sub_step = 0.5;
      double sub_step_sum = 0.0;
      //
      const double SwellIndex = (*data_Plasticity)(2);
      const double CompressIndex = (*data_Plasticity)(1);
      const double M_s = (*data_Plasticity)(0);
      const double M2 = M_s*M_s;
      const double pmin = (*data_Plasticity)(9);
      //
      double vartheta=0.0;
      double suc=0.0;
      double dsuc=0.0;

      p = q = q_tr = p_c = vep = 0.;
      NPStep = 1;
      //
      int dim = 2;
      ns =4;
      if(Dep->Cols()>4)
      {
         dim = 3;
         ns = 6;
      }
      bool isLoop = true;                         // Used only to avoid warnings with .net

      // Get the total effective plastic strain
      ep = (*ele_val->pStrain)(GPiGPj);

      // Get the preconsolidation pressure of the previous time step
      p_cn = (*ele_val->prep0)(GPiGPj);

      // Get the void ratio of the previous time step
      e_0 = (*ele_val->e_i)(GPiGPj);
      // Stress of previous load step
      dev = -(dStrain[0]+dStrain[1]+dStrain[2]);
      for(i=0; i<ns; i++)
      {
         TryStress[i] = (*ele_val->Stress)(i,GPiGPj);
         dStress0[i] = dStress1[i] = 0.;
      }
      //TEST
      int test_c = 0;
      // Begin substeps
      //WW bool OK = true; //OK411
      for(;;)                                     //WW while(OK)
      {
         /*
         if(test_c>2)
          test_c = test_c;
         */
         //TEST
         test_c++;
         //if(test_c>20)
         //  test_c = test_c;
         if(test_c>100)
            break;

         for(int kk=0; kk<2; kk++)
         {
            double *dsig = NULL;
            double  de_vsw = 0.;
            if(kk==0)
            {
               dsig = dStress0;
               for(i=0; i<ns; i++)
                  dsig[i] = TryStress[i];
            }
            else
            {
               dsig = dStress1;
               for(i=0; i<ns; i++)
                  dsig[i] = TryStress[i] + dStress0[0];
            }
            p = -( dsig[0]+dsig[1]+dsig[2])/3.0;
            if(fabs(p)<pmin)
               p = pmin;
            // Volume strain increment
            vartheta=(1.0+e_0)/(CompressIndex-SwellIndex);
            // if(fabs(dev)<MKleinsteZahl)
            if(SwellingPressureType==3)
            {
               suc = (*data_Youngs)(7);
               dsuc = (*data_Youngs)(8);
                                                  //0: sign is nagive -
               de_vsw = TEPSwellingParameter(p)*dsuc/(suc+1.0e5);
               for(i=0; i<3; i++)
                  dStrain[i] += de_vsw;
               K = (1.0+e_0)*fabs(p)/(*data_Youngs)(6);
               if(K<10.e6) K = 10.e6;
               vartheta=(1.0+e_0)/(CompressIndex-(*data_Youngs)(6));
            }
            else
               K = (1.0+e_0)*fabs(p)/SwellIndex;
            // else
            //    K = (1.0+e_0)*fabs(p)*(exp(PoissonRatio*dev/SwellIndex)-1.0)/dev;

            G = 1.5*K*(1-2.0*PoissonRatio)/(1+PoissonRatio);
            Lambda = K-2.0*G/3.0;
            ElasticConsitutive(dim, Dep);
            //
            for(i=0; i<ns; i++)
               dsig[i] = 0.0;
            Dep->multi(dStrain, dsig, sub_step);
            // Recover strain increment
            if(SwellingPressureType==3)
            {
               for(i=0; i<3; i++)
                  dStrain[i] -= de_vsw;
            }
         }
         // Stress estimation
         // double norm_ds = 0.;
         // double norm_s = 0.;
         for(i=0; i<ns; i++)
         {
            TryStress[i] += 0.5*(dStress0[i]+dStress1[i]);
                                                  //DevStress as temporary buffer
            DevStress[i] = dStress0[i]-dStress1[i];
         }
         double norms = StressNorm(TryStress, dim);
         double q_f = 1.0;
         double R_n = 1.0;
         if(norms>DBL_EPSILON)
         {
            R_n = 0.5*StressNorm(DevStress, dim)/norms;
            if(R_n>DBL_EPSILON)
               q_f = 0.95*sqrt(s_tol/R_n);
         }
         else
            R_n = 0.0;
         if(q_f<0.2) q_f = 0.2;
         if(q_f>1.2) q_f = 1.2;
         sub_step *= q_f;
         if((sub_step_sum+sub_step)>1.0)
            sub_step = 1.0-sub_step_sum;
         // Ckeck convergence
         if(R_n<s_tol)                            //accepted
            sub_step_sum += sub_step;
         else
         {
            for(i=0; i<ns; i++)
               TryStress[i] -= 0.5*(dStress0[i]+dStress1[i]);
            continue;
         }
         //

         for(i=0; i<ns; i++)
            DevStress[i] = TryStress[i];

         p_tr = -DeviatoricStress(DevStress)/3.0;

         q_tr = sqrt(1.5*TensorMutiplication2(DevStress, DevStress, dim));

         // If yield, integrate the stress
         F = q_tr*q_tr/M2 + p_tr*(p_tr-p_cn);

         p = p_tr;
         q = q_tr;
         p_c = p_cn;
         vep = 0.0;
         vep0 = 0.0;

         dampFac=1.0;
         NPStep = 0;

         F0 = F;
         if(pcs_deformation==1) F=-1.0;
         if((*data_Plasticity)(3)<MKleinsteZahl)  // p_c0=0
            F = -1.0;
         //TEST CAM-CLAY
         if(p_tr<0)
            F = -1.0;

         if(F>f_tol&&!PreLoad)                    // in yield status
         {
            // Local Newton-Raphson procedure to compute the volume plastic strain
            vep = 0.0;

            //Associative flow rule
            while(isLoop)                         // Newton step for the plastic multiplier
            {
               NPStep++;

               alpha1 = (2.0*p-p_c)/(1.0+(2.0*K+vartheta*p_c)*vep);
               alpha2 = -q/(vep+M2/(G*6.0));
               alpha3 = vartheta*p_c*alpha1;
               alpha1 *= -K;

               dfdp = 2.0*p-p_c;
               dfdq = 2.0*q/M2;

               Jac = alpha1*dfdp+dfdq*alpha2-p*alpha3;

               while(isLoop)                      // Damp
               {

                  vep = vep0 - dampFac*F/Jac;
                  p_c = 0.0;

                  dampFac = 1.0;
                  while(isLoop)                   // Newton step for p_c
                  {
                     NPStep++;
                     if(NPStep>MaxI)
                     {
                        //    printf("\n Too much iteration in Newton step the integration of Cam-Clay \n");
                        //TEST	abort();
                        //test
                        break;
                     }

                     //
                     alpha1 = vartheta*vep/(1.0+2.0*vep*K);
                     alpha2 = p_cn*exp(alpha1*(2.0*p_tr-p_c));
                     beta1 = -alpha1*alpha2-1.0;  //dG(p_c)
                     alpha2 -= p_c;               //G(p_c)
                     p_c -= alpha2/beta1;
                     if(fabs(alpha2)<s_tol) break;
                     if(p_c<0.0) break;
                  }

                  p = (p_tr+vep*K*p_c)/(1.0+2.0*vep*K);
                  q = q_tr/(1.0+6.0*G*vep/M2);

                  F = q*q/M2 + p*(p-p_c);
                  if(F>0.0) break;
                  if(fabs(F/F0)<s_tol) break;
                  dampFac = 0.8;
               }
               vep0 = vep;
               if(fabs(F/F0)<s_tol) break;
            }
         }
         //End substep
         // Plastic strain
         //    alpha1 = 6.0*q*q/(M2*M2*(2.0*p-p_c)*(2.0*p-p_c));
         //    ep += fabs(vep)*sqrt(2.0*(1.0/9.0+alpha1)/3.0);
         ep += 3.0*fabs(vep)*q/M2;
         if(fabs(sub_step_sum-1.0)<DBL_EPSILON)
            break;
      }
      //-------------------------------------------------------------
      // Consistent tangential matrix

      if(Update<1&&NPStep>0)
      {
         alpha1 = 1.0+2.0*K*vep+p_c*vartheta*vep; //a
                                                  //a1
         double a1 = (1.0+p_c*vartheta*vep)/alpha1;
         double a2 = -(2.0*p-p_c)/alpha1;         //a2
         double a3 = 2.0*p_c*vartheta*vep/alpha1; //a3
                                                  //a4
         double a4 = vartheta*p_c*(2.0*p-p_c)/(K*alpha1);
         double a5 = sqrt(1.5)/(1.0+6.0*G*vep/M2);//a5
                                                  //a6
         double a6 =  -3.0*q/(1.0+6.0*G*vep/M2)/M2;

         alpha1 = -4.0*G*q*a6/M2-K*((2.0*a2-a4)*p-a2*p_c);
         double b1 = -K*((a3-2.0*a1)*p+a1*p_c)/alpha1;
         double b2 = 4.0*G*a5*q/(alpha1*M2);

         gamma1 = 2.0*G*q/q_tr;
         gamma2 = K*(a1+a2*b1)-gamma1/3.0;
         gamma3 = -K*a2*b2;
         gamma4 = -2.0*fac*G*a6*b1;
         gamma5 = 2.0*G*fac*(a5+a6*b2)-gamma1;
         //
         // Normalize the stress
         for(i=0; i<ns; i++)
         {
            TryStress[i] = DevStress[i];
            TryStress[i] /=q_tr*fac ;
         }

         // Row 1
         beta1 = gamma2+gamma4*TryStress[0];
         beta2 = gamma3+gamma5*TryStress[0];
         (*Dep)(0,0) = gamma1 + beta1 +beta2*TryStress[0];
         (*Dep)(0,1) =          beta1 +beta2*TryStress[1];
         (*Dep)(0,2) =          beta1 +beta2*TryStress[2];
         (*Dep)(0,3) =                 beta2*TryStress[3];
         if(dim==3)
         {
            (*Dep)(0,4) =              beta2*TryStress[4];
            (*Dep)(0,5) =              beta2*TryStress[5];
         }

         // Row 2
         beta1 = gamma2+gamma4*TryStress[1];
         beta2 = gamma3+gamma5*TryStress[1];
         (*Dep)(1,0) =         beta1 +beta2*TryStress[0];
         (*Dep)(1,1) = gamma1 +beta1 +beta2*TryStress[1];
         (*Dep)(1,2) =         beta1 +beta2*TryStress[2];
         (*Dep)(1,3) =                beta2*TryStress[3];
         if(dim==3)
         {
            (*Dep)(1,4) =             beta2*TryStress[4];
            (*Dep)(1,5) =             beta2*TryStress[5];
         }

         // Row 3
         beta1 = gamma2+gamma4*TryStress[2];
         beta2 = gamma3+gamma5*TryStress[2];
         (*Dep)(2,0) =          beta1 +beta2*TryStress[0];
         (*Dep)(2,1) =          beta1 +beta2*TryStress[1];
         (*Dep)(2,2) =   gamma1+beta1 +beta2*TryStress[2];
         (*Dep)(2,3) =                 beta2*TryStress[3];
         if(dim==3)
         {
            (*Dep)(2,4) =              beta2*TryStress[4];
            (*Dep)(2,5) =              beta2*TryStress[5];
         }

         // Row 4
         beta1 = gamma4*TryStress[3];
         beta2 = gamma5*TryStress[3];
         (*Dep)(3,0) = beta1 + beta2*TryStress[0];
         (*Dep)(3,1) = beta1 + beta2*TryStress[1];
         (*Dep)(3,2) = beta1 + beta2*TryStress[2];
         (*Dep)(3,3) =  gamma1+beta2*TryStress[3];
         if(dim==3)
         {
            (*Dep)(3,4) =      beta2*TryStress[4];
            (*Dep)(3,5) =      beta2*TryStress[5];
            // End row 4

            // Row 5
            beta1 = gamma4*TryStress[4];
            beta2 = gamma5*TryStress[4];
            (*Dep)(4,0) = beta1 + beta2*TryStress[0];
            (*Dep)(4,1) = beta1 + beta2*TryStress[1];
            (*Dep)(4,2) = beta1 + beta2*TryStress[2];
            (*Dep)(4,3) =         beta2*TryStress[3];
            (*Dep)(4,4) = gamma1+ beta2*TryStress[4];
            (*Dep)(4,5) =         beta2*TryStress[5];

            // Row 6
            beta1 = gamma4*TryStress[5];
            beta2 = gamma5*TryStress[5];
            (*Dep)(5,0) = beta1 + beta2*TryStress[0];
            (*Dep)(5,1) = beta1 + beta2*TryStress[1];
            (*Dep)(5,2) = beta1 + beta2*TryStress[2];
            (*Dep)(5,3) =         beta2*TryStress[3];
            (*Dep)(5,4) =         beta2*TryStress[4];
            (*Dep)(5,5) = gamma1+ beta2*TryStress[5];

         }
         //-------------------------------------------------------------

         // Update stresses
         for(i=0; i<ns; i++)
         {
            DevStress[i] /=1.0+6.0*G*vep/M2;
            TryStress[i] = DevStress[i];
         }

         // True stress
         for(i=0; i<3; i++)
            TryStress[i] -= p;
      }
      else
      {
         if(Update<1)
            ElasticConsitutive(dim, Dep);
      }

      for(i=0; i<ns; i++)
         dStrain[i] = TryStress[i];

      // Save the current stresses
      if(Update>0)
      {
         if((*data_Plasticity)(3)<MKleinsteZahl)  // p_c0=0
         {
            p_c = p+q*q/(M2*p);
            (*ele_val->prep0)(GPiGPj) = p_c;
         }
         if(ep>0.0)
         {
            var1 = dev*(1.+e_0);
            e_0 -= var1;
            (*ele_val->pStrain)(GPiGPj) = ep;
            (*ele_val->e_i)(GPiGPj) = e_0;
            (*ele_val->prep0)(GPiGPj) = p_c;
         }
      }

      // For the case of the initial stress being given, the contribution of
      // the initial stress to the right hand side should be taken into account
      // If the initial stress is not accounted into the final stress
      //
      //  for(i=0; i<3; i++)
      //    dStrain[i] -= (*data_Plasticity)(6+i); // Initial stress
      //
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain deduced by creep
   Programing:
   12/2005 WW
   last modified:
   **************************************************************************/
   void CSolidProperties::AddStain_by_Creep(const int ns, double *stress_n,
      double *dstrain, double temperature)
   {
      int i, dim;
      double norn_S, fac=0.0;
      DeviatoricStress(stress_n);
      dim = 2;
      if(ns>4) dim =3;
      norn_S = sqrt(3.0*TensorMutiplication2(stress_n, stress_n, dim)/2.0);
      if(norn_S< DBL_MIN) return;
      switch(Creep_mode)
      {
         case 1:
            //  fac = pow(1.5, (*data_Creep)(1)+1.0)*(*data_Creep)(0)*pow(norn_S, (*data_Creep)(1)-1.0)*dt;
            // fac = pow(2.0/3.0, (*data_Creep)(1))*(*data_Creep)(0)*pow(norn_S, (*data_Creep)(1))*dt;
            fac = (*data_Creep)(0)*pow(norn_S, (*data_Creep)(1))*dt;
            break;
         case 2:
            // gas constant = R = 8.314472(15) J ?K-1 ?mol-1
            // ec= A*exp(-G/RT)s^n
            fac = 1.5*dt*(*data_Creep)(0)*exp(-(*data_Creep)(2)/(8.314472*(temperature+273.15)))*
               pow(norn_S, (*data_Creep)(1));
            break;

      }
      for(i=0; i<ns; i++)
         dstrain[i] -= fac*stress_n[i]/norn_S;
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain induced by HL creep model
   Programing:
   10/2008 UJG/WW
   last modified:
   **************************************************************************/
   void  CSolidProperties::CleanTrBuffer_HL_ODS()
   {
      for(int i=0; i<6; i++)
         (*data_Creep)(i,1) = 0.0;
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain induced by HL creep model
   Programing:
   10/2008 UJG/WW
   last modified:
   **************************************************************************/
   void  CSolidProperties::AccumulateEtr_HL_ODS(const ElementValue_DM *ele_val, const int nGS)
   {
      int i, ns;
      ns = ele_val->xi->Size();
      //
      for(i=0; i<ns; i++)
         (*ele_val->xi)(i) += (*data_Creep)(i, 1)/(double)nGS;

   }
   /**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain induced by HL creep model
   Programing:
   01/2009 UJG/WW
   last modified:
   **************************************************************************/
   void CSolidProperties::CalcYoungs_SVV(const double strain_v)
   {
      double nv = Poisson_Ratio();
      E = (*data_Youngs)(0)/(1.+(*data_Youngs)(1)*pow(strain_v, (*data_Youngs)(2)));
      Lambda = E * nv / ((1. + nv) * (1. - 2. * nv));
      G = 0.5 * E / (1. + nv);
      K=(3.0*Lambda+2.0*G)/3.0;
   }
   /**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain induced by HL creep model
   Programing:
   10/2008 UJG/WW
   last modified:
   **************************************************************************/
   void CSolidProperties::AddStain_by_HL_ODS(const ElementValue_DM *ele_val, double *stress_n,
      double *dstrain, double temperature)
   {
      int i, ns, dim;
      double norn_S, norm_str;
      static double epsilon_tr[6];
      ns = ele_val->xi->Size();
      dim = 2;
      if(ns>4) dim =3;
      DeviatoricStress(stress_n);
      norn_S = sqrt(1.5*TensorMutiplication2(stress_n, stress_n, dim));
      //
      for(i=0; i<ns; i++)
         epsilon_tr[i] = (*ele_val->xi)(i);
      //
      norm_str = sqrt(2.0*TensorMutiplication2(epsilon_tr, epsilon_tr, dim)/3.0);
      double max_etr = norn_S/(*data_Creep)(6, 0)*exp((*data_Creep)(4, 0)*norn_S);
      double eta_k = (*data_Creep)(3, 0)*exp((*data_Creep)(5, 0)*norn_S);
      double eta_m = (*data_Creep)(0, 0)*exp((*data_Creep)(1, 0)*norn_S)
         *exp((temperature+273.16)*(*data_Creep)(2, 0));
      if(max_etr<DBL_EPSILON)
         return;
      for(i=0; i<ns; i++)
      {
         (*data_Creep)(i, 1) += 1.5*dt*(1-norm_str/max_etr)*stress_n[i]/eta_k;
         dstrain[i] -= 1.5*dt*((1-norm_str/max_etr)/eta_k+1/eta_m)*stress_n[i];
      }
   }

   /**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2005 CMCD From fluid properties
   last modified:
   **************************************************************************/
   void CSolidProperties::CalPrimaryVariable(vector<string>& pcs_name_vector)
   {

      CRFProcess* m_pcs = NULL;
      int nidx0,nidx1;
      if(!Fem_Ele_Std)                            //OK
         return;

      for(int i=0;i<(int)pcs_name_vector.size();i++)
      {
         m_pcs = PCSGet(pcs_name_vector[i],true);
         if (!m_pcs) return;                      //MX
         nidx0 = m_pcs->GetNodeValueIndex(pcs_name_vector[i]);
         nidx1 = nidx0+1;

         if(mode==0)                              // Gauss point values
         {
            primary_variable_t0[i]= Fem_Ele_Std->interpolate(nidx0,m_pcs);
            primary_variable_t1[i]= Fem_Ele_Std->interpolate(nidx1,m_pcs);
            primary_variable[i] = (1.-Fem_Ele_Std->pcs->m_num->ls_theta)*Fem_Ele_Std->interpolate(nidx0,m_pcs)
               + Fem_Ele_Std->pcs->m_num->ls_theta*Fem_Ele_Std->interpolate(nidx1,m_pcs);

         }
         else if(mode==2)                         // Element average value
         {
            primary_variable[i] = (1.-Fem_Ele_Std->pcs->m_num->ls_theta)*Fem_Ele_Std->elemnt_average(nidx0,m_pcs)
               + Fem_Ele_Std->pcs->m_num->ls_theta*Fem_Ele_Std->elemnt_average(nidx1,m_pcs);
            primary_variable_t0[i]= Fem_Ele_Std->elemnt_average(nidx0,m_pcs);
            primary_variable_t1[i]= Fem_Ele_Std->elemnt_average(nidx1,m_pcs);
         }
      }
   }

   /**************************************************************************
   FEMLib-Method:
   01/2006 OK Implementation
   05/2009 OK DENSITY
   09/2009 OK Bugfix, write only existing data
   **************************************************************************/
   void CSolidProperties::Write(std::fstream* msp_file)
   {
      //----------------------------------------------------------------------
      // KEYWORD
      *msp_file << "#SOLID_PROPERTIES" << std::endl;
      //-----------------------------------------------------------------------
      //NAME
      if(name.length()>0)
      {
         *msp_file << " $NAME" << std::endl;
         *msp_file << "  ";
         *msp_file << name << std::endl;
      }
      //-----------------------------------------------------------------------
      // GEO_TYPE
      //-----------------------------------------------------------------------
      // DIMENSION
      //-----------------------------------------------------------------------
      // PROPERTIES
      //.......................................................................
      *msp_file << " $DENSITY" << std::endl;
      *msp_file << "  " << Density_mode;
      if(data_Density)                            //OK410
         *msp_file << " " << (*data_Density)(0) << std::endl;
      else
         *msp_file << " Warning: no density data" << std::endl;
      //.......................................................................
      // Elasticity properties
      *msp_file << " $ELASTICITY" << std::endl;
      if(Poisson_Ratio())                         //OK410
         *msp_file << "  POISSION " << Poisson_Ratio() << std::endl;
      *msp_file << "  YOUNGS_MODULUS" << std::endl;
      *msp_file << "  " << Youngs_mode;
      if(data_Youngs)                             //OK410
         *msp_file << " " << (*data_Youngs)(0) << std::endl;
      //.......................................................................
      // Thermal properties
      *msp_file << " $THERMAL" << std::endl;
      if(ThermalExpansion>=0)
      {
         *msp_file << "  EXPANSION" << std::endl;
         *msp_file << "  " << ThermalExpansion << std::endl;
      }
      if(Capacity_mode>0)
      {
         *msp_file << "  CAPACITY" << std::endl;
         *msp_file << "  " << Capacity_mode;
         *msp_file << " " << (*data_Capacity)(0) << std::endl;
      }
      if(this->Conductivity_mode>0)               //OK410
      {
         *msp_file << "  CONDUCTIVITY" << std::endl;
         *msp_file << "  " << Conductivity_mode;
         *msp_file << " " << (*data_Conductivity)(0) << std::endl;
      }
      //-----------------------------------------------------------------------
   }
   /**************************************************************************
   FEMLib-Method:
   03/2008 WW Implementation
   **************************************************************************/
   void CSolidProperties::TEPSwellingParameter_kis(const double suction)
   {
      double val = 0.;
      double alf_i = (*data_Youngs)(1);
      // k_(s)
      // if(suction<=-0.999*1e-6/alf_i)
      val = 1.0+alf_i*suction;
      // else
      //   val = 0.001;
      if(val<0.0) val = 0.001;
      // Swelling index. Kappa
      (*data_Youngs)(6)= val*(*data_Youngs)(0);
   }

   /**************************************************************************
   FEMLib-Method:
   03/2008 WW Implementation
   **************************************************************************/
   double CSolidProperties::TEPSwellingParameter(const double mean_stress)
   {
      double val = 0.;
      double alf_sp = (*data_Youngs)(3);
      double pref = (*data_Youngs)(5);
      double e0 = (*data_Plasticity)(4);
      double suction = (*data_Youngs)(7);
      // k_(s)
      TEPSwellingParameter_kis(suction);
      //
      if(mean_stress<1.0e-20)
         val = 1.0 + alf_sp*log(1.0e-20/pref);
      else if(mean_stress>=pref*exp(-1.0/alf_sp))
      {
         val = 0.;
         return val;
      }
      else
         val = 1. + alf_sp*log(mean_stress/pref);
      return val*(*data_Youngs)(2)*exp((*data_Youngs)(4)*suction)/(3.+3.*e0);
   }

}                                                 // end namespace


/////////////////////////////////////////////////////////////////////////////

/**************************************************************************
FEMLib-Method:
Task: Master read function
Programing:
08/2004 OK Implementation for fluid properties
08/2004 WW Modification for solid properties
01/2005 OK Boolean type
01/2005 OK Destruct before read
**************************************************************************/
bool MSPRead(std::string file_base_name)
{
   //----------------------------------------------------------------------
   //OK  MSPDelete();
   //----------------------------------------------------------------------
   SolidProp::CSolidProperties *m_msp = NULL;
   char line[MAX_ZEILE];
   std::string sub_line;
   std::string line_string;
   std::ios::pos_type position;
   //========================================================================
   // File handling
   std::string msp_file_name = file_base_name + MSP_FILE_EXTENSION;
   std::ifstream msp_file (msp_file_name.data(),std::ios::in);
   if (!msp_file.good())
      return false;
   msp_file.seekg(0L,std::ios::beg);
   //========================================================================
   // Keyword loop
   std::cout << "MSPRead" << std::endl;
   while (!msp_file.eof())
   {
      msp_file.getline(line,MAX_ZEILE);
      line_string = line;
      if(line_string.find("#STOP")!=string::npos)
         return true;
      //----------------------------------------------------------------------
                                                  // keyword found
      if(line_string.find("#SOLID_PROPERTIES")!=std::string::npos)
      {
         m_msp = new SolidProp::CSolidProperties();
         m_msp->file_base_name = file_base_name;
         position = m_msp->Read(&msp_file);
         msp_vector.push_back(m_msp);
         msp_file.seekg(position,std::ios::beg);
      }                                           // keyword found
   }                                              // eof
   return true;
   //========================================================================
}


/**************************************************************************
ROCKFLOW - Funktion: TensorMutiplication2

Aufgabe:
   Calculate tensor mutiplication: a_{ij}b_{ij}

Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E
   const double *s1: Stress tensor in an array (s1_11, s1_22, s1_12, s1_33)
   const double *s2: Stress tensor in an array (s2_11, s2_22, s2_12, s2_33)

Ergebnis:
- double - The value of third stress invariant

Programmaenderungen:
04/2003   WW  Erste Version

**************************************************************************/
double TensorMutiplication2(const double *s1, const double *s2, const int Dim)
{
   switch(Dim)
   {
      case 2:
         return  s1[0]*s2[0]+s1[1]*s2[1]+s1[2]*s2[2]+2.0*s1[3]*s2[3];
         break;
      case 3:
         return  s1[0]*s2[0]+s1[1]*s2[1]+s1[2]*s2[2]+2.0*s1[3]*s2[3]
            +2.0*s1[4]*s2[4]+2.0*s1[5]*s2[5];
         break;
   }
   return 0.0;                                    // To avoid warnings
}


/**************************************************************************
GeoSys: Norm of stress
   06/2008   WW  Programming

**************************************************************************/
double StressNorm(const double *s, const int Dim)
{
   double val = 0.0;
   double mean_s = (s[0]+s[1]+s[2])/3.0;
   double s1 = s[0]-mean_s;
   double s2 = s[1]-mean_s;
   double s3 = s[2]-mean_s;
   val = s1*s1+s2*s2+s3*s3+2.0*s[3]*s[3];
   if(Dim==3)
      val += +2.0*s[4]*s[4]+2.0*s[5]*s[5];
   return val;
}


/**************************************************************************
 ROCKFLOW - Funktion: TensorMutiplication3

Aufgabe:
   Calculate tensor mutiplication: a_{ij}b_{jk}c_{ki} for plane strain problem

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E
   const double *s1: Stress tensor in an array (s1_11, s1_22, s1_12, s1_33)
   const double *s2: Stress tensor in an array (s2_11, s2_22, s2_12, s2_33)
   const double *s3: Stress tensor in an array (s3_11, s3_22, s3_12, s3_33)

Ergebnis:
- double - The value of third stress invariant

Programmaenderungen:
04/2003   WW  Erste Version

**************************************************************************/
double TensorMutiplication3(const double *s1, const double *s2, const double *s3,
const int Dim)
{
   switch(Dim)
   {
      case 2:
         return  (
            s1[0]*(s2[0]*s3[0]+s2[3]*s3[3])       // s1_11*(s2_11*s3_11+s2_12*s3_21)
            +     s1[3]*(s2[0]*s3[3]+s2[3]*s3[1]) // s1_12*(s2_11*s3_12+s2_12*s3_22)
            +     s1[3]*(s2[3]*s3[0]+s2[1]*s3[3]) // s1_21*(s2_21*s3_11+s2_22*s3_21)
            +     s1[1]*(s2[3]*s3[3]+s2[1]*s3[1]) // s1_22*(s2_21*s3_12+s2_22*s3_22)
            +     s1[2]*s2[2]*s3[2]) /3.0;        // s33*s33*s33
         break;
      case 3:
         return  (
         // s1_11*(s2_11*s3_11+s2_12*s3_21+s2_13*s3_31)
            s1[0]*(s2[0]*s3[0]+s2[3]*s3[3]+s2[4]*s3[4])
         // s1_12*(s2_11*s3_12+s2_12*s3_22+s2_13*s3_32)
            +     s1[3]*(s2[0]*s3[3]+s2[3]*s3[1]+s2[4]*s3[5])
         // s1_13*(s2_11*s3_13+s2_12*s3_23+s2_13*s3_33)
            +     s1[4]*(s2[0]*s3[4]+s2[3]*s3[5]+s2[4]*s3[2])
         // s1_21*(s2_21*s3_11+s2_22*s3_21+s2_23*s3_31)
            +     s1[3]*(s2[3]*s3[0]+s2[1]*s3[3]+s2[5]*s3[4])
         // s1_22*(s2_21*s3_12+s2_22*s3_22+s2_23*s3_32)
            +     s1[1]*(s2[3]*s3[3]+s2[1]*s3[1]+s2[5]*s3[5])
         // s1_23*(s2_21*s3_13+s2_22*s3_23+s2_23*s3_33)
            +     s1[5]*(s2[3]*s3[4]+s2[1]*s3[5]+s2[5]*s3[2])
         // s1_31*(s2_31*s3_11+s2_32*s3_21+s2_33*s3_31)
            +     s1[4]*(s2[4]*s3[0]+s2[5]*s3[3]+s2[2]*s3[3])
         // s1_32*(s2_31*s3_12+s2_32*s3_22+s2_33*s3_32)
            +     s1[5]*(s2[4]*s3[3]+s2[5]*s3[1]+s2[2]*s3[5])
         // s1_33*(s2_31*s3_13+s2_32*s3_23+s2_33*s3_33)
            +     s1[2]*(s2[4]*s3[4]+s2[5]*s3[5]+s2[2]*s3[2]))/3.0;
         break;
   }
   return 0.0;                                    // To avoid warnings
}


/**************************************************************************
  ROCKFLOW - Funktion: DeviatoricStress

   Aufgabe:
   Computing the deviatoric stresses
  Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
     double *Stress                   :   Stresses

   Ergebnis:
   - double* - The first stress invariant

Programmaenderungen:
08/2003   WW  Erste Version (for 2D)

**************************************************************************/
double DeviatoricStress(double *Stress)
{
   int i;
   double I1=Stress[0]+Stress[1]+Stress[2];
   for(i=0; i<3; i++)
      Stress[i] -= I1/3.0;

   return I1;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
last modified:
**************************************************************************/
void MSPDelete()
{
   long i;
   int no_msp =(int)msp_vector.size();
   for(i=0;i<no_msp;i++)
   {
      delete msp_vector[i];
   }
   msp_vector.clear();
}


/**************************************************************************
FEMLib-Method:
01/2006 OK Implementation
**************************************************************************/
void MSPWrite(std::string base_file_name)
{
   CSolidProperties* m_msp = NULL;
   //----------------------------------------------------------------------
   // File handling
   std::fstream msp_file;
   std::string msp_file_name = base_file_name + MSP_FILE_EXTENSION;
   msp_file.open(msp_file_name.data(),ios::trunc|ios::out);
   msp_file.setf(ios::scientific,ios::floatfield);
   msp_file.precision(12);
   if (!msp_file.good()) return;
   //----------------------------------------------------------------------
   msp_file << "GeoSys-MSP: Material Solid Properties -------------" << std::endl;
   //----------------------------------------------------------------------
   for(int i=0;i<(int)msp_vector.size();i++)
   {
      m_msp = msp_vector[i];
      m_msp->Write(&msp_file);
   }
   msp_file << "#STOP";
   msp_file.close();
   //----------------------------------------------------------------------
}


/**************************************************************************
FEMLib-Method:
07/2007 OK Implementation
**************************************************************************/
void MSPStandardKeywords()
{
   msp_key_word_vector.clear();
   string in;
   in = "POISSON_RATIO";
   msp_key_word_vector.push_back(in);
   in = "YOUNGS_MODULUS";
   msp_key_word_vector.push_back(in);
}


/**************************************************************************
FEMLib-Method:
07/2007 OK Implementation
**************************************************************************/
CSolidProperties* MSPGet(std::string mat_name)
{
   CSolidProperties *m_msp = NULL;
   for(int i=0;i<(int)msp_vector.size();i++)
   {
      m_msp = msp_vector[i];
      if(mat_name.compare(m_msp->name)==0)
         return m_msp;
   }
   return NULL;
}
