/**************************************************************************
FEMLib - Object: MAT-FP
Task: class implementation
Programing:
08/2004 OK Implementation
last modified:
**************************************************************************/
#ifndef rf_mfp_new_INC
#define rf_mfp_new_INC

//#include "rfmat_cp.h"//AKS
class CompProperties;

/*---------------------------------------------------------------*/
//namespace FiniteElement {class CElement;}
//using FiniteElement::CElement;
namespace FiniteElement {class CFiniteElementStd;}
using FiniteElement::CFiniteElementStd;
class CFluidProperties
{
   private:
      // State variables
      double primary_variable[10];                //WW
      double primary_variable_t0[10];             //CMCD
      double primary_variable_t1[10];             //CMCD
      CRFProcess *mfp_pcs;
      // CElement *m_element;
      //
      bool cal_gravity;                           //YD/WW
      // FEM
                                                  //WW
      friend class  FiniteElement::CFiniteElementStd;
      // PCS
      CRFProcess *m_pcs;                          //OK4704
   public:

      int fluid_id;                               // specification of substance (NB JUN 09)
      double rhoc;                                //critical_density; //NB
      double Tc;                                  //critical_temperature;
      double pc;                                  //critical_pressure;
      double Tt;                                  //triple_point_temperature;
      double pt;                                  //triple_point_pressure;
      double Rs;                                  //specific_gas_constant;
      double Ru;                                  //universal_gas_constant;
      double omega;                               // azentric factor for Peng-Robinson EOS
      double molar_mass;
      std::string name;
      std::string fluid_name;                     //NB4801
      // Limits and coefficients for free Helmholtz Energy, NB JUN 09
      int limit[5];
      double k [2][8];
      double K [14][56];
      // compressibility
      int compressibility_model_pressure;         //NB
      int compressibility_model_temperature;      //NB
      int compressibility_pressure;               //NB
      int compressibility_temperature;            //NB

      int phase;

      // FEM
      CFiniteElementStd *Fem_Ele_Std;
      long node;                                  //OK4704
      // Density
      int density_model;
      double rho_0;
      double drho_dp;
      double drho_dT;

      double drho_dC;
      std::string rho_fct_name;
      // Viscosity
      int viscosity_model;
      double viscosity;
      double my_0;
      double dmy_dp;
      double dmy_dT;
      double dmy_dC;
      std::string my_fct_name;
      // Thermal properties

      double specific_heat_capacity;
      std::string heat_capacity_fct_name;
      int heat_conductivity_model;
      double heat_conductivity;
      std::string heat_conductivity_fct_name;
      int heat_diffusion_model;                   //AKS
      double temperature_buffer;                  //YD, shifted to public JOD
      int heat_capacity_model;                    //YD, shifted to public JOD
      // Electrical properties
      int electric_conductivity_model;
      int electric_conductivity_num_val;
      double *electric_conductivity_val;
      // Chemical properties
      std::string dif_fct_name;
      int diffusion_model;                        /* SB:p2 */
      double diffusion;                           /*SB:2p */
      // State variables
      double p_0;
      double T_0;
      double C_0;
      double Z;
      double T_Latent1, T_Latent2, latent_heat;
      int heat_phase_change_curve;
      // IO
      std::string file_base_name;
      int mode;
      // PCS  YD
      std::vector<std::string>density_pcs_name_vector;
      std::vector<std::string>viscosity_pcs_name_vector;
      std::vector<std::string>specific_heat_capacity_pcs_name_vector;
      std::vector<std::string>heat_conductivity_pcs_name_vector;
      std::vector<std::string>enthalpy_pcs_name_vector;
                                                  //AKS
      std::vector<CompProperties*>component_vector;
      // DAT
      //string dat_type_name;
      //string dat_name;
      // FCT
      std::string fct_pcs_name;
      std::string fct_name;
      //
      CFluidProperties(void);
      ~CFluidProperties(void);
      std::ios::pos_type Read(std::ifstream*);
      void Write(std::ofstream*);
      void CalPrimaryVariable(std::vector<std::string>& pcs_name_vector);
      // Add an argument: double* variables = NULL. 28.05.2008 WW
      double Density(double *variables = NULL);
      double drhodP (double P, double T);
                                                  //AKS
      double MixtureSubProperity(int properties, long idx_elem, double p, double T);
      double drhodT (double P, double T);
      double Viscosity(double *variables = NULL); //OK4709
                                                  //NB Jan09
      double SpecificHeatCapacity(double *variables = NULL);
      void therm_prop(std::string caption);       //NB 4.9.05
      double PhaseChange();                       // JOD
      double HeatConductivity(double *variables = NULL);
      double CalcEnthalpy(double temperature);
      //WW double Enthalpy(int,double);
      //WW double EnthalpyPhase(long,int,double*,double);
      //WW double MassFraction(long number,int comp,double*gp,double theta, CFiniteElementStd* assem=NULL);
      //WW double InternalEnergy(long number,double*gp,double theta);
      //WW double DensityTemperatureDependence(long number,int comp,double*gp,double theta);
      double vaporDensity(const double T);        //WW
                                                  //WW
      double vaporDensity_derivative(const double T);
      bool CheckGravityCalculation() const {return cal_gravity;}
      int GetHeatCapacityModel() const            //YD
      {
         return heat_capacity_model;
      }
      // Derivations of free Helmholtz energy, NB JUN 09
      double phi_r_d (double rho, double T, int c);
      double phi_r_tt (double rho, double T, int c);
      double phi_0_t (double T,int c);
      double phi_r_t (double rho, double T,int c);
      double phi_r_dt (double rho, double T, int c);
      double phi_r_dd (double rho, double T, int c);
      double phi_0_tt (double T, int c);
                                                  //AKS
      double CalCopressibility(long idx_elem, double p, double T);

   private:
      double GasViscosity_Reichenberg_1971(double,double);
                                                  //AKS
      double GasViscosity_Chung_1988(long idx_elem, double,double);
      double MATCalcFluidDensityMethod8(double p, double T, double C);
      double LiquidViscosity_Yaws_1976(double);
      double LiquidViscosity_Marsily_1986(double);
      double LiquidViscosity_NN(double,double);
      double LiquidViscosity_CMCD(double p,double T,double C);
      double MATCalcHeatConductivityMethod2(double p, double T, double C);
      double MATCalcFluidHeatCapacityMethod2(double p, double T, double C);
};

extern std::vector<CFluidProperties*>mfp_vector;
extern bool MFPRead(std::string);
extern void MFPWrite(std::string);
#define MFP_FILE_EXTENSION ".mfp"
//WW extern double MFPCalcVapourPressure(double);
                                                  //WW
extern double MFPCalcFluidsHeatCapacity(CFiniteElementStd* assem=NULL);
extern double MFPCalcFluidsHeatConductivity(long index,double*gp,double theta, CFiniteElementStd* assem=NULL);
extern void MFPDelete();
                                                  //OK/YD
extern CFluidProperties* MFPGet(const std::string&);
extern CFluidProperties* MFPGet(int);             //NB JUN 09
                                                  //NB AUG 09
double MFPGetNodeValue(long,const std::string&,int);
#endif
