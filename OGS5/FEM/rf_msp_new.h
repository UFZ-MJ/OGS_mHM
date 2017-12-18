/**************************************************************************
FEMLib - Object: MAT-SP
Task: class implementation
Programing:
08/2004 WW Implementation
last modified:
**************************************************************************/
#ifndef rf_msp_new_INC
#define rf_msp_new_INC

// C++ STL
//#include <fstream>
//#include <string>
//#include <vector>

#define MSP_FILE_EXTENSION ".msp"

namespace FiniteElement
{
   class element; class CFiniteElementVec;
   class CFiniteElementStd; class ElementValue_DM;
}

namespace Math_Group {class Matrix;}
namespace process{class CRFProcessDeformation;}

#if defined(WIN32)
class CMATGroupEditorDataEdit;                    //WW
#endif
namespace SolidProp
{

   using FiniteElement::CFiniteElementVec;
   using FiniteElement::ElementValue_DM;
   using FiniteElement::CFiniteElementStd;
   using FiniteElement::ElementValue;
   using Math_Group::Matrix;
   using process::CRFProcessDeformation;
   using ::CRFProcess;
   /*---------------------------------------------------------------*/
   class CSolidProperties
   {
      private:
         // Material parameters
         double PoissonRatio;
         int Youngs_mode;
         int excavation;                          //12.2009. WW
         bool excavated;                          //12.2009. To be ..... WW
         Matrix *data_Youngs;
         double ThermalExpansion;
         //
         double s_tol;                            //16.06.2008 WW
         double f_tol;                            //16.06.2008 WW
         double biot_const;
         double grav_const;                       //WW
         Matrix *data_Density;
         //
         Matrix *data_Capacity;
         Matrix *data_Conductivity;
         //
         Matrix *data_Plasticity;
         Matrix *data_Creep;
         //
         int Density_mode;
         //
         int Capacity_mode;
         int Conductivity_mode;
         int Plasticity_type;
         double primary_variable[10];             //CMCD
         double primary_variable_t0[10];          //CMCD
         double primary_variable_t1[10];          //CMCD
         // Creep property
         // 1. Stationary Norton model
         int Creep_mode;
         //
         bool axisymmetry;

         int mode;                                //CMCD
         // Swelling pressure
         int SwellingPressureType;
         double Max_SwellingPressure;
         //
         std::string CurveVariable_Conductivity;
         int CurveVariableType_Conductivity;
         // Secondary data
         // Elasticity
         double E;                                // Youngs moduls calculated from data_Youngs
         double Lambda;
         double G;                                // Shear stress modulus
         double K;                                // Bulk modulus

         // Rotation matrices and their transpose: UJG 25.11.2009
         Matrix *Crotm;                           // If this is needed by permaebility calculation, we keep it. Otherwise remove it. (To do, UJG/WW)
         Matrix *D_tran;

         // Plasticity
         double dl2;
         // 2. Single yield surface
         Matrix *d2G_dSdS;
         Matrix *d2G_dSdM;
         Matrix *LocalJacobi;                     // To store local Jacobi matrix
         Matrix *inv_Jac;                         // To store the inverse of the  Jacobi matrix
         Matrix *sumA_Matrix;
         double *rhs_l;                           // To store local unknowns of 15
         double *x_l;                             // To store local unknowns of 15
         int  *Li;
         void AllocateMemoryforSYS();
         void ResizeMatricesSYS(const int Dim);

         // Direct stress integration for Drucker-Prager
         double *devS;
         double *dFds;
         double *dGds;
         double *D_dFds;
         double *D_dGds;
         // Mini linear solver
         void Gauss_Elimination(const int DimE, Matrix& AA, int *L,  double *xx);
         void Gauss_Back(const int DimE, Matrix& AA, double * rhs, int *L, double *xx);
         // Thermal properties
         int thermal_conductivity_tensor_type;
         double thermal_conductivity_tensor[9];
         std::string thermal_conductivity_tensor_type_name;
         // Handles. May be used by GUI
         std::string solid_name;
         //-------------------------------------------------------------
         // Numeric
         double CalulateValue(const Matrix *data, const double x) const;
         double Kronecker(const int ii, const int jj);

         // Friends that can access to this data explicitly
         friend bool MSPRead(std::string file_base_name);
         friend void MSPWrite(std::string);
         //WW
         friend class FiniteElement::CFiniteElementVec;
         friend class FiniteElement::CFiniteElementStd;
         friend class FiniteElement::ElementValue;
         friend class process::CRFProcessDeformation;
         friend class ::CRFProcess;
#if defined(WIN32)                          //15.03.2008 WW
         friend class ::CMATGroupEditorDataEdit;
#endif
         //WW
      public:
         //
         CSolidProperties();
         ~CSolidProperties();

         std::ios::pos_type Read(std::ifstream*);
                                                  //CMCD
         FiniteElement::CFiniteElementStd *Fem_Ele_Std;
         std::string name;
         CRFProcess*m_pcs;                        //NW
         // IO
         std::string file_base_name;
         // Output
         void Write(std::fstream*);
                                                  //CMCD
         void CalPrimaryVariable(std::vector<std::string>& pcs_name_vector);

         //-------------------------------------------------------------
         // Access to data
         //-------------------------------------------------------------
         // 1. Density
         double Density(double refence = 0.0);
         // 2. Thermal
         double Heat_Capacity(double refence = 0.0);
         // Boiling model
         double Heat_Capacity(double temperature, double porosity, double Sat);
         int GetCapacityModel() const {return Capacity_mode;}
         int GetConductModel() const {return Conductivity_mode;}
         bool CheckTemperature_in_PhaseChange(const double T0, const double T1);
         double Enthalpy(double temperature, const double latent_factor );
         double Heat_Conductivity(double refence = 0.0);
         void HeatConductivityTensor(const int dim, double* tensor, int group);
         //   int GetCapacityMode() {return Capacity_mode;};  ??
         // 3. Elasticity
#ifdef RFW_FRACTURE
         double Youngs_Modulus(CElem* elem, double refence = 0.0);
                                                  //RFW, for fracture calc
         double Get_Youngs_Min_Aperture(CElem *elem);
#endif
#ifndef RFW_FRACTURE
         double Youngs_Modulus(double refence = 0.0);
#endif
         void SetYoungsModulus(const double El)  {(*data_Youngs)(0)=El;}
         double Poisson_Ratio() const {return PoissonRatio;}
         void CalcYoungs_SVV(const double strain_v);
         double Thermal_Expansion() const {return ThermalExpansion;}
         // 4. Plasticity
         int Plastictity() const {return Plasticity_type;}
         double GetPlasticParameter(const int index) {return (*data_Plasticity)(index);}
         // 5. Creep
         int CreepModel() const {return Creep_mode;}
         double GetCreepParameter(const int index) {return (*data_Creep)(index);}
         // Initilize density
         void NullDensity();

         //-------------------------------------------------------------
         // Manipulators of data
         //-------------------------------------------------------------
         // 1. Elasticity
#ifndef RFW_FRACTURE
         void Calculate_Lame_Constant();
#endif
#ifdef RFW_FRACTURE
         void  Calculate_Lame_Constant(CElem* elem);
#endif
         // For thermal elastic model
         void ElasticConsitutive(const int Dimension, Matrix *D_e) const;
         // For transverse isotropic linear elasticity: UJG 24.11.2009
         void ElasticConstitutiveTransverseIsotropic(const int Dimension);
         Matrix* getD_tran() const {return D_tran;}
         void CalculateTransformMatrixFromNormalVector(const int Dimension);
         // 2. Plasticity
         // 2.1 Drucker-Prager
         double GetAngleCoefficent_DP(const double Angle);
         double GetYieldCoefficent_DP(const double Angle);
         void CalulateCoefficent_DP();
         bool StressIntegrationDP(const int GPiGPj, const ElementValue_DM *ele_val,
            double *TryStress, double& dPhi, const int Update);
         void ConsistentTangentialDP(Matrix *Dep, const double dPhi, const int Dim);
         bool DirectStressIntegrationDP(const int GPiGPj,
            const ElementValue_DM *ele_val, double *TryStress, const int Update);
         void TangentialDP(Matrix *Dep);
         // 2.2 Single yield surface model
         void dF_dNStress(double *dFdS, const double *DevS, const double *S_Invariants,
            const double *MatN1, const int LengthStrs);
         void dF_dStress(double *dFdS, const double *RotV, const double *S_Invariants,
            const double *MatN1, const int LengthStrs);
         void dF_dMat(double *dFdM, const double *S_Invariants, const double *MatN1);
         void dG_dNStress(double *dGdS, const double *DevS, const double *S_Invariants,
            const double *MatN1, const int LengthStrs);
         void dG__dNStress_dNStress(const double *DevS, const double *S_Invariants,
            const double *MatN1, const int LengthStrs );
         void dG__dStress_dStress(const double *DevS, const double *RotV,
            const double *S_Invariants, const double *MatN1,
            const int LengthStrs);
         void dG_dSTress_dMat(const double *DevS, const double *S_Invariants,
            const double *MatN1, const int LengthStrs);
         void dfun2(const double *DevS, const double *RotV, const double *S_Invariants,
            const double *MatN1, const int LengthStrs);

         int CalStress_and_TangentialMatrix_SYS(const int GPiGPj, const ElementValue_DM *ele_val,
            const Matrix *De,  Matrix *D_ep, double *dStress,
            const int Update);
         // 2.2 Cam-clay model
         void CalStress_and_TangentialMatrix_CC(const int GPiGPj,
            const ElementValue_DM *ele_val, double *dStrain,  Matrix *Dep, const int Update);
         // Substep integration. 16.05.2008 WW
         void CalStress_and_TangentialMatrix_CC_SubStep(const int GPiGPj,
            const ElementValue_DM *ele_val, double *dStrain,  Matrix *Dep, const int Update);
         // Parameter function for thermal elatic model. Last modifed on 15.03.2008 //WW
         double TEPSwellingParameter(const double mean_stress);
         void TEPSwellingParameter_kis(const double suction);
         // Strain inrement by creep
         void AddStain_by_Creep(const int ns, double *stress_n, double *dstrain, double temperature=0.0);
         void AddStain_by_HL_ODS(const ElementValue_DM *ele_val, double *stress_n, double *dstrain, double temperature=30);
         void CleanTrBuffer_HL_ODS();
         void AccumulateEtr_HL_ODS(const ElementValue_DM *ele_val, const int nGS);

         // Plasticity
         // 1. Drucker-Prager
         double Al;
         double Xi;
         double Y0;
         double BetaN;
         double Hard;
         double Hard_Loc;

                                                  //WW
         std::vector<std::string>  capacity_pcs_name_vector;
                                                  //WW
         std::vector<std::string>  conductivity_pcs_name_vector;
   };

}                                                 // end namespace
extern std::vector<SolidProp::CSolidProperties*> msp_vector;
extern bool MSPRead(std::string file_base_name);
extern void MSPWrite(std::string);
extern void MSPDelete();
                                                  //OK
extern std::vector<std::string> msp_key_word_vector;
extern void MSPStandardKeywords();                //OK
                                                  //OK
extern SolidProp::CSolidProperties* MSPGet(std::string);

extern double StressNorm(const double *s, const int Dim);
extern double TensorMutiplication2(const double *s1, const double *s2, const int Dim);
extern double TensorMutiplication3(const double *s1, const double *s2, const double *s3, const int Dim);
extern double DeviatoricStress(double *Stress);
#endif
