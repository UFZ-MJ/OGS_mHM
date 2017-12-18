/*
   Class element declaration
   class for finite element.
   Designed and programmed by WW/OK, 12/2004
   modified by TF, 10/2010
*/

#ifndef fem_std_INC
#define fem_std_INC

//#include "FEMEnums.h"

#include "fem_ele.h"
// Problems
#include "rf_mfp_new.h"
#include "rf_msp_new.h"
#include "rf_out_new.h"                           //OK

//-----------------------------------------------------
// Process type
//L: Liquid flow
//U: Unconfined flow
//G: Groundwater flow
//T: Two-phase flow
//C: Componental flow
//H: heat transport
//M: Mass transport
//O: Overland flow
//R: Richards flow
//F: Fluid momentum
//A: Gas flow
enum EnumProcessType { L, U, G, T, C, H, M, O, R, F, A, V, P};
//-----------------------------------------------------

namespace process {class CRFProcessDeformation;}
using SolidProp::CSolidProperties;
// Predeclared classes  01/07, WW
class CMediumProperties;
class CSolidProperties;
class CFluidProperties;

class CRFProcess;
namespace FiniteElement
{
   using Math_Group::SymMatrix;
   using Math_Group::Matrix;
   using Math_Group::Vec;
   using process::CRFProcessDeformation;
   using ::CRFProcess;

   class CFiniteElementStd : public CElement
   {
      public:
         CFiniteElementStd(CRFProcess *Pcs, const int C_Sys_Flad, const int order=1);
         ~CFiniteElementStd();

         // Set material data
         void SetMaterial(const int phase=0);
         // Set memory for local matrices
         void SetMemory();
         // Set variable  YD
         void SetVariable();

         // Set coupling information
         void ConfigureCoupling(CRFProcess* pcs, const int *Shift, bool dyn=false);

         // Element claculation
         // 1. Mass matrix
         void CalcMass();
         void CalcMass2();
         void CalcMassPSGLOBAL();                 // PCH
         // 2. Lumped mass matrix
         void CalcLumpedMass();
         void CalcLumpedMass2();
         void CalcLumpedMassPSGLOBAL();           // PCH
         // 3. Laplace matrix
         void CalcLaplace();
         // 4. Gravity term
         void CalcGravity();
         // 5. Strain coupling matrix
         void CalcStrainCoupling();
         // 6. Thermal coupling
         void CalcRHS_by_ThermalDiffusion();
         // 7. Advection matrix
         void CalcAdvection();
         // 8. Storage matrix
         void CalcStorage();
         // 9. Content matrix
         void CalcContent();
         //
         void CalcSatution();                     //WW
         //
                                                  //25.08.2008. WW
         void CalcEnergyNorm(const double *x_n1, double &err_norm0, double &err_normn);
                                                  //25.09.2008. WW
         void CalcEnergyNorm_Dual(const double *x_n1, double &err_norm0, double &err_normn);
         //
         void CalcNodeMatParatemer();             //WW
         // Assembly
         void Assembly();
         void Assembly(int option, int dimension);// PCH for Fluid Momentum
         void Cal_Velocity();
         void Cal_Velocity_2();                   //CB this is to provide velocity only at the element center of gravity
         void  Cal_GP_Velocity_FM(int *i_ind);    //SB 4900 interpolate node velocities to Gauss point velocities
         //
                                                  //OK
         void AssembleParabolicEquationRHSVector();

         // CVFEM functions for overland flow   JOD
                                                  // to move
         void GetOverlandBasisFunctionMatrix_Line();
                                                  // to move
         void GetOverlandBasisFunctionMatrix_Quad();
         void CalcOverlandCoefficients(double* head, double* axx, double* ayy, double* ast);
         void CalcOverlandCoefficientsLine(double* head, double* axx, double* ast );
         void CalcOverlandCoefficientsQuad(double* head, double* axx, double* ayy, double* ast );
         void CalcOverlandCoefficientsTri(double* head, double* axx, double* ayy, double* ast );
         void CalcOverlandNLTERMS(double* H, double* HaaOld, double* swval, double* swold);
         void CalcOverlandNLTERMSRills(double* H, double* HaaOld, double* swval, double* swold);
         void CalcOverlandNLTERMSChannel(double* H, double* HaaOld, double* swval, double* swold);
         void CalcOverlandCKWR(double* head, double* ckwr, int* iups);
         void CalcOverlandCKWRatNodes(int i, int j, double* head, double* ckwr, int* iups);
         void CalcOverlandResidual(double* head, double* swval, double* swold, double ast, double* residuall, double** amat);
         double CalcOverlandJacobiNodes(int i, int j, double *depth, double *depth_keep, double akrw, double axx, double ayy, double** amatij, double* sumjac);
         void CalcOverlandUpwindedCoefficients(double** amat, double* ckwr, double axx, double ayy);
         //
                                                  //CB added by CB: 090507
         inline void UpwindAlphaMass(double *alpha);
                                                  //CB added by CB: 090507
         inline void UpwindUnitCoord(int p, int point, int ind);
         int UpwindElement(int option, int phase);// PCH
                                                  //CB added by CB: 090507
         inline void UpwindSummandMass(const int gp, int& gp_r, int& gp_s, int& gp_t, double *alpha, double *summand);
                                                  //NW
         inline double CalcSUPGCoefficient(double*vel,int ip);
                                                  //NW
         inline void CalcSUPGWeightingFunction(double*vel, int ip, double &tau, double*v_dN);
                                                  //NW
         inline double CalcSUPGEffectiveElemenetLength(double *vel);
         // Gauss value
         void ExtropolateGauss(CRFProcess *m_pcs, const int idof);
         //

      private:
         bool newton_raphson;                     //24.05.2007 WW
         long index;
         int dof_index;                           //24.02.2007 WW
         // Column index in the node value table
         int idx0, idx1, idxS, idxSn0, idxSn1, idx3;
         int idxp0,idxp1, idxp20, idxp21;
         int phase;
         int comp;                                // Component
         int LocalShift;                          // For RHS
         // Danymic
         int idx_vel_disp[3], idx_pres;
         // Velocity
         int idx_vel[3];                          //WW
         // Material properties
         double mat[9];
         double *eqs_rhs;                         //For DDC WW
         bool heat_phase_change;
         ::CRFProcess *cpl_pcs;                   // Pointer to coupled process. WW

         char pcsT;
         //     /**
         //      * process type, \sa enum ProcessType
         //      */
         //     ProcessType _pcs_type; // TF

         bool dynamic;
         CRFProcess *mfp_pcs;
         CSolidProperties *SolidProp;
         CFluidProperties *FluidProp;
         CFluidProperties *GasProp;
         CMediumProperties *MediaProp;
         CMediumProperties *MediaProp1;           // Matrix for the dual model. YD/WW
         CSolidProperties *SolidProp1;            // Matrix for the dual model. YD/WW
         CRFProcess *pcs;
         process::CRFProcessDeformation *dm_pcs;
         bool flag_cpl_pcs;                       //OK
         //-------------------------------------------------------
         // Auxillarary matrices
         Matrix *StiffMatrix;
         Matrix *AuxMatrix;
         Matrix *AuxMatrix1;
         // Gravity matrix;
         //25.2.2007.WW  SymMatrix *GravityMatrix;
         // Gauss point value. Buffers. // Some changes. 27.2.2007 WW
         double TG, TG0, PG, PG0, PG2,PG20, drho_gw_dT;
         double Sw, rhow, poro, dSdp;
         double rho_gw, rho_ga, rho_g, p_gw, M_g, tort;
         //
         double edlluse[16];
         double edttuse[16];

         // Local matrices
         Matrix *Mass;                            //CB symMatrix *Mass; // unsymmetric in case of upwinding
         Matrix *Mass2;
         Matrix *Laplace;
         Matrix *Advection;                       //SB4200
         Matrix *Storage;                         //SB4200
         Matrix *Content;                         //SB4209
         Matrix *StrainCoupling;
         Vec       *RHS;
         //-------------------------------------------------------
         void SetHighOrderNodes();                // 25.2.2007 WW
         // Primary as water head
         bool HEAD_Flag;
         //
         void Config();
         //
         inline double CalCoefMass();
                                                  // 25.2.2007 WW
         inline double CalCoefMass2(int dof_index);
                                                  // 03.3.2009 PCH
         inline double CalCoefMassPSGLOBAL(int dof_index);
         inline void CalCoefLaplace(bool Gravity, int ip=0);
                                                  // 10 2008 PCH
         inline void CalCoefLaplaceMultiphase(int phase, int ip=0);
         inline void CalCoefLaplace2(bool Gravity, int dof_index);
                                                  // 03 2009 PCH
         inline void CalCoefLaplacePSGLOBAL(bool Gravity, int dof_index);
         inline double CalCoefAdvection();        //SB4200 OK/CMCD
         inline double CalCoefStorage();          //SB4200
         inline double CalCoefContent();
         inline double CalCoefStrainCouping();
         inline double  CalcCoefDualTransfer();
                                                  // 27.2.2007 WW
         inline double CalCoef_RHS_T_MPhase(int dof_index);
                                                  // 27.2.2007 WW
         inline double CalCoef_RHS_M_MPhase(int dof_index);
         inline double CalCoef_RHS_PSGLOBAL(int dof_index);
                                                  //  NB
         inline double CalCoef_RHS_T_PSGlobal(int dof_index);
                                                  // 03.2007 PCH
         inline void CalCoef_RHS_Pc(int dof_index);
                                                  //AKS
         inline double CalCoef_RHS_AIR_FLOW(int dof_index);
                                                  //AKS
         inline double CalCoef_RHS_HEAT_TRANSPORT(int dof_index);
                                                  //AKS
         inline double CalCoef_RHS_HEAT_TRANSPORT2(int dof_index);
         inline void CalNodalEnthalpy();
         //-----------------------------------------------------
         // Process type
         //L: Liquid flow
         //U: Unconfined flow
         //G: Groundwater flow
         //T: Two-phase flow
         //C: Componental flow
         //H: heat transport
         //M: Mass transport
         //O: Liquid flow
         //R: Richards flow
         //A: Gas flow
         //F: Fluid Momentum
         EnumProcessType PcsType;
         //-----------------------------------------------------
         // Local Assembly
         // Assembly of parabolic equation
         void AssembleParabolicEquation();        //OK4104
                                                  //SB4200
         void AssembleMixedHyperbolicParabolicEquation();
         void AssembleParabolicEquationNewton();
                                                  // JOD
         void AssembleParabolicEquationNewtonJacobian(double** jacob, double* Haa, double* HaaOld, double axx, double ayy, double** amat, double ast, double* swold, double* residuall, int* iups);
         inline void Assemble_strainCPL();        // Assembly of strain coupling
         void AssembleMassMatrix(int option);     // PCH
         // Assembly of RHS by Darcy's gravity term
         void Assemble_Gravity();
         void Assemble_Gravity_Multiphase();
         // Assembly of RHS by temperature for m-phase flow 27.2.2007 WW
         void Assemble_RHS_T_MPhaseFlow();
         // Assembly of RHS by deformation. 27.2.2007 WW
         void Assemble_RHS_M();
         void Assemble_RHS_Pc();                  // 03.2009 PCH
         void Assemble_RHS_AIR_FLOW();            //AKS
         void Assemble_RHS_HEAT_TRANSPORT();      //AKS
         void Assemble_RHS_HEAT_TRANSPORT2();     //AKS
         void Assemble_RHS_T_PSGlobal();          // Assembly of RHS by temperature for PSGlobal
         void AssembleRHS(int dimension);         // PCH
         void Assemble_DualTransfer();
         bool check_matrices;                     //OK4104
         void AssembleRHSVector();                //OK
         void AssembleCapillaryEffect();          // PCH
                                                  // PCH for debugging
         void PrintTheSetOfElementMatrices(std::string mark);
         // Friend classes, 01/07, WW
         friend class ::CMediumProperties;
         friend class SolidProp::CSolidProperties;
         friend class ::CFluidProperties;
         // Friend functions. WW
         friend double ::MFPCalcFluidsHeatCapacity(CFiniteElementStd* assem);

         // Auxillarary vectors for node values
         // Vector of local node values, e.g. pressure, temperature.
         // Assume maximium element nodes is 20
         //double OldMatrix[64]; // For grid adapting
         double *NodalVal;
         double *NodalVal0;                       //?? NodalValueSaturation, NodalValueTemperature; ...
         double *NodalVal1;
         double *NodalVal2;
         double *NodalVal3;
         double *NodalVal4;
         double *NodalValC;
         double *NodalValC1;
         double *NodalVal_Sat;
         double *NodalVal_SatNW;
         double *NodalVal_p2;
         double *NodalVal_p20;                    //AKS
         //
         double *weight_func;                     //NW
         //
         friend class ::CRFProcess;
   };

   // Vector for storing element values WW
   class ElementValue
   {
      public:
         ElementValue(CRFProcess* m_pcs, CElem* ele);
         ~ElementValue();
         void getIPvalue_vec(const int IP, double * vec);
         void GetEleVelocity(double * vec);
      private:
         // Friend class
         friend class ::CRFProcess;
         friend class FiniteElement::CFiniteElementStd;
         friend class ::COutput;                  //OK
         // Process
         CRFProcess *pcs;
         // Data
         Matrix Velocity;
         Matrix Velocity_g;
   };
}                                                 // end namespace


/*------------------------------------------------------------------
   Finite element calculation for standard PDE.
   12.12.2004 WW
------------------------------------------------------------------*/
#endif
