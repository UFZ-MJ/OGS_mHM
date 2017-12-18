/**************************************************************************
ROCKFLOW - Object: Process PCS_Deformation
Task:
Programing:
07/2003 WW Implementation
**************************************************************************/
#ifndef pcs_dm_INC
#define pcs_dm_INC

#include "Configure.h"

#include <vector>

#include "rf_pcs.h"

// Strong discontinuity
extern bool Localizing;                           // for tracing localization
typedef struct
{
   int ElementIndex;
   int NumInterFace;                              // Number of intersection faces
   // Local indeces of intersection faces (3D)
   int *InterFace;
} DisElement;
extern std::vector<DisElement*> LastElement;      // Last discontinuity element correponding to SeedElement
extern std::vector<long> ElementOnPath;           // Element on the discontinuity path

namespace FiniteElement {class CFiniteElementVec; }
using FiniteElement::CFiniteElementVec;
class CPARDomain;

namespace process
{

   // Elasto-plastic Deformation
   class CRFProcessDeformation : public CRFProcess
   {
      public:
         CRFProcessDeformation();
         virtual ~CRFProcessDeformation();

         void Initialization();

         // Assemble system equation
         void GlobalAssembly();

                                                  // overloaded
         double Execute(const int CouplingIterations=0);

         // Aux. Memory
         double* GetAuxArray() const {return ARRAY;};

         void ScalingNodeForce(const double SFactor);
         void InitGauss();
         //
         void SetInitialGuess_EQS_VEC();
         void UpdateIterativeStep(const double damp, const int Type);
         void InitializeNewtonSteps(const int type);
         double NormOfUpdatedNewton();
         void StoreLastSolution(const int ty=0);
         void RecoverSolution(const int ty=0);
         double NormOfDisp();
#ifndef NEW_EQS
         double NormOfUnkonwn();
#endif
         // Stress
         // For partitioned HM coupled scheme
         void ResetCouplingStep();
         void ResetTimeStep();
         //
         void UpdateStress();
         void UpdateInitialStress(bool ZeroInitialS);
         void Extropolation_GaussValue();

         // Excavation computation
         void ReleaseLoadingByExcavation();
         void CreateInitialState4Excavation();

         // Dynamic
         bool CalcBC_or_SecondaryVariable_Dynamics(bool BC = false);
         // Calculate scaling factor for load increment
         double CaclMaxiumLoadRatio() ;

         // Write stresses
         void WriteGaussPointStress();
         void ReadGaussPointStress();
      public:                                     //OK
         CFiniteElementVec *fem_dm;
      private:
         void InitialMBuffer();
         double *ARRAY;
         int *unknown_vector_indeces;
         long *unknown_node_numbers;
         int counter;
         double InitialNorm;
         double InitialNormU;
         double InitialNormU0;

         //
         double error_k0;
         // Domain decompisition
         inline void DomainAssembly(CPARDomain* m_dom);

         // For strong discontinuity approach
         void Trace_Discontinuity();
         long MarkBifurcatedNeighbor(const int PathIndex);
   };
}                                                 // end namespace


extern void CalStressInvariants(const long Node_Inex,double *StressInv);
                                                  // For visualization
extern void CalMaxiumStressInvariants(double *StressInv);
extern double LoadFactor;
extern double Tolerance_global_Newton;
extern double Tolerance_Local_Newton;
extern int enhanced_strain_dm;
extern int number_of_load_steps;
extern int problem_dimension_dm;
extern int PreLoad;
extern bool GravityForce;
#endif
