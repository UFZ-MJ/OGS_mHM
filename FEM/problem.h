/*
  Class to define a problem to be solved.
  Programming WW 08.07.2008
*/
#ifndef problem_INC
#define problem_INC

#include "Configure.h"

#include <vector>
class CRFProcess;

// GEOLIB
#include "GEOObjects.h"

//---------------------------------------------------------------------
//Pointers to member functions
class Problem;
typedef  double (Problem::*ProblemMemFn)(void);
#define Call_Member_FN(object,ptrToMember)  ((object)->*(ptrToMember))
//---------------------------------------------------------------------
class Problem
{
   public:
      Problem(char* filename = NULL);
      ~Problem();
      void Euler_TimeDiscretize();
      void RosenBrock_TimeDiscretize() {};
      //
      void SetActiveProcesses();
      void PCSRestart();
      //
      bool CouplingLoop();
      void PostCouplingLoop();
      // Copy u_n for auto time stepping
      double* GetBufferArray() {return buffer_array;};

      /**
       * get the geometric objects stored in GEOLIB::GEOObjects
       * @return a pointer to an instance of class GEOLIB::GEOObjects
       */
      const GEOLIB::GEOObjects* getGeoObj () const;
      /**
       * Get the name of the project. The name is used by GEOLIB::GEOObjects
       * to access the geometric data.
       * @return the name to acces geometric data
       */
      const std::string& getGeoObjName () const;

#ifdef BRNS
      // BRNS-Coupling: For writing spatially resolved reaction rates at the final iteration,
      // we need to get the timing information.

      double getCurrentTime();
      double getEndTime();
#endif                                         //BRNS

   private:
      // Time:
      double start_time;
      double end_time;
      double current_time;
      double *buffer_array;
      int step_control_type;
      // Mixed time step WW
      double dt0;                                 // Save the original time step size

      // Controls
      int loop_index;
      int max_coupling_iterations;
      size_t max_time_steps;
      double coupling_tolerance;
      //
      int lop_coupling_iterations;
      bool CalcVelocities;
      bool conducted;

      // Print flag
      bool print_result;
      // Processes

      std::vector<CRFProcess*> total_processes;
      std::vector<CRFProcess*> transport_processes;
      std::vector<CRFProcess*> multiphase_processes;
      ProblemMemFn *active_processes;
      std::vector<int> active_process_index;
      std::vector<int> coupled_process_index;
      //
      bool *exe_flag;
      inline int AssignProcessIndex(CRFProcess *m_pcs, bool activefunc = true);
      //
      void PCSCreate();
      // Perform processes:
      inline double LiquidFlow();
      inline double RichardsFlow();
      inline double TwoPhaseFlow();
      inline double MultiPhaseFlow();
      inline double PS_Global();                  // 03 2009 PCH
      inline double GroundWaterFlow();
      inline double ComponentalFlow();
      inline double OverlandFlow();
      inline double AirFlow();
      inline double HeatTransport();
      inline double FluidMomentum();
      inline double RandomWalker();
      inline double MassTrasport();
      inline double Deformation();
      // Accessory
      void LOPExecuteRegionalRichardsFlow(CRFProcess*m_pcs_global);
      void LOPCalcELEResultants();
      inline void ASMCalcNodeWDepth(CRFProcess *m_pcs);
      void PCSCalcSecondaryVariables();
      bool Check();                               //OK

      /**
       * pointer to an instance of class GEOObjects,
       * that manages geometric entities
       */
      GEOLIB::GEOObjects* _geo_obj;               // TF
      /**
       * project/file name for geometry file,
       * used to access data in data manager GEOObjects
       */
      std::string _geo_name;                      // TF
};

extern Problem *aproblem;
extern bool MODCreate();                          //OK
#endif
