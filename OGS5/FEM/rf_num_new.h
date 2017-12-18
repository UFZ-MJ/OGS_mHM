/**************************************************************************
FEMLib - Object: NUM
Task: class implementation
Programing:
11/2004 OK Implementation
last modified:
**************************************************************************/
#ifndef rf_num_new_INC
#define rf_num_new_INC

#include "Configure.h"

#define NUM_FILE_EXTENSION ".num"
// C++ STL
#include <fstream>
#include <string>
#include <vector>
#include "prototyp.h"
//----------------------------------------------------------------
class CNumerics
{
   private:
      // cf. Computational Geomachanics pp.62 WW
      double *DynamicDamping;
   public:
      // method
      std::string method_name;                    //OK
      // PCS
      std::string pcs_type_name;
      // RENUMBER
      int renumber_method;
      int renumber_parameter;
      // LS - Linear Solver
      int ls_method;
      int ls_max_iterations;
      int ls_error_method;
      double ls_error_tolerance;
      double ls_theta;
      int ls_precond;
      int ls_storage_method;
      // LS - Linear Solver
      std::string nls_method_name;
      int nls_method;                             // Picard or Newton
      int nls_error_method;                       //WW
      int nls_max_iterations;
      double nls_error_tolerance;
      double nls_error_tolerance_local;
      double nls_relaxation;
      // CPL WW
      double cpl_tolerance;
      int  cpl_iterations;
      std::string  cpl_variable;                  // MB
      // ELE
      int ele_gauss_points;                       // probably element-type-wise
      int ele_mass_lumping;
      int ele_upwind_method;                      //CB
      double ele_upwinding;
      int ele_supg_method;                        //NW
      int ele_supg_method_length;                 //NW
      int ele_supg_method_diffusivity;            //NW
      // Deformation
      int GravityProfile;
      // LAGRANGE method //OK
      double lag_quality;
      int lag_max_steps;
      double lag_local_eps;
      int lag_time_weighting;
      double lag_min_weight;
      int lag_use_matrix;
      int lag_vel_method;
      //
      // Dynamics
      bool CheckDynamic();
      double GetDynamicDamping_beta1 () const {return DynamicDamping[0];}
      double GetDynamicDamping_beta2 () const {return DynamicDamping[1];}
      double GetDynamicDamping_bbeta () const {return DynamicDamping[2];}
      //
      CNumerics(std::string);
      ~CNumerics(void);
      std::ios::pos_type Read(std::ifstream*);
      void Write(std::fstream*);
};

extern std::vector<CNumerics*>num_vector;
extern bool NUMRead(std::string);
extern void NUMWrite(std::string);
extern void NUMDelete();
extern CNumerics* NUMGet(std::string);

//////////////////////////////////////////////////////////////////////////
// SOLVER
//////////////////////////////////////////////////////////////////////////
typedef struct
{
   char *name;
   long type;
   long maxiter;
   double eps;
   long nom;
   long precond;
   long store;
   long criterium;
   double theta;
   long repeat;
   double time;
   long kind;
   long level;

} LINEAR_SOLVER_PROPERTIES;

typedef struct
{
   char pcs_type_name[80];
   void *matrix;
   double *b;
   double *x;
   long dim;
   long store;
   LINEAR_SOLVER_PROPERTIES *lsp;
   char *lsp_name;
   double *xx;
   double *r;
   double *memory;
   int memory_number;
   double **new_memory;
   void (*init_function) ();
   void (*assemble_function) (double *, double *, double);
   long assemble_index;
   long level;
   IntFuncDXDXL LinearSolver;
   long master_iter;
   long iter_sum;
   char *system_time_name;
   char *system_time_assemble_function_name;
   char *name_group_ls;
   char *name_ls;
   long number_ls;
   long num_of_unknowns_ls;
   //OK UNKNOWN_LINEAR_SOLVER **unknown_ls;
   int unknown_vector_dimension;                  /* nodal degree of freedom */
   int *unknown_vector_indeces;                   /* pointer of field unknown_vector_index[unknown_vector_dimension]   */
   long *unknown_node_numbers;                    /* pointer of field unknown_node_numbers[unknown_vector_dimension]   */
   int *unknown_update_methods;                   /* pointer of field unknown_update_methods[unknown_vector_dimension] */
} LINEAR_SOLVER;

#ifdef USE_MPI                                    //WW
extern LINEAR_SOLVER *InitVectorLinearSolver(LINEAR_SOLVER*);
#endif
#ifndef NEW_EQS                                   //WW 07.11.2008
extern LINEAR_SOLVER *InitLinearSolver(LINEAR_SOLVER*);
//
extern void SetLinearSolverType(LINEAR_SOLVER*,CNumerics*);
extern LINEAR_SOLVER *InitializeLinearSolver(LINEAR_SOLVER*,CNumerics*);
extern LINEAR_SOLVER *InitMemoryLinearSolver(LINEAR_SOLVER*,int);
extern LINEAR_SOLVER *SetMemoryZeroLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER *SetZeroLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER *DestroyLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER *DestroyMemoryLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER *SetLinearSolver(LINEAR_SOLVER*);
extern LINEAR_SOLVER *CreateLinearSolver(long store, long dim);
extern LINEAR_SOLVER *CreateLinearSolverDim(long store, int unknown_vector_dimension, long dim);
extern void ConfigSolverProperties(void);
//
extern double CalcNormOfRHS(LINEAR_SOLVER*);
//
extern int GetUnknownVectorDimensionLinearSolver(LINEAR_SOLVER*);
#endif                                            //ifndef NEW_EQS //WW 07.11.2008

//////////////////////////////////////////////////////////////////////////
// NUM
//////////////////////////////////////////////////////////////////////////
extern double GetNumericalTimeCollocation(char*name);
extern int GetNumericsGaussPoints(int element_dimension);
extern double NUMCalcIterationError(double *new_iteration, double *old_iteration, double *reference, long length, int method);
#endif
