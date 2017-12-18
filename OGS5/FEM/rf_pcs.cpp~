/**************************************************************************
ROCKFLOW - Object: Process PCS
Programing:
02/2003 OK Implementation
  /2003 WW CRFProcessDeformation
11/2003 OK re-organized
07/2004 OK PCS2
02/2005 WW/OK Element Assemblier and output
12/2007 WW Classes of sparse matrix (jagged diagonal storage) and linear solver
           and parellelisation of them
02/2008 PCH OpenMP parallelization for Lis matrix solver
**************************************************************************/
#include "Configure.h"
#include "FEMEnums.h"

/*--------------------- MPI Parallel  -------------------*/
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
#include <mpi.h>
#endif
/*--------------------- MPI Parallel  -------------------*/

/*--------------------- OpenMP Parallel ------------------*/
#if defined(LIS)
#include "lis.h"
#include <omp.h>
#endif
/*--------------------- OpenMP Parallel ------------------*/

#include "makros.h"
// C
#ifndef __APPLE__
#include <malloc.h>
#endif

// C++
#include <cfloat>
#include <iostream>
#include <iomanip>                                //WW
//#include <algorithm> // header of transform. WW
// GEOLib
//#include "geo_ply.h"
/*------------------------------------------------------------------------*/
/* MshLib */
//#include "msh_elem.h"
//#include "msh_lib.h"
/*-----------------------------------------------------------------------*/
/* Objects */
#include "rf_pcs.h"
#include "pcs_dm.h"
#ifndef NEW_EQS                                   //WW. 07.11.2008
#include "solver.h"                               // ConfigRenumberProperties
#endif
#include "rf_st_new.h"                            // ST
//#include "rf_bc_new.h" // ST
//#include "rf_mmp_new.h" // MAT
#include "rf_ic_new.h"                            // IC
#include "fem_ele_std.h"                          // ELE
//#include "msh_lib.h" // ELE
//#include "rf_tim_new.h"
//#include "rf_out_new.h"
#include "rfmat_cp.h"
//#include "rf_mfp_new.h" // MFP
//#include "rf_num_new.h"
//#include "gs_project.h"
#include "rf_fct.h"
//#include "femlib.h"
#include "rf_node.h"

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/* Tools */
#ifndef NEW_EQS                                   //WW. 06.11.2008
#include "matrix.h"
#endif
#ifdef MFC                                        //WW
#include "rf_fluid_momentum.h"
#endif
/* Tools */
#include "mathlib.h"
//#include "files0.h"
//#include "par_ddc.h"
#include "tools.h"
//#include "rf_pcs.h"
#include "files0.h"
#ifdef GEM_REACT
// GEMS chemical solver
#include "rf_REACT_GEM.h"
REACT_GEM *m_vec_GEM;
#endif
#ifdef BRNS
// BRNS chemical solver
#include "rf_REACT_BRNS.h"
REACT_BRNS *m_vec_BRNS;
#endif

// New EQS
#ifdef NEW_EQS
#include "equation_class.h"
#endif
#include "problem.h"

using namespace std;

/*-------------------- ITPACKV    ---------------------------*/
extern void transM2toM6(void);
/*-------------------- ITPACKV    ---------------------------*/
/*-------------------- JAD    ---------------------------*/
extern void transM2toM5(void);
/*-------------------- JAD    ---------------------------*/
/*-----------------------------------------------------------------------*/
/* LOP */
//16.12.2008. WW #include "rf_apl.h" // Loop...
//16.12.2008. WW #include "loop_pcs.h"
extern VoidFuncVoid LOPCalcSecondaryVariables_USER;
//------------------------------------------------------------------------
// PCS
VoidXFuncVoidX PCSDestroyELEMatrices[PCS_NUMBER_MAX];
//------------------------------------------------------------------------
// Globals, to be checked
int pcs_no_fluid_phases = 0;
int pcs_no_components = 0;
bool pcs_monolithic_flow = false;
int dm_pcs_number = 0;
int pcs_deformation = 0;
int dm_number_of_primary_nvals = 2;
bool show_onces_adp = true;
bool show_onces_mod = true;
bool show_onces_mod_flow = true;
bool show_onces_density = true;
int memory_opt = 0;
int problem_2d_plane_dm;
int anz_nval = 0;
int anz_nval0 = 0;                                //WW
//
int size_eval=0;                                  //WW

NvalInfo *nval_data = NULL;
int anz_eval = 0;
EvalInfo *eval_data = NULL;
string project_title("New project");              //OK41

bool hasAnyProcessDeactivatedSubdomains = false;  //NW

//--------------------------------------------------------
// Coupling Flag. WW
bool T_Process = false;
bool H_Process = false;
bool M_Process = false;
bool RD_Process = false;
bool MH_Process = false;                          // MH monolithic scheme
bool MASS_TRANSPORT_Process = false;
bool FLUID_MOMENTUM_Process = false;
bool RANDOM_WALK_Process = false;
bool pcs_created = false;

namespace process{class CRFProcessDeformation;}
using process::CRFProcessDeformation;
using Mesh_Group::CNode;
using Mesh_Group::CElem;
using FiniteElement::ElementValue;
using Math_Group::vec;

#define noCHECK_EQS
#define noCHECK_ST_GROUP
#define noCHECK_BC_GROUP

extern size_t max_dim;                            //OK411 todo

//////////////////////////////////////////////////////////////////////////
// PCS vector
//////////////////////////////////////////////////////////////////////////
// It is better to have space between data type and data name. WW
vector<LINEAR_SOLVER *> PCS_Solver;               //WW
vector<CRFProcess*> pcs_vector;
vector<double*> ele_val_vector;                   //PCH
//vector<string> ele_val_name_vector; // PCH
template <class T> T *resize(T *array, size_t old_size, size_t new_size);
//////////////////////////////////////////////////////////////////////////
// Construction / destruction
//////////////////////////////////////////////////////////////////////////

/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2003 OK Implementation
02/2005 WW Local elment assembly (all protected members)
last modified:
**************************************************************************/
CRFProcess::CRFProcess(void) :
_problem (NULL), fem(NULL), Memory_Type(0),
Write_Matrix(false), matrix_file(NULL),
                                                  //WW  //WW
WriteSourceNBC_RHS(0), ele_val_name_vector (std::vector<std::string>())
{
   TempArry = NULL;
   //SB:GS4  pcs_component_number=0; //SB: counter for transport components
   pcs_component_number = pcs_no_components -1;
   //----------------------------------------------------------------------
   // NUM
   pcs_num_name[0] = NULL;
   pcs_num_name[1] = NULL;
   pcs_nonlinear_iterations = 1;
   pcs_nonlinear_iteration_tolerance = 1.0e8;
   pcs_coupling_iterations = 1;
   pcs_sol_name = NULL;
   m_num = NULL;
   cpl_type_name = "PARTITIONED";                 //OK
   num_type_name = "FEM";                         //OK
   rwpt_app = 0;                                  // PCH Application types for RWPT such as Cell Dynamics, Crypto, etc.
#ifndef NEW_EQS                                //WW 07.11.2008
   eqs = NULL;                                    //WW
#endif
   dof = 1;                                       //WW
   //----------------------------------------------------------------------
   // ELE
   pcs_number_of_evals = 0;
   NumDeactivated_SubDomains = 0;
   Deactivated_SubDomain = NULL;
   //----------------------------------------------------------------------
   //
   mobile_nodes_flag = -1;
   //----------------------------------------------------------------------
   // USER
   PCSSetIC_USER = NULL;
   //----------------------------------------------------------------------
   // TIM
   tim_type_name = "TRANSIENT";                   //OK
   time_unit_factor = 1.0;
   timebuffer = 1.0e-5;                           //WW
   //_pcs_type_name.empty();
   adaption = false;                              //HS 03.2008
   //----------------------------------------------------------------------
   // CPL
   for(int i=0; i<10; i++)
      Shift[i] = 0;
   selected = true;                               //OK
   // MSH OK
   m_msh = NULL;
   // Reload solutions
   reload=-1;
   nwrite_restart=1;                              // kg44 write every timestep is default
   pcs_nval_data = NULL;
   pcs_eval_data = NULL;
   non_linear = false;                            //OK/CMCD
   cal_integration_point_value = false;           //WW
   continuum = 0;
   //adaption = false; JOD removed
   compute_domain_face_normal = false;            //WW
   use_velocities_for_transport = false;
   //
   additioanl2ndvar_print = -1;                   //WW
   flow_pcs_type=0;                               //CB default: liquid flow, Sat = 1
   //----------------------------------------------------------------------
   m_bCheck = false;                              //OK
   m_bCheckOBJ = false;                           //OK
   m_bCheckNOD = false;                           //OK
   m_bCheckELE = false;                           //OK
   m_bCheckEQS = false;                           //OK
   //
   write_boundary_condition = false;              //15.01.2008. WW
   accepted = true;                               //25.08.2008. WW
   accept_steps = 0;                              //27.08.1008. WW
   reject_steps = 0;                              //27.08.1008. WW
   ML_Cap = 0;                                    // 23.01.2009 PCH
   PartialPS = 0;                                 // 16.02 2009 PCH

#ifdef USE_MPI                                 //WW
   cpu_time_assembly = 0;
#endif
   // New equation and solver WW
#ifdef NEW_EQS
   eqs_new = NULL;
   configured_in_nonlinearloop = false;
#endif
   flag_couple_GEMS = 0;                          // 11.2009 HS
}


void CRFProcess::setProblemObjectPointer (Problem* problem)
{
   _problem = problem;
}


const Problem* CRFProcess::getProblemObjectPointer () const
{
   return _problem;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2003 OK Implementation
01/2005 WW Local assemblier as a member
last modified:
**************************************************************************/
CRFProcess::~CRFProcess(void)
{
   long i;
   //----------------------------------------------------------------------
   // Finite element
   if(fem) delete fem;                            //WW
   fem = NULL;
   //----------------------------------------------------------------------
   // ELE: Element matrices
   ElementMatrix *eleMatrix = NULL;
   ElementValue* gp_ele = NULL;
   if(Ele_Matrices.size()>0)
   {
      for (i=0;i<(long)Ele_Matrices.size();i++)
      {
         eleMatrix = Ele_Matrices[i];
         delete eleMatrix;
         eleMatrix = NULL;
      }
      Ele_Matrices.clear();
   }
   //----------------------------------------------------------------------
   // ELE: Element Gauss point values
   if(ele_gp_value.size()>0)
   {
      for(i=0;i<(long)ele_gp_value.size();i++)
      {
         gp_ele = ele_gp_value[i];
         delete gp_ele;
         gp_ele = NULL;
      }
      ele_gp_value.clear();
   }
   //----------------------------------------------------------------------
   // OUT: Matrix output
   if (matrix_file)
   {
      matrix_file->close();
      delete matrix_file;
   }
   //----------------------------------------------------------------------
   // NOD: Release memory of node values
   for(i=0;i<(int)nod_val_vector.size();i++)
   {
      delete [] nod_val_vector[i];                // Add []. WW
      nod_val_vector[i] = NULL;
   }
   nod_val_vector.clear();
   //----------------------------------------------------------------------
   // ST:
   CNodeValue* m_nod_val = NULL;
   for(i=0;i<(int)st_node_value.size();i++)
   {
      m_nod_val = st_node_value[i];
      //OK delete st_node_value[i];
      //OK st_node_value[i] = NULL;
      if(m_nod_val->check_me)                     //OK
      {
         m_nod_val->check_me = false;
         delete m_nod_val;
         m_nod_val = NULL;
      }
   }
   st_node_value.clear();
   //----------------------------------------------------------------------
   for(i=0; i<(int)bc_node_value.size(); i++)
   {
      delete bc_node_value[i];
      bc_node_value[i] = NULL;
   }
   bc_node_value.clear();
   //----------------------------------------------------------------------
   //_pcs_type_name.clear();
   //----------------------------------------------------------------------
   // CON
   continuum_vector.clear();
   //Haibing 13112006------------------------------------------------------
   for(i=0;i<(long)ele_val_vector.size();i++)
      delete[] ele_val_vector[i];
   ele_val_vector.clear();
   //----------------------------------------------------------------------
   if(Deactivated_SubDomain)                      //05.09.2007 WW
   {
      delete [] Deactivated_SubDomain;
      Deactivated_SubDomain = NULL;
   }
}


/**************************************************************************
FEMLib-Method:
Task:    Gauss point values for CFEMSH
Programing:
08/2005 WW Implementation
**************************************************************************/
void CRFProcess::AllocateMemGPoint()
{
   //	if (_pcs_type_name.find("FLOW") == 0)
   //		return;
   const size_t mesh_ele_vector_size (m_msh->ele_vector.size());
   for (size_t i = 0; i < mesh_ele_vector_size; i++)
   {
      ele_gp_value.push_back(new ElementValue(this, m_msh->ele_vector[i]));
   }
}


/**************************************************************************
FEMLib-Method:
Task:    This function is a part of the monolithic scheme
         and it is used to assign pcs name to IC, ST, BC, TIM and OUT. object
Programing:
07/2005 WW Implementation
10/2010 TF cody style improvements
**************************************************************************/
void CRFProcess::SetOBJNames()
{
   // IC
   const size_t ic_vector_size(ic_vector.size());
   for (size_t i = 0; i < ic_vector_size; i++)
   {
      ic_vector[i]->setProcessType(this->getProcessType());
   }

   // BC
   std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
   while (p_bc != bc_list.end())
   {
      (*p_bc)->setProcessType(this->getProcessType());
      ++p_bc;
   }

   // ST
   const size_t st_vector_size (st_vector.size());
   for (size_t i = 0; i < st_vector_size; i++)
   {
      st_vector[i]->setProcessType(this->getProcessType());
   }

   // TIM
   const size_t time_vector_size (time_vector.size());
   for (size_t i = 0; i < time_vector_size; i++)
   {
      //		Tim = time_vector[i];
      //		Tim->pcs_type_name = _pcs_type_name;
      time_vector[i]->pcs_type_name = convertProcessTypeToString (this->getProcessType());
   }

   // OUT
   // OK4216
   const size_t out_vector_size (out_vector.size());
   for (size_t i = 0; i < out_vector_size; i++)
   {
      //		m_out = out_vector[i];
      //		m_out->_pcs_type_name = _pcs_type_name;
      out_vector[i]->setProcessType(this->getProcessType());
      //m_out->pcs_pv_name = pcs_primary_function_name[0];//CMCD
      //string temp = pcs_primary_function_name[0];
   }
}


/**************************************************************************
PCSLib-Method:
10/2002 OK Implementation
04/2004 WW Modification for 3D problems
02/2005 WW New fem calculator
04/2005 OK MSHCreateNOD2ELERelations
07/2005 WW Geometry element objects
02/2006 WW Removed memory leaking
01/2006 YD MMP for each PCS
04/2006 WW Unique linear solver for all processes if they share the same mesh
06/2006 WW Rearrange incorporation of BC and ST. Set BC and ST for domain decomposition
**************************************************************************/
void CRFProcess::Create()
{
   // we need the string representation of process type at some points
   std::string pcs_type_name(
      convertProcessTypeToString(this->getProcessType()));

   // Element matrix output. WW
   if (Write_Matrix)
   {
      std::cout << "->Write Matrix" << '\n';
#if defined(USE_MPI)
      char stro[32];
      sprintf(stro, "%d",myrank);
      string m_file_name = FileName +"_"+pcs_type_name+(string)stro+"_element_matrix.txt";
#else
      std::string m_file_name = FileName + "_" + pcs_type_name
         + "_element_matrix.txt";
#endif
      matrix_file = new std::fstream(m_file_name.c_str(), ios::trunc | ios::out);
      if (!matrix_file->good())
         std::cout << "Warning in GlobalAssembly: Matrix files are not found"
            << std::endl;
   }
   //----------------------------------------------------------------------------
   if (m_msh)                                     //OK->MB please shift to Config()
   {
      //		if (_pcs_type_name.compare("GROUNDWATER_FLOW") == 0)
      if (this->getProcessType() == GROUNDWATER_FLOW)
         MSHDefineMobile(this);
   }
   //----------------------------------------------------------------------------
   int DOF = GetPrimaryVNumber();                 //OK should be PCS member variable
   //----------------------------------------------------------------------------
   // MMP - create mmp groups for each process //YD
   std::cout << "->Create MMP" << std::endl;
   CMediumPropertiesGroup *m_mmp_group = NULL;
   int continua = 1;                              //WW
   if (RD_Process)
      continua = 2;
   m_mmp_group = MMPGetGroup(pcs_type_name);
   for (int i = 0; i < continua; i++)
   {
      if (!m_mmp_group)
      {
         m_mmp_group = new CMediumPropertiesGroup();
         m_mmp_group->pcs_type_name = pcs_type_name;
         m_mmp_group->Set(this);
         mmp_group_list.push_back(m_mmp_group);
      }
   }
   m_mmp_group = NULL;
   //----------------------------------------------------------------------------
   // NUM_NEW
   std::cout << "->Create NUM" << std::endl;
   //	if (pcs_type_name.compare("RANDOM_WALK")) { // PCH RWPT does not need this.
   if (this->getProcessType() != RANDOM_WALK)     // PCH RWPT does not need this.
   {
      CNumerics *m_num_tmp = NULL;
      size_t no_numerics(num_vector.size());
      for (size_t i = 0; i < no_numerics; i++)
      {
         m_num_tmp = num_vector[i];

         if ((pcs_type_name.compare(m_num_tmp->pcs_type_name) == 0)
            || (m_num_tmp->pcs_type_name.compare(
            pcs_primary_function_name[0]) == 0))
         {
            m_num = m_num_tmp;
            break;
         }
      }
   }
   if (!m_num)
   {
      std::cout << "Warning in CRFProcess::Create() - no NUM data" << std::endl;
      m_num = new CNumerics(pcs_type_name);       //OK
      //		m_num = m_num_tmp;
   }
   else
   {
      pcs_nonlinear_iterations = m_num->nls_max_iterations;
      pcs_nonlinear_iteration_tolerance = m_num->nls_error_tolerance;
      if (pcs_nonlinear_iterations > 1)           //WW
         non_linear = true;
   }
   //----------------------------------------------------------------------------
   // EQS - create equation system
   //WW CreateEQS();
   std::cout << "->Create EQS" << '\n';
#ifdef NEW_EQS
   size_t k;
   for(k=0; k<fem_msh_vector.size(); k++)
   {
      if(m_msh==fem_msh_vector[k])
         break;
   }
   if(type==4||type==41)
      eqs_new = EQS_Vector[2*k+1];
   else
      eqs_new = EQS_Vector[2*k];
#else
   //WW  phase=1;
   // create EQS
   if (type == 4)
   {
      eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF, DOF
         * m_msh->GetNodesNumber(true));
      InitializeLinearSolver(eqs, m_num);
      PCS_Solver.push_back(eqs);                  //WW
   }
   else if (type == 41)
   {
      if (num_type_name.find("EXCAVATION") != string::npos)
         eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF - 1, DOF
            * m_msh->GetNodesNumber(true));
      else
         eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF,
            (DOF - 1) * m_msh->GetNodesNumber(true)
            + m_msh->GetNodesNumber(false));
      InitializeLinearSolver(eqs, m_num);
      PCS_Solver.push_back(eqs);                  //WW
   }
   else
   {
      // If there is a solver exsiting. WW
      CRFProcess *m_pcs = NULL;
      for (size_t i = 0; i < pcs_vector.size(); i++)
      {
         m_pcs = pcs_vector[i];
         if (m_pcs && m_pcs->eqs)
         {
            //				if (m_pcs->_pcs_type_name.find("DEFORMATION") == string::npos)
            if (!isDeformationProcess(m_pcs->getProcessType()))
               break;
         }
      }
      // If unique mesh
      if (m_pcs && m_pcs->eqs && (fem_msh_vector.size() == 1))
         eqs = m_pcs->eqs;
      else
      {
         eqs = CreateLinearSolver(m_num->ls_storage_method,
            m_msh->GetNodesNumber(false) * DOF);
         InitializeLinearSolver(eqs, m_num);
         PCS_Solver.push_back(eqs);
      }
   }
#endif                                         // If NEW_EQS
   // Set solver properties: EQS<->SOL
   // Internen Speicher allokieren
   // Speicher initialisieren

   //----------------------------------------------------------------------------
   // Time unit factor //WW
   std::cout << "->Create TIM" << '\n';
   //CTimeDiscretization* Tim = TIMGet(_pcs_type_name);
   Tim = TIMGet(pcs_type_name);
   if (!Tim)
   {
                                                  //21.08.2008. WW
      Tim = new CTimeDiscretization(*time_vector[0], pcs_type_name);
      time_vector.push_back(Tim);                 //21.08.2008. WW
   }
   if (Tim->time_unit.find("MINUTE") != std::string::npos)
      time_unit_factor = 60.0;
   else if (Tim->time_unit.find("HOUR") != std::string::npos)
      time_unit_factor = 3600.0;
   else if (Tim->time_unit.find("DAY") != std::string::npos)
      time_unit_factor = 86400.0;
   else if (Tim->time_unit.find("MONTH") != std::string::npos)
      time_unit_factor = 2592000.0;
   else if (Tim->time_unit.find("YEAR") != std::string::npos)
      time_unit_factor = 31536000;

   //
   if (type == 4 || type == 41)
      m_msh->SwitchOnQuadraticNodes(true);
   else
      m_msh->SwitchOnQuadraticNodes(false);
   //
   if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC")
      != string::npos)                            //WW
   {
      setBC_danymic_problems();
      setST_danymic_problems();
   }
   else
   {
      // BC - create BC groups for each process
      cout << "->Create BC" << '\n';
      CBoundaryConditionsGroup *m_bc_group = NULL;
      for (int i = 0; i < DOF; i++)
      {
         //OKm_bc_group = BCGetGroup(_pcs_type_name,pcs_primary_function_name[i]);
         //OKif(!m_bc_group){
         BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
         m_bc_group = new CBoundaryConditionsGroup();
                                                  //OK
         m_bc_group->setProcessTypeName(pcs_type_name);
         m_bc_group->setProcessPrimaryVariableName(
            pcs_primary_function_name[i]);        //OK
         m_bc_group->Set(this, Shift[i]);
         bc_group_list.push_back(m_bc_group);     //Useless, to be removed. WW
         m_bc_group = NULL;
         //OK}
      }
      if (bc_node_value.size() < 1)               //WW
         cout << "Warning: no boundary conditions specified for "
            << pcs_type_name << endl;

      // ST - create ST groups for each process
      cout << "->Create ST" << '\n';
      CSourceTermGroup *m_st_group = NULL;

      if (WriteSourceNBC_RHS == 2)                // Read from file
         ReadRHS_of_ST_NeumannBC();
      else                                        // WW // Calculate directly
      {
         for (int i = 0; i < DOF; i++)
         {
            //OK m_st_group = m_st_group->Get(pcs_primary_function_name[i]);
            m_st_group = STGetGroup(pcs_type_name,
               pcs_primary_function_name[i]);
            if (!m_st_group)
            {
               m_st_group = new CSourceTermGroup();
                                                  //OK
               m_st_group->pcs_type_name = pcs_type_name;
                                                  //OK
               m_st_group->pcs_pv_name = pcs_primary_function_name[i];
               m_st_group->Set(this, Shift[i]);
                                                  //Useless, to be removed. WW
               st_group_list.push_back(m_st_group);
            }
         }
         if (WriteSourceNBC_RHS == 1)             // WW
            WriteRHS_of_ST_NeumannBC();
      }
      m_st_group = NULL;
   }
   // Write BC/ST nodes for vsualization.WW
   if (write_boundary_condition && WriteSourceNBC_RHS != 2)
      WriteBC();

   // ELE - config and create element values
   cout << "->Config ELE values" << '\n';
   AllocateMemGPoint();

   // ELE - config element matrices
   // NOD - config and create node values
   cout << "->Config NOD values" << '\n';
   double* nod_values = NULL;
   double* ele_values = NULL;                     // PCH

   number_of_nvals = 2* DOF + pcs_number_of_secondary_nvals;
   for (int i = 0; i < pcs_number_of_primary_nvals; i++)
   {
                                                  // new time
      nod_val_name_vector.push_back(pcs_primary_function_name[i]);
                                                  // old time //need this MB!
      nod_val_name_vector.push_back(pcs_primary_function_name[i]);
   }
   for (int i = 0; i < pcs_number_of_secondary_nvals; i++)
                                                  // new time
      nod_val_name_vector.push_back(pcs_secondary_function_name[i]);
   //
   long m_msh_nod_vector_size = m_msh->NodesNumber_Quadratic;
   for (long j = 0; j < m_msh_nod_vector_size; j++)
   {
      nod_values = new double[number_of_nvals];
      for (int i = 0; i < number_of_nvals; i++)
         nod_values[i] = 0.0;
      nod_val_vector.push_back(nod_values);
   }
   // Create element values - PCH
   int number_of_evals = 2* pcs_number_of_evals ; //PCH, increase memory
   if (number_of_evals > 0)                       // WW added this "if" condition
   {
      for (int i = 0; i < pcs_number_of_evals; i++)
      {
                                                  // new time
         ele_val_name_vector.push_back(pcs_eval_name[i]);
                                                  // old time
         ele_val_name_vector.push_back(pcs_eval_name[i]);
      }
      size_t m_msh_ele_vector_size(m_msh->ele_vector.size());
      if (ele_val_vector.size() == 0)
      {
         for (size_t j = 0; j < m_msh_ele_vector_size; j++)
         {
            ele_values = new double[number_of_evals];
            size_eval += number_of_evals;         //WW
            for (int i = 0; i < number_of_evals; i++)
               ele_values[i] = 0.0;
            ele_val_vector.push_back(ele_values);
         }
      }
      else
      {
         for (size_t j = 0; j < m_msh_ele_vector_size; j++)
         {
            ele_values = ele_val_vector[j];
            ele_values = resize(ele_values, size_eval, size_eval
               + number_of_evals);
            size_eval += number_of_evals;
            ele_val_vector[j] = ele_values;
         }
      }
   }
   //
   if (reload >= 2 && type != 4 && type != 41)
   {
                                                  // PCH
      cout << "Reloading the primary variables... " << endl;
      ReadSolution();                             //WW
   }

   if (reload < 2)                                // PCH: If reload is set, no need to have ICs
   {
      // IC
      cout << "->Assign IC" << '\n';
      SetIC();
   }
   else
   {
      // Bypassing IC
      cout << "RELOAD is set to be " << reload << ". So, bypassing IC's"
         << endl;
   }

   if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC")
      != string::npos)                            //WW
      setIC_danymic_problems();

   // Keep all local matrices in the memory
   if (type != 55)                                //Not for fluid momentum. WW
   {
      if (Memory_Type != 0)
         AllocateLocalMatrixMemory();
      if (type == 4 || type == 41)
      {
         // Set initialization function
         //      CRFProcessDeformation *dm_pcs = (CRFProcessDeformation *) (this);
         CRFProcessDeformation *dm_pcs =
            static_cast<CRFProcessDeformation *> (this);
         dm_pcs->Initialization();
      } else                                      // Initialize FEM calculator
      {
         int Axisymm = 1;                         // ani-axisymmetry
         if (m_msh->isAxisymmetry())
            Axisymm = -1;                         // Axisymmetry is true
         fem = new CFiniteElementStd(this, Axisymm
            * m_msh->GetCoordinateFlag());
         fem->SetGaussPointNumber(m_num->ele_gauss_points);
      }
   }

   // Initialize the system equations
   if (PCSSetIC_USER)
      PCSSetIC_USER(pcs_type_number);

   if (compute_domain_face_normal)                //WW
      m_msh->FaceNormal();
}


/**************************************************************************
FEMLib-Method:
Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
Programing:
12/2005 WW
03/2006 WW Write as acsi
04/2006 WW
last modified:
**************************************************************************/
inline void CRFProcess::WriteRHS_of_ST_NeumannBC()
{
   std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));
   std::string m_file_name = FileName + "_" + pcs_type_name + "_ST_RHS.asc";
   std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
   if (!os.good())
   {
      cout << "Failure to open file: " << m_file_name << endl;
      abort();
   }

   os << "$PCS_TYPE  " << endl;

   os << pcs_type_name << endl;
   os << "geo_node_number  ";
   os << "msh_node_number  ";
   os << "CurveIndex ";
   os << "node_value ";
   os << endl;
   os.setf(std::ios::scientific, std::ios::floatfield);
   os.precision(14);
   const size_t st_node_value_size(st_node_value.size());
   os << st_node_value_size << endl;
   for (size_t i = 0; i < st_node_value_size; i++)
      st_node_value[i]->Write(os);
   os.close();
}


/**************************************************************************
FEMLib-Method:
Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
Programing:
03/2006 WW
last modified: 04/2006
**************************************************************************/
inline void CRFProcess::ReadRHS_of_ST_NeumannBC()
{
   std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));
   std::string m_file_name = FileName + "_" + pcs_type_name + "_ST_RHS.asc";
   std::ifstream is(m_file_name.c_str(), std::ios::in);
   if (!is.good())
   {
      cout << "File " << m_file_name << " is not found" << endl;
      abort();
   }

   std::string s_buffer;
   getline(is, s_buffer);
   getline(is, s_buffer);
   getline(is, s_buffer);
   size_t size;
   is >> size >> ws;
   st_node_value.clear();
   for (size_t i = 0; i < size; i++)
   {
      CNodeValue *cnodev = new CNodeValue();
      cnodev->Read(is);
      st_node_value.push_back(cnodev);
   }
   is.close();
}


/**************************************************************************
FEMLib-Method:
Task: Write the solution
Programing:
04/2006 WW
last modified:
**************************************************************************/
void CRFProcess:: WriteSolution()
{
   if ( reload==2 || reload<=0 ) return;
                                                  //kg44 write out only between nwrite_restart timesteps
   if ( ( aktueller_zeitschritt % nwrite_restart  ) > 0 ) return;

   std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));
   std::string m_file_name = FileName +"_"+pcs_type_name+"_"+pcs_primary_function_name[0]+"_primary_value.asc";
   std::ofstream os ( m_file_name.c_str(), ios::trunc|ios::out );
   if (!os.good() )
   {
      cout << "Failure to open file: "<<m_file_name << endl;
      abort();
   }

   os.precision(15);                              // 15 digits accuracy seems enough? more fields are filled up with random numbers!
   os.setf(std::ios_base::scientific,std::ios_base::floatfield);

   int j;
   int *idx(new int [2*pcs_number_of_primary_nvals]);
   for ( j=0; j<pcs_number_of_primary_nvals; j++ )
   {
      idx[j] = GetNodeValueIndex ( pcs_primary_function_name[j] );
      idx[j+pcs_number_of_primary_nvals] = idx[j]+1;
   }
   for (long i=0; i<m_msh->GetNodesNumber ( false ); i++ )
   {
      for ( j=0; j<2*pcs_number_of_primary_nvals; j++ )
         os<<GetNodeValue ( i,idx[j] ) <<"  ";
      os<<endl;
   }
   os.close();
   cout << "Write solutions for timestep " << aktueller_zeitschritt << " into file " << m_file_name << endl;
   delete [] idx;
}


/**************************************************************************
FEMLib-Method:
Task: Write the solution
Programing:
04/2006 WW
last modified:
**************************************************************************/
void CRFProcess:: ReadSolution()
{
   std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));
   std::string m_file_name = FileName +"_"+pcs_type_name+"_"+pcs_primary_function_name[0]+"_primary_value.asc";
   std::ifstream is ( m_file_name.c_str(), ios::in );
   if (!is.good())
   {
      cout << "Failure to open file: "<<m_file_name << endl;
      abort();
   }
   int j;

   int *idx (new int [2*pcs_number_of_primary_nvals]);
   double *val (new double [2*pcs_number_of_primary_nvals]);

   for (j=0; j<pcs_number_of_primary_nvals; j++ )
   {
      idx[j] = GetNodeValueIndex ( pcs_primary_function_name[j] );
      idx[j+pcs_number_of_primary_nvals] = idx[j]+1;
   }
   for (long i=0; i<m_msh->GetNodesNumber ( false ); i++ )
   {
      for (j=0; j<2*pcs_number_of_primary_nvals; j++ )
         is>>val[j];
      is>>ws;
      for (j=0; j<2*pcs_number_of_primary_nvals; j++ )
         SetNodeValue ( i,idx[j], val[j] );
   }
   is.close();
   delete [] idx ;
   delete [] val ;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 WW Set coupling data
last modified:
**************************************************************************/
void CRFProcess::setIC_danymic_problems()
{
   const char *function_name[7];
   int i, j, nv;
   nv = 0;
   if(max_dim==1)                                 // 2D
   {
      nv = 5;
      function_name[0] = "DISPLACEMENT_X1";
      function_name[1] = "DISPLACEMENT_Y1";
      function_name[2] = "VELOCITY_DM_X";
      function_name[3] = "VELOCITY_DM_Y";
      function_name[4] = "PRESSURE1";
   }
   else                                           //3D
   {
      nv =7;
      function_name[0] = "DISPLACEMENT_X1";
      function_name[1] = "DISPLACEMENT_Y1";
      function_name[2] = "DISPLACEMENT_Z1";
      function_name[3] = "VELOCITY_DM_X";
      function_name[4] = "VELOCITY_DM_Y";
      function_name[5] = "VELOCITY_DM_Z";
      function_name[6] = "PRESSURE1";
   }

   CInitialCondition *m_ic = NULL;
   long no_ics =(long)ic_vector.size();
   int nidx;
   for(i=0;i<nv;i++)
   {
      nidx = GetNodeValueIndex(function_name[i]);
      for(j=0;j<no_ics;j++)
      {
         m_ic = ic_vector[j];
         if(m_ic->getProcessPrimaryVariable() == convertPrimaryVariable(function_name[i]))
            m_ic->Set(nidx);
      }
   }

}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 WW Set coupling data
last modified:
**************************************************************************/
void CRFProcess::setST_danymic_problems()
{
   const char *function_name[7];
   size_t nv = 0;
   if (max_dim == 1)                              // 2D
   {
      nv = 5;
      function_name[0] = "DISPLACEMENT_X1";
      function_name[1] = "DISPLACEMENT_Y1";
      function_name[2] = "VELOCITY_DM_X";
      function_name[3] = "VELOCITY_DM_Y";
      function_name[4] = "PRESSURE1";
   }                                              //3D
   else
   {
      nv = 7;
      function_name[0] = "DISPLACEMENT_X1";
      function_name[1] = "DISPLACEMENT_Y1";
      function_name[2] = "DISPLACEMENT_Z1";
      function_name[3] = "VELOCITY_DM_X";
      function_name[4] = "VELOCITY_DM_Y";
      function_name[5] = "VELOCITY_DM_Z";
      function_name[6] = "PRESSURE1";
   }

   // ST - create ST groups for each process
   CSourceTermGroup *m_st_group = NULL;
   std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));
   for (size_t i = 0; i < nv; i++)
   {
      m_st_group = STGetGroup(pcs_type_name, function_name[i]);
      if (!m_st_group)
      {
         m_st_group = new CSourceTermGroup();
         m_st_group->pcs_type_name = pcs_type_name;
         m_st_group->pcs_pv_name = function_name[i];
         m_st_group->Set(this, Shift[i], function_name[i]);
         st_group_list.push_back(m_st_group);     //Useless, to be removed. WW
      }
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 WW Set coupling data
last modified:
**************************************************************************/
void CRFProcess::setBC_danymic_problems()
{
   const char *function_name[7];
   size_t nv = 0;
   if (max_dim == 1)                              // 2D
   {
      nv = 5;
      function_name[0] = "DISPLACEMENT_X1";
      function_name[1] = "DISPLACEMENT_Y1";
      function_name[2] = "VELOCITY_DM_X";
      function_name[3] = "VELOCITY_DM_Y";
      function_name[4] = "PRESSURE1";
   }                                              //3D
   else
   {
      nv = 7;
      function_name[0] = "DISPLACEMENT_X1";
      function_name[1] = "DISPLACEMENT_Y1";
      function_name[2] = "DISPLACEMENT_Z1";
      function_name[3] = "VELOCITY_DM_X";
      function_name[4] = "VELOCITY_DM_Y";
      function_name[5] = "VELOCITY_DM_Z";
      function_name[6] = "PRESSURE1";
   }

   cout << "->Create BC" << '\n';
   CBoundaryConditionsGroup *m_bc_group = NULL;
   std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));
   for (size_t i = 0; i < nv; i++)
   {
      BCGroupDelete(pcs_type_name, function_name[i]);
      m_bc_group = new CBoundaryConditionsGroup();
                                                  //OK
      m_bc_group->setProcessTypeName(pcs_type_name);
                                                  //OK
      m_bc_group->setProcessPrimaryVariableName(function_name[i]);
      m_bc_group->Set(this, Shift[i], function_name[i]);
      bc_group_list.push_back(m_bc_group);        //Useless, to be removed. WW
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
02/2005 WW Set coupling data
last modified:
**************************************************************************/
void CRFProcess::ConfigureCouplingForLocalAssemblier()
{
   bool Dyn = false;
   if(pcs_type_name_vector.size()&&pcs_type_name_vector[0].find("DYNAMIC")!=string::npos)
      Dyn = true;
   if(fem)
      fem->ConfigureCoupling(this, Shift, Dyn);
}


/**************************************************************************
FEMLib-Method:
06/2003 OK Implementation
        WW 2nd version, PCS_Solver
**************************************************************************/
void PCSDestroyAllProcesses(void)
{
   CRFProcess *m_process = NULL;
   long i;
   int j;
   //----------------------------------------------------------------------
   // SOLver
#ifdef NEW_EQS                                 //WW
#if defined(USE_MPI)
   for(j=0;j<(int)EQS_Vector.size();j+=2)         //WW
#else
      for(j=0;j<(int)EQS_Vector.size();j++)       //WW
#endif
   {
      if(EQS_Vector[j])  delete EQS_Vector[j];
      EQS_Vector[j] = NULL;
#if defined(USE_MPI)
      EQS_Vector[j+1] = NULL;
#endif
   }
#else
   // SOLDelete()
   LINEAR_SOLVER *eqs;
   for(j=0;j<(int)PCS_Solver.size();j++)
   {
      eqs = PCS_Solver[j];
      if(eqs->unknown_vector_indeces)
         eqs->unknown_vector_indeces = \
            (int*) Free(eqs->unknown_vector_indeces);
      if(eqs->unknown_node_numbers)
         eqs->unknown_node_numbers = \
            (long*) Free(eqs->unknown_node_numbers);
      if(eqs->unknown_update_methods)
         eqs->unknown_update_methods = \
            (int*) Free(eqs->unknown_update_methods);
      eqs = DestroyLinearSolver(eqs);
   }
   PCS_Solver.clear();                            //WW
#endif
   //----------------------------------------------------------------------
   // PCS
   for(j=0;j<(int)pcs_vector.size();j++)
   {
      m_process = pcs_vector[j];
#ifdef USE_MPI                              //WW
      //  if(myrank==0)
      m_process->Print_CPU_time_byAssembly();
#endif
      if(m_process->pcs_nval_data)
         m_process->pcs_nval_data = (PCS_NVAL_DATA *) Free(m_process->pcs_nval_data);
      if(m_process->pcs_eval_data)
         m_process->pcs_eval_data = (PCS_EVAL_DATA *) Free(m_process->pcs_eval_data);
#ifdef PCS_NOD
      for (i=0;i<NodeListSize();i++)
      {
         k = GetNode(i);
         k->values[m_process->pcs_number] = (double *)Free(k->values[m_process->pcs_number]);
      }
#endif
      if(m_process->TempArry)                     //MX
         m_process->TempArry = (double *) Free(m_process->TempArry);
      delete(m_process);
   }
   pcs_vector.clear();                            //WW
   //----------------------------------------------------------------------
   // MSH
   for(i=0;i<(long)fem_msh_vector.size();i++)
   {
      if(fem_msh_vector[i]) delete fem_msh_vector[i];
      fem_msh_vector[i] = NULL;
   }
   fem_msh_vector.clear();
   //----------------------------------------------------------------------

   // DOM WW
#if defined(USE_MPI)
   //if(myrank==0)
   dom_vector[myrank]->PrintEQS_CPUtime();        //WW
#endif
   for(i=0;i<(long)dom_vector.size();i++)
   {
      if(dom_vector[i]) delete dom_vector[i];
      dom_vector[i] = NULL;
   }
   dom_vector.clear();
   //----------------------------------------------------------------------
   // ELE
   for(i=0;i<(long)ele_val_vector.size();i++)
      delete ele_val_vector[i];
   ele_val_vector.clear();
   //----------------------------------------------------------------------
   // IC ICDelete()
   for(i=0;i<(long)ic_vector.size();i++)
      delete ic_vector[i];
   ic_vector.clear();
   //----------------------------------------------------------------------
   MSPDelete();                                   //WW
   BCDelete();                                    //WW
   ICDelete();                                    //HS
   BCGroupDelete();                               //HS
   STDelete();                                    //WW
   STGroupsDelete();                              //HS
   GEOLIB_Clear_GeoLib_Data();                    //HS
   //......................................................................
   TIMDelete();                                   //OK
   OUTDelete();
   NUMDelete();
   MFPDelete();
   MSPDelete();
   MMPDelete();
   MMPGroupDelete();
   MCPDelete();
   //----------------------------------------------------------------------
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2004 OK Implementation
08/2004 WW Read the deformation process
           Check the comment key '//' in .pcs
last modified:
10/2010 TF changed process type handling from string to enum
**************************************************************************/
bool PCSRead(std::string file_base_name)
{
   //----------------------------------------------------------------------
   char line[MAX_ZEILE];
   int indexCh1a, indexCh2a;
   std::string CommentK("//");
   std::string line_string;
   std::string pcs_file_name;
   std::ios::pos_type position;
   //========================================================================
   // File handling
   pcs_file_name = file_base_name + PCS_FILE_EXTENSION;
   std::ifstream pcs_file(pcs_file_name.data(), ios::in);
   if (!pcs_file.good())
   {
      cout << "Warning: no PCS data *.pcs file is missing" << endl;
      return false;
   }

   // rewind the file
   pcs_file.clear();
   pcs_file.seekg(0, std::ios::beg);
   //========================================================================
   // Keyword loop
   std::cout << "PCSRead ... " << std::flush;
   while (!pcs_file.eof())
   {
      pcs_file.getline(line, MAX_ZEILE);
      line_string = line;
      line_string = GetUncommentedLine(line_string);
      if (line_string.find("#STOP") != string::npos)
         break;
      indexCh1a = (int) line_string.find_first_of(CommentK.c_str());
      indexCh2a = (int) line_string.find("#PROCESS");
      //----------------------------------------------------------------------
                                                  // keyword found
      if (indexCh2a > indexCh1a && (indexCh1a == -1))
      {
         CRFProcess *m_pcs = new CRFProcess();
         m_pcs->file_name_base = file_base_name;  //OK
         position = m_pcs->Read(&pcs_file);
         //..................................................................
         // pcs_type_name_vector
         if (m_pcs->pcs_type_name_vector.size() > 1)
         {
            string pname = m_pcs->pcs_type_name_vector[0]
               + m_pcs->pcs_type_name_vector[1];
            m_pcs->pcs_type_name_vector.pop_back();
            if (pname.find("FLOW") != string::npos && pname.find(
               "DEFORMATION") != string::npos)
            {
               //					m_pcs->_pcs_type_name = "DEFORMATION_FLOW";
               m_pcs->setProcessType (DEFORMATION_FLOW);
               MH_Process = true;                 // MH monolithic scheme
               if (pname.find("DYNAMIC") != string::npos)
                  m_pcs->pcs_type_name_vector[0] = "DYNAMIC";
            }
         }

         m_pcs->pcs_number = pcs_vector.size();
         //RelocateDeformationProcess(m_pcs);
         //			if (m_pcs->_pcs_type_name.find("DEFORMATION") != string::npos) { // TF
         if (isDeformationProcess(m_pcs->getProcessType()))
         {
            pcs_vector.push_back(m_pcs->CopyPCStoDM_PCS());
            pcs_vector[pcs_vector.size() - 1]->pcs_number = pcs_vector.size();
            delete m_pcs;
         } else
         pcs_vector.push_back(m_pcs);

         pcs_file.seekg(position, std::ios::beg);
      }                                           // keyword found
   }                                              // eof

   std::cout << "done, read " << pcs_vector.size() << " processes" << std::endl;

   return true;
}


/**************************************************************************
FEMLib-Method:
Task: Copy data to dm_pcs from PCS read function
Programing:
06/2007 OK/WW Implementation
10/2010 TF many improvements
**************************************************************************/
CRFProcess* CRFProcess::CopyPCStoDM_PCS()
{
   // Numerics
   if(num_type_name.compare("STRONG_DISCONTINUITY")==0)
      enhanced_strain_dm=1;

   CRFProcessDeformation *dm_pcs (new CRFProcessDeformation());
   dm_pcs->setProcessType (this->getProcessType());
   dm_pcs->pcs_type_name_vector.push_back(pcs_type_name_vector[0].data());
   dm_pcs->Write_Matrix = Write_Matrix;
   dm_pcs->WriteSourceNBC_RHS = WriteSourceNBC_RHS;
   dm_pcs->num_type_name = num_type_name;
   dm_pcs-> Memory_Type = Memory_Type;
   dm_pcs->NumDeactivated_SubDomains = NumDeactivated_SubDomains;
   dm_pcs->reload = reload;
                                                  //WW
   dm_pcs->write_boundary_condition =  write_boundary_condition;
   if(!dm_pcs->Deactivated_SubDomain)
      dm_pcs->Deactivated_SubDomain = new int[NumDeactivated_SubDomains];
   for(int i=0; i<NumDeactivated_SubDomains; i++)
      dm_pcs->Deactivated_SubDomain[i] = Deactivated_SubDomain[i];
   pcs_deformation = 1;
   //
   return dynamic_cast<CRFProcess *> (dm_pcs);
}


/**************************************************************************
FEMLib-Method:
Task: PCS read function
Programing:
06/2004 OK Implementation
08/2004 WW Read deformation process
11/2004 OK file streaming
12/2005 OK MSH_TYPE
01/2006 OK GEO_TYPE
**************************************************************************/
std::ios::pos_type CRFProcess::Read(std::ifstream *pcs_file)
{
   char line[MAX_ZEILE];
   string line_string;
   string CommentK("//");
   string hash("#");
   bool new_keyword = false;
   bool new_subkeyword = false;
   ios::pos_type position;
   ios::pos_type position_subkeyword;
   std::stringstream line_stream;
   saturation_switch = false;                     // JOD for Richards
   //----------------------------------------------------------------------
   while (!new_keyword)
   {
      position = pcs_file->tellg();
      pcs_file->getline(line, MAX_ZEILE);
      line_string = line;
      if (line_string.find(hash) != string::npos)
      {
         new_keyword = true;
         break;
      }
      //....................................................................
      // WW Comment line
      if (line_string.find_first_of(CommentK.c_str()) != string::npos)
         return position;
      //SB check for comment sign ;
      line_string = GetUncommentedLine(line_string);
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$PCS_TYPE") != string::npos)
      {
         while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
         {
            position = pcs_file->tellg();
            line_string = GetLineFromFile1(pcs_file);
            if (line_string.find("#") != string::npos)
            {
               return position;
            }
            if (line_string.find("$") != string::npos)
            {
               new_subkeyword = true;
               break;
            }
            line_stream.str(line_string);
            std::string pcs_type_name;
            line_stream >> pcs_type_name;
            this->setProcessType (convertProcessType(pcs_type_name));
            line_stream.clear();

            //				if (_pcs_type_name.find("FLOW") != string::npos) {
            if (isFlowProcess (this->getProcessType()))
            {
               pcs_no_fluid_phases++;
               H_Process = true;
            }
            //				if (_pcs_type_name.compare("PS_GLOBAL") == 0) {
            if (this->getProcessType() == PS_GLOBAL)
            {
               H_Process = true;
            }
            //				if (_pcs_type_name.compare("FLUID_FLOW") == 0) {
            //					_pcs_type_name = "LIQUID_FLOW";
            //				}
            //				if (_pcs_type_name.find("DEFORMATION") != string::npos) {
            if (isDeformationProcess(this->getProcessType()))
            {
               M_Process = true;
               //if only "DEFORMATION", do not create solver for FLUID process. WW
               if (pcs_no_fluid_phases < 1)
                  pcs_no_fluid_phases = 1;
            }
            //				if (_pcs_type_name.compare("MASS_TRANSPORT") == 0) {
            if (this->getProcessType() == MASS_TRANSPORT)
            {
               H_Process = true;
               MASS_TRANSPORT_Process = true;
               pcs_no_components++;
               this->setProcessPrimaryVariable(CONCENTRATION);
            }
            //				if (_pcs_type_name.find("HEAT") != string::npos)
            if (this->getProcessType() == HEAT_TRANSPORT)
               T_Process = true;

            pcs_type_name_vector.push_back(pcs_type_name);

            //				if (_pcs_type_name.compare("FLUID_MOMENTUM") == 0) {
            if (this->getProcessType() == FLUID_MOMENTUM)
            {
               FLUID_MOMENTUM_Process = true;
            }
            //				if (_pcs_type_name.compare("RANDOM_WALK") == 0) {
            if (this->getProcessType() == RANDOM_WALK)
            {
               RANDOM_WALK_Process = true;
            }
         }
         /*
          *pcs_file >> _pcs_type_name;
          if(_pcs_type_name.find("FLOW")!=string::npos){
          pcs_no_fluid_phases++;
          }
          if(_pcs_type_name.compare("FLUID_FLOW")==0){
          _pcs_type_name = "LIQUID_FLOW";
          }
          if(_pcs_type_name.find("DEFORMATION")!=string::npos)
          {
          //if only "DEFORMATION", do not create solver for FLUID process. WW
         if(pcs_no_fluid_phases<1)
         pcs_no_fluid_phases = 1;
         }
         if(_pcs_type_name.compare("MASS_TRANSPORT")==0)
         pcs_no_components++;
         pcs_file->ignore(MAX_ZEILE,'\n');
         continue;
         */
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$NUM_TYPE") != string::npos)
      {
         *pcs_file >> num_type_name;
         pcs_file->ignore(MAX_ZEILE, '\n');
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$CPL_TYPE") != string::npos)
      {
         *pcs_file >> cpl_type_name;
         if (cpl_type_name.compare("MONOLITHIC") == 0)
         {
            pcs_monolithic_flow = true;
            pcs_deformation = 11;
         }
         pcs_file->ignore(MAX_ZEILE, '\n');
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$TIM_TYPE") != string::npos)
      {
         *pcs_file >> tim_type_name;
         pcs_file->ignore(MAX_ZEILE, '\n');
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$APP_TYPE") != string::npos)
      {
         *pcs_file >> rwpt_app;
         pcs_file->ignore(MAX_ZEILE, '\n');
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$PRIMARY_VARIABLE") != string::npos)
      {
         *pcs_file >> primary_variable_name;
         pcs_file->ignore(MAX_ZEILE, '\n');
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$ELEMENT_MATRIX_OUTPUT") != string::npos)
      {
         *pcs_file >> Write_Matrix;               //WW
         pcs_file->ignore(MAX_ZEILE, '\n');
         continue;
      }
      //....................................................................
                                                  //WW
      if (line_string.find("$BOUNDARY_CONDITION_OUTPUT") != string::npos)
      {
         write_boundary_condition = true;
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$ST_RHS") != string::npos)
      {
         *pcs_file >> WriteSourceNBC_RHS;         //WW
         pcs_file->ignore(MAX_ZEILE, '\n');
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$MEMORY_TYPE") != string::npos)
      {
         *pcs_file >> Memory_Type;                //WW
         pcs_file->ignore(MAX_ZEILE, '\n');
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$RELOAD") != string::npos)
      {
         *pcs_file >> reload;                     //WW
         if (reload == 1 || reload == 3)
            *pcs_file >> nwrite_restart;          //kg44 read number of timesteps between writing restart files
         pcs_file->ignore(MAX_ZEILE, '\n');
         continue;
      }
                                                  // subkeyword found
      if (line_string.find("$DEACTIVATED_SUBDOMAIN") != string::npos)
      {
                                                  //WW
         *pcs_file >> NumDeactivated_SubDomains >> ws;
         Deactivated_SubDomain = new int[NumDeactivated_SubDomains];
         for (int i = 0; i < NumDeactivated_SubDomains; i++)
            *pcs_file >> Deactivated_SubDomain[i] >> ws;
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$MSH_TYPE") != string::npos)
      {
         *pcs_file >> msh_type_name >> ws;
         continue;
      }
      //....................................................................
      //		if (line_string.find("$GEO_TYPE") != string::npos) { //OK
      //			*pcs_file >> geo_type >> geo_type_name >> ws;
      //			continue;
      //		}
      //
      //....................................................................
                                                  // subkeyword found
      if (line_string.find("$MEDIUM_TYPE") != string::npos)
      {
         while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
         {
            position_subkeyword = pcs_file->tellg();
            *pcs_file >> line_string;
            if (line_string.size() == 0)
               break;
            if (line_string.find("#") != string::npos)
            {
               new_keyword = true;
               break;
            }
            if (line_string.find("$") != string::npos)
            {
               new_subkeyword = true;
               break;
            }
            if (line_string.find("CONTINUUM") != string::npos)
            {
               *pcs_file >> line_string;
                                                  //WW
               double w_m = strtod(line_string.data(), NULL);
               continuum_vector.push_back(w_m);
                                                  //WW
               continuum_vector.push_back(1.0 - w_m);
               break;                             //WW
            }
            pcs_file->ignore(MAX_ZEILE, '\n');
         }
         continue;
      }
                                                  //OK
      if (line_string.find("$SATURATION_SWITCH") != string::npos)
      {
         saturation_switch = true;
         ;
         continue;
      }
                                                  //SB4900
      if (line_string.find("$USE_VELOCITIES_FOR_TRANSPORT") != string::npos)
      {
                                                  // Only for fluid momentum process
         if (this->getProcessType () == FLUID_MOMENTUM)
            use_velocities_for_transport = true;
         continue;
      }
      //....................................................................
      //....................................................................
   }
   //----------------------------------------------------------------------
   return position;
}


/**************************************************************************
FEMLib-Method:
01/2004 OK Implementation
08/2004 WW Read the deformation process
           Check the comment key '//' in .pcs
06/2009 OK Write only if existing
**************************************************************************/
void PCSWrite(string file_base_name)
{
   if((int)pcs_vector.size()<1)
      return;
   //----------------------------------------------------------------------
   // File handling
   string pcs_file_name = file_base_name + PCS_FILE_EXTENSION;
   fstream pcs_file (pcs_file_name.data(),ios::trunc|ios::out);
   pcs_file.clear();
   //----------------------------------------------------------------------
   // PCS loop
   cout << "PCSWrite" << endl;
   CRFProcess* m_pcs = NULL;
   for(int i=0;i<(int)pcs_vector.size();i++)
   {
      m_pcs = pcs_vector[i];
      pcs_file << "#PROCESS" << endl;
      m_pcs->Write(&pcs_file);
   }
   //----------------------------------------------------------------------
   pcs_file << "#STOP" << endl;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2004 OK Implementation
12/2005 OK MSH_TYPE
last modified:
**************************************************************************/
void CRFProcess::Write(std::fstream *pcs_file)
{
   *pcs_file << " $PCS_TYPE" << endl;
   *pcs_file << "  " << convertProcessTypeToString(this->getProcessType()) << endl;

   *pcs_file << " $NUM_TYPE" << endl;
   *pcs_file << "  " << num_type_name << endl;

   *pcs_file << " $CPL_TYPE" << endl;
   *pcs_file << "  " << cpl_type_name << endl;

   *pcs_file << " $TIM_TYPE" << endl;
   *pcs_file << "  " << tim_type_name << endl;

   if (msh_type_name.size() > 0)
   {
      *pcs_file << " $MSH_TYPE" << endl;
      *pcs_file << "  " << msh_type_name << endl;
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
12/2004 OK Implementation
last modified:
10/2010 changed access to process type
**************************************************************************/
CRFProcess* PCSGet(const std::string &pcs_type_name)
{
   ProcessType pcs_type(convertProcessType(pcs_type_name));
   for (size_t i = 0; i < pcs_vector.size(); i++)
   {
      //		m_pcs = pcs_vector[i];
      //		if(m_pcs->pcs_type_name.compare(pcs_type_name)==0) { TF
      if (pcs_vector[i]->getProcessType() == pcs_type)
      {
         return pcs_vector[i];
      }
   }
   return NULL;
}


CRFProcess* PCSGet(ProcessType pcs_type)
{
   for (size_t i = 0; i < pcs_vector.size(); i++)
   {
      if (pcs_vector[i]->getProcessType() == pcs_type)
      {
         return pcs_vector[i];
      }
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
11/2008 TK New Version with Primary Variable Comparision
last modified:
10/2010 TF changed access to process type
**************************************************************************/
CRFProcess* PCSGetNew(const string &pcs_type_name, const string &primary_variable_name)
{
   CRFProcess *m_pcs_return = NULL;

   ProcessType pcs_type (convertProcessType (pcs_type_name));

   int matches = 0;
   for (size_t i = 0; i < pcs_vector.size(); i++)
   {
      CRFProcess *pcs = pcs_vector[i];
      //		if (pcs->pcs_type_name.compare(pcs_type_name) == 0) { TF
      if (pcs->getProcessType () == pcs_type)
      {
         for (size_t j = 0; j < pcs->GetPrimaryVNumber(); j++)
         {
            if (primary_variable_name.compare(pcs->GetPrimaryVName(j)) == 0)
            {
               m_pcs_return = pcs;
               matches++;
               if (matches > 1)
                  return NULL;
            }
         }
      }
   }
   if (matches == 0)
      return NULL;
   else return m_pcs_return;
}


//////////////////////////////////////////////////////////////////////////
// Access
//////////////////////////////////////////////////////////////////////////

// OK->SB please try Get function
CRFProcess *CRFProcess::GetProcessByFunctionName(char *name)
{
   CRFProcess *m_process = NULL;
   /* Tests */
   if (!name) return m_process;
   int i;
   int no_processes =(int)pcs_vector.size();
   for(i=0;i<no_processes;i++)
   {
      m_process = pcs_vector[i];
      if ( strcmp(StrUp(m_process->pcs_primary_function_name[0]),StrUp(name))==0 )
      {
         break;
      }
   }
   return m_process;
}


//SB: new 3912
CRFProcess *CRFProcess::GetProcessByNumber(int number)
{
   CRFProcess *m_process = NULL;
   /* Tests */
   if (number < 1) return m_process;
   int i;
   int no_processes =(int)pcs_vector.size();
   for(i=0;i<no_processes;i++)
   {
      m_process = pcs_vector[i];
      if ( m_process->pcs_number == number )
      {
         break;
      }
   }
   return m_process;
}


//////////////////////////////////////////////////////////////////////////
// Configuration
//////////////////////////////////////////////////////////////////////////

/**************************************************************************
FEMLib-Method:
Task:
Programing:
02/2003 OK Implementation
08/2003 WW Modified to fit monolithic scheme
02/2005 OK Unsaturated flow (Richards model)
02/2005 MB string
05/2005 WW/DL Dymanic problem
01/2006 YD Dual Richards
OKToDo switch to char
03/2009 PCH PS_GLOBAL
get rid of type
**************************************************************************/
void CRFProcess::Config(void)
{
   std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
   // Set mesh pointer to corresponding mesh
   m_msh = FEMGet(pcs_type_name);
   if (!m_msh)
   {
      cout << "Error in CRFProcess::Config - no MSH data" << endl;
      return;
   }
   CheckMarkedElement();                          //WW

   if (continuum_vector.size() == 0)              // YD
      continuum_vector.push_back(1.0);

   //	if (_pcs_type_name.compare("LIQUID_FLOW") == 0) {
   if (this->getProcessType() == LIQUID_FLOW || this->getProcessType() == FLUID_FLOW)
   {
      std::cout << "CRFProcess::Config LIQUID_FLOW" << std::endl;
      type = 1;
      ConfigLiquidFlow();
   }
   //	if (_pcs_type_name.compare("GROUNDWATER_FLOW") == 0) {
   if (this->getProcessType() == GROUNDWATER_FLOW)
   {
      type = 1;
      ConfigGroundwaterFlow();
   }
   //	if (_pcs_type_name.compare("RICHARDS_FLOW") == 0) {
   if (this->getProcessType() == RICHARDS_FLOW)
   {
      if (continuum_vector.size() > 1)
      {
         RD_Process = true;
         type = 22;
      } else
      type = 14;
      ConfigUnsaturatedFlow();
   }
   //	if (_pcs_type_name.compare("OVERLAND_FLOW") == 0) {
   if (this->getProcessType() == OVERLAND_FLOW)
   {
      type = 66;
      max_dim = 1;
      ConfigGroundwaterFlow();
   }
   //	if (_pcs_type_name.compare("AIR_FLOW") == 0) { //OK
   if (this->getProcessType() == AIR_FLOW)        //OK
   {
      type = 5;
      ConfigGasFlow();
   }
   //	if (_pcs_type_name.compare("TWO_PHASE_FLOW") == 0) {
   if (this->getProcessType() == TWO_PHASE_FLOW)
   {
      type = 12;
      ConfigMultiphaseFlow();
   }
   //	if (_pcs_type_name.compare("COMPONENTAL_FLOW") == 0) {
   //	if (COMPONENTAL_FLOW) {
   //		type = 11;
   //		ConfigNonIsothermalFlow();
   //	}
   //	if (_pcs_type_name.compare("HEAT_TRANSPORT") == 0) {
   if (this->getProcessType() == HEAT_TRANSPORT)
   {
      type = 3;
      ConfigHeatTransport();
   }
   //	if (_pcs_type_name.compare("MASS_TRANSPORT") == 0) {
   if (this->getProcessType() == MASS_TRANSPORT)
   {
      type = 2;
      ConfigMassTransport();
   }
   //	if (_pcs_type_name.find("DEFORMATION") != string::npos)
   if (isDeformationProcess(this->getProcessType()))
      ConfigDeformation();
   //	if (_pcs_type_name.find("FLUID_MOMENTUM") != string::npos
   if (this->getProcessType() == FLUID_MOMENTUM)
   {
      type = 55;                                  //WW
      ConfigFluidMomentum();
   }
   //	if (_pcs_type_name.find("RANDOM_WALK") != string::npos) {
   if (this->getProcessType() == RANDOM_WALK)
   {
      type = 55;                                  //WW
      ConfigRandomWalk();
   }
   //	if (_pcs_type_name.find("MULTI_PHASE_FLOW") != string::npos) {//24.02.2007 WW
   if (this->getProcessType() == MULTI_PHASE_FLOW)//24.02.2007 WW
   {
      type = 1212;
      ConfigMultiPhaseFlow();
   }
   //	if (_pcs_type_name.find("PS_GLOBAL") != string::npos) {//24.02.2007 WW
   if (this->getProcessType() == PS_GLOBAL)       //24.02.2007 WW
   {
      type = 1313;
      ConfigPS_Global();
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2003 OK Implementation
        WW Splitted for processes
last modified:
02/2005 MB Pressure version for LIQUID Flow
**************************************************************************/
void CRFProcess::ConfigLiquidFlow()
{
   pcs_num_name[0] = "PRESSURE0";
   pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
   // NOD values
   pcs_number_of_primary_nvals = 1;
   pcs_number_of_secondary_nvals = 0;
   pcs_primary_function_name[0] = "PRESSURE1";
   pcs_primary_function_unit[0] = "Pa";
   // ELE values
   pcs_number_of_evals = 6;
   pcs_eval_name[0] = "VOLUME";
   pcs_eval_unit[0] = "m3";
   pcs_eval_name[1] = "VELOCITY1_X";
   pcs_eval_unit[1] = "m/s";
   pcs_eval_name[2] = "VELOCITY1_Y";
   pcs_eval_unit[2] = "m/s";
   pcs_eval_name[3] = "VELOCITY1_Z";
   pcs_eval_unit[3] = "m/s";
   pcs_eval_name[4] = "POROSITY";                 //MX, test for n=n(c), 04.2005
   pcs_eval_unit[4] = "-";
   pcs_eval_name[5] = "PERMEABILITY";             //JTARON 2010 -- need this for index call of heterogeneous permeability
   pcs_eval_unit[5] = "m2";
   //----------------------------------------------------------------------
   // Secondary variables
   pcs_number_of_secondary_nvals = 0;             //WW
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "HEAD";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;               //WW
   //WW
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;               //WW
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;               //WW
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;               //WW

   //WW / TF
   // Output material parameters
   configMaterialParameters();
}


/**************************************************************************
FEMLib-Method:
03/2003 OK Implementation
        WW Splitted for processes
02/2005 MB head version for GroundwaterFlow
08/2006 OK FLUX
**************************************************************************/
void CRFProcess::ConfigGroundwaterFlow()
{
   pcs_num_name[0] = "HEAD";
   pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_HEAD";
   // NOD values
   pcs_number_of_primary_nvals = 1;
   pcs_primary_function_name[0] = "HEAD";
   pcs_primary_function_unit[0] = "m";
   // ELE values
   pcs_number_of_evals = 6;
   pcs_eval_name[0] = "VOLUME";
   pcs_eval_unit[0] = "m3";
   pcs_eval_name[1] = "VELOCITY1_X";
   pcs_eval_unit[1] = "m/s";
   pcs_eval_name[2] = "VELOCITY1_Y";
   pcs_eval_unit[2] = "m/s";
   pcs_eval_name[3] = "VELOCITY1_Z";
   pcs_eval_unit[3] = "m/s";
   pcs_eval_name[4] = "PERMEABILITY";
   pcs_eval_unit[4] = "m^2";
   pcs_eval_name[5] = "POROSITY";
   pcs_eval_unit[5] = "-";
   //----------------------------------------------------------------------
   // Secondary variables
   pcs_number_of_secondary_nvals = 5;
   pcs_secondary_function_name[0] = "FLUX";
   pcs_secondary_function_unit[0] = "m3/s";
   pcs_secondary_function_timelevel[0] = 1;
   pcs_secondary_function_name[1] = "WDEPTH";
   pcs_secondary_function_unit[1] = "m";
   pcs_secondary_function_timelevel[1] = 1;
   pcs_secondary_function_name[2] = "COUPLING";   //JOD
   pcs_secondary_function_unit[2] = "m/s";
   pcs_secondary_function_timelevel[2] = 0;
   pcs_secondary_function_name[3] = "COUPLING";   //JOD
   pcs_secondary_function_unit[3] = "m/s";
   pcs_secondary_function_timelevel[3] = 1;
   pcs_secondary_function_name[4] = "STORE";      //JOD  for subtiming 4.7.10
   pcs_secondary_function_unit[4] = "m";
   pcs_secondary_function_timelevel[4] = 1;

   pcs_number_of_secondary_nvals = 5;             //WW
   //WW
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;               //WW
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;               //WW
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;               //WW
   //----------------------------------------------------------------------
   //WW / TF
   // Output material parameters
   configMaterialParameters();

   if (m_msh)
      m_msh->DefineMobileNodes(this);
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
10/2004 OK Implementation
last modified:
**************************************************************************/
void CRFProcess::ConfigGasFlow()
{
   //----------------------------------------------------------------------
   // Primary variables - NOD values
   pcs_number_of_primary_nvals = 1;
   pcs_number_of_secondary_nvals = 0;
   pcs_primary_function_name[0] = "PRESSURE1";
   pcs_primary_function_unit[0] = "Pa";
   //----------------------------------------------------------------------
   // Secondary variables - NOD values
   pcs_number_of_secondary_nvals = 1;
   pcs_secondary_function_name[0] = "NOD_MASS_FLUX";
   pcs_secondary_function_unit[0] = "kg/s";
   //----------------------------------------------------------------------
   // ELE values
   pcs_number_of_evals = 3;
   pcs_eval_name[0] = "VELOCITY1_X";
   pcs_eval_unit[0] = "m/s";
   pcs_eval_name[1] = "VELOCITY1_Y";
   pcs_eval_unit[1] = "m/s";
   pcs_eval_name[2] = "VELOCITY1_Z";
   pcs_eval_unit[2] = "m/s";
   //----------------------------------------------------------------------
   // NUM
   pcs_num_name[0] = "PRESSURE0";
   pcs_sol_name    = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2003 OK Implementation
        WW Splitted for processes
last modified:
**************************************************************************/
void CRFProcess::ConfigMultiphaseFlow()
{
   switch(pcs_type_number)
   {
      case 0:
         pcs_num_name[0]= "PRESSURE0";
         pcs_sol_name   = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
         break;
      case 1:
         pcs_num_name[0]= "SATURATION0";
         pcs_sol_name   = "LINEAR_SOLVER_PROPERTIES_SATURATION1";
         break;
   }
   //----------------------------------------------------------------------
   // NOD Primary variables
   pcs_number_of_primary_nvals=1;
   switch(pcs_type_number)
   {
      case 0:
         pcs_primary_function_name[0] = "PRESSURE1";
         pcs_primary_function_unit[0] = "Pa";
         break;
      case 1:
         pcs_primary_function_name[0] = "SATURATION2";
         pcs_primary_function_unit[0] = "m3/m3";
         break;
   }
   //----------------------------------------------------------------------
   // NOD Secondary variables
   pcs_number_of_secondary_nvals = 6;
   switch(pcs_type_number)
   {
      case 0:
         pcs_secondary_function_name[0]      = "PRESSURE2";
         pcs_secondary_function_unit[0]      = "Pa";
         pcs_secondary_function_timelevel[0] = 0;
         pcs_secondary_function_name[1]      = "PRESSURE2";
         pcs_secondary_function_unit[1]      = "Pa";
         pcs_secondary_function_timelevel[1] = 1;
         pcs_secondary_function_name[2]      = "PRESSURE_CAP";
         pcs_secondary_function_unit[2]      = "Pa";
         pcs_secondary_function_timelevel[2] = 0;
         pcs_secondary_function_name[3]      = "FLUX";
         pcs_secondary_function_unit[3]      = "m3/s";
         pcs_secondary_function_timelevel[3] = 0;
         pcs_secondary_function_name[4]      = "DENSITY1";
         pcs_secondary_function_unit[4]      = "kg/m3";
         pcs_secondary_function_timelevel[4] = 1;
         pcs_secondary_function_name[5]      = "VISCOSITY1";
         pcs_secondary_function_unit[5]      = "Pa s";
         pcs_secondary_function_timelevel[5] = 1;
         break;
      case 1:
         pcs_secondary_function_name[0]      = "SATURATION1";
         pcs_secondary_function_timelevel[0] = 0;
         pcs_secondary_function_unit[0]      = "m3/m3";
         pcs_secondary_function_name[1]      = "SATURATION1";
         pcs_secondary_function_timelevel[1] = 1;
         pcs_secondary_function_unit[1]      = "m3/m3";
         pcs_secondary_function_name[2]      = "PRESSURE_CAP";
         pcs_secondary_function_unit[2]      = "Pa";
         pcs_secondary_function_timelevel[2] = 1;
         pcs_secondary_function_name[3]      = "FLUX";
         pcs_secondary_function_unit[3]      = "m3/s";
         pcs_secondary_function_timelevel[3] = 1;
         pcs_secondary_function_name[4]      = "DENSITY2";
         pcs_secondary_function_unit[4]      = "kg/m3";
         pcs_secondary_function_timelevel[4] = 1;
         pcs_secondary_function_name[5]      = "VISCOSITY2";
         pcs_secondary_function_unit[5]      = "Pa s";
         pcs_secondary_function_timelevel[5] = 1;
         break;
   }
   //----------------------------------------------------------------------
   // ELE values
   pcs_number_of_evals = 7;
   switch(pcs_type_number)
   {
      case 0:
         pcs_eval_name[0] = "VELOCITY1_X";
         pcs_eval_unit[0] = "m/s";
         pcs_eval_name[1] = "VELOCITY1_Y";
         pcs_eval_unit[1] = "m/s";
         pcs_eval_name[2] = "VELOCITY1_Z";
         pcs_eval_unit[2] = "m/s";
         pcs_eval_name[3] = "POROSITY1";          //MX 03.2005
         pcs_eval_unit[3] = "-";
         pcs_eval_name[4] = "POROSITY1_IL";       //MX 03.2005
         pcs_eval_unit[4] = "-";
         pcs_eval_name[5] = "PERMEABILITY1";      //MX 03.2005
         pcs_eval_unit[5] = "-";
         pcs_eval_name[6] = "POROSITY1_SW";       //MX 03.2005
         pcs_eval_unit[6] = "-";
         break;
      case 1:
         pcs_eval_name[0] = "VELOCITY2_X";
         pcs_eval_unit[0] = "m/s";
         pcs_eval_name[1] = "VELOCITY2_Y";
         pcs_eval_unit[1] = "m/s";
         pcs_eval_name[2] = "VELOCITY2_Z";
         pcs_eval_unit[2] = "m/s";
         pcs_eval_name[3] = "POROSITY";           //MX 03.2005
         pcs_eval_unit[3] = "-";
         pcs_eval_name[4] = "POROSITY_IL";        //MX 03.2005
         pcs_eval_unit[4] = "-";
         pcs_eval_name[5] = "PERMEABILITY";       //MX 03.2005
         pcs_eval_unit[5] = "-";
         pcs_eval_name[6] = "POROSITY_SW";        //MX 03.2005
         pcs_eval_unit[6] = "-";
         break;
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2003 OK Implementation
        WW Splitted for processes
last modified:
**************************************************************************/
void CRFProcess::ConfigNonIsothermalFlow()
{
   //----------------------------------------------------------------------
   // Primary variables
   pcs_number_of_primary_nvals=1;
   switch(pcs_type_number)
   {
      case 0:
         pcs_num_name[0]= "PRESSURE0";
         pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
         pcs_primary_function_name[0] = "PRESSURE1";
         pcs_primary_function_unit[0] = "Pa";
         break;
      case 1:
         pcs_num_name[0]= "SATURATION0";
         pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_SATURATION1";
         pcs_primary_function_name[0] = "SATURATION2";
         pcs_primary_function_unit[0] = "m3/m3";
         break;
   }
   //----------------------------------------------------------------------
   // Secondary variables
   pcs_number_of_secondary_nvals = 6;
   switch(pcs_type_number)
   {
      case 0:
         pcs_secondary_function_name[0] = "PRESSURE2";
         pcs_secondary_function_timelevel[0] = 0;
         pcs_secondary_function_unit[0] = "Pa";
         pcs_secondary_function_name[1] = "PRESSURE2";
         pcs_secondary_function_timelevel[1] = 1;
         pcs_secondary_function_unit[1] = "Pa";
         pcs_secondary_function_name[2] = "MASS_FRACTION1";
         pcs_secondary_function_timelevel[2] = 0;
         pcs_secondary_function_unit[2] = "kg/kg";
         pcs_secondary_function_name[3] = "MASS_FRACTION1";
         pcs_secondary_function_timelevel[3] = 1;
         pcs_secondary_function_unit[3] = "kg/kg";
         pcs_secondary_function_name[4] = "PRESSURE_CAP";
         pcs_secondary_function_timelevel[4] = 0;
         pcs_secondary_function_unit[4] = "Pa";
         pcs_secondary_function_name[5] = "DENSITY1";
         pcs_secondary_function_timelevel[5] = 1;
         pcs_secondary_function_unit[5] = "kg/m3";
         break;
      case 1:
         pcs_secondary_function_name[0] = "SATURATION1";
         pcs_secondary_function_timelevel[0] = 0;
         pcs_secondary_function_unit[0] = "m3/m3";
         pcs_secondary_function_name[1] = "SATURATION1";
         pcs_secondary_function_timelevel[1] = 1;
         pcs_secondary_function_unit[1] = "m3/m3";
         pcs_secondary_function_name[2] = "MASS_FRACTION2";
         pcs_secondary_function_timelevel[2] = 0;
         pcs_secondary_function_unit[2] = "kg/kg";
         pcs_secondary_function_name[3] = "MASS_FRACTION2";
         pcs_secondary_function_timelevel[3] = 1;
         pcs_secondary_function_unit[3] = "kg/kg";
         pcs_secondary_function_name[4] = "PRESSURE_CAP";
         pcs_secondary_function_timelevel[4] = 1;
         pcs_secondary_function_unit[4] = "Pa";
         pcs_secondary_function_name[5] = "DENSITY2";
         pcs_secondary_function_timelevel[5] = 1;
         pcs_secondary_function_unit[5] = "kg/m3";
         break;
   }
   // Node
   pcs_number_of_primary_nvals=1;
   // ELE values
   pcs_number_of_evals = 14;
   pcs_eval_name[0]  = "COMP_FLUX";
   pcs_eval_name[1]  = "POROSITY";
   pcs_eval_name[2]  = "PERMEABILITY";
   pcs_eval_name[3]  = "VELOCITY1_X";
   pcs_eval_name[4]  = "VELOCITY1_Y";
   pcs_eval_name[5]  = "VELOCITY1_Z";
   pcs_eval_name[6]  = "VELOCITY2_X";
   pcs_eval_name[7]  = "VELOCITY2_Y";
   pcs_eval_name[8]  = "VELOCITY2_Z";
   pcs_eval_name[9]  = "POROSITY_IL";
   pcs_eval_name[10] = "VoidRatio";
   pcs_eval_name[11] = "PorosityChange";
   pcs_eval_name[12] = "n_sw_Rate";
   pcs_eval_name[13]  = "POROSITY_SW";
   pcs_eval_unit[0]  = "kg/s";
   pcs_eval_unit[1]  = "m3/m3";
   pcs_eval_unit[2]  = "m2";
   pcs_eval_unit[3]  = "m/s";
   pcs_eval_unit[4]  = "m/s";
   pcs_eval_unit[5]  = "m/s";
   pcs_eval_unit[6]  = "m/s";
   pcs_eval_unit[7]  = "m/s";
   pcs_eval_unit[8]  = "m/s";
   pcs_eval_unit[9]  = "-";
   pcs_eval_unit[10] = "-";
   pcs_eval_unit[11] = "-";
   pcs_eval_unit[12] = "-";
   pcs_eval_unit[13] = "-";
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2003 OK Implementation
        WW Splitted for processes
last modified:
**************************************************************************/
//void CRFProcess::ConfigNonIsothermalFlowRichards()

/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2004 SB Implementation
        WW Splitted for processes
01/2006 OK Tests
08/2006 OK FLUX
**************************************************************************/
void CRFProcess::ConfigMassTransport()
{
   long comp=1;
   /* count transport processes */
   pcs_component_number++;
   comp = pcs_component_number;
   // 1 NOD values
   // 1.1 primary variables
   pcs_number_of_primary_nvals = 1;
   pcs_primary_function_name[0]= new char[80];
   //  sprintf(pcs_primary_function_name[0], "%s%li","CONCENTRATION",comp);
   //----------------------------------------------------------------------
   // Tests
   int size;
   size = (int)cp_vec.size();
   int comb;                                      //OK411
   comb = pcs_component_number;

   if((int)cp_vec.size()<pcs_component_number+1)
   {
      cout << "Error in CRFProcess::ConfigMassTransport - not enough MCP data" << endl;
      return;
   }
   //----------------------------------------------------------------------
   pcs_primary_function_name[0] = cp_vec[pcs_component_number]->compname.c_str();
   //sprintf(pcs_primary_function_name[0], "%s", cp_vec[pcs_component_number]->compname.c_str());
   pcs_primary_function_unit[0] = "kg/m3";        //SB
   /* SB: Eintrag component name in Ausgabestruktur */ //SB:todo : just one phase todo:name
   /*
         pcs_primary_function_name[0] = GetTracerCompName(0,this->pcs_component_number-1);
         name_initial_condition_tracer_component = pcs_primary_function_name[0];
         pcs_ic_name_mass = pcs_primary_function_name[0];
   */
   // 1.2 secondary variables
   pcs_number_of_secondary_nvals = 2;             //SB3909
   pcs_secondary_function_name[0]= new char[80];
   char pcs_secondary_function_name_tmp [80];
   sprintf(pcs_secondary_function_name_tmp, "%s%li","MASS_FLUX_",comp);
   pcs_secondary_function_name[0] = pcs_secondary_function_name_tmp;
   //      pcs_secondary_function_name[0] = "MASS_FLUX1";
   pcs_secondary_function_unit[0] = "kg/m3/s";
   pcs_secondary_function_timelevel[0] = 0;
   pcs_secondary_function_name[1]= new char[80];
   sprintf(pcs_secondary_function_name_tmp, "%s%li","MASS_FLUX_",comp);
   pcs_secondary_function_name[1] = pcs_secondary_function_name_tmp;
   pcs_secondary_function_unit[1] = "kg/m3/s";
   pcs_secondary_function_timelevel[1] = 1;
   //KG44 added secondary function for adaptive time stepping
   if (adaption)
   {
      pcs_number_of_secondary_nvals = 3;
      pcs_secondary_function_name[2]= new char[80];
      sprintf(pcs_secondary_function_name_tmp, "%s%li","CONC_BACK_",comp);
      pcs_secondary_function_name[2] = pcs_secondary_function_name_tmp;
      pcs_secondary_function_unit[2] = "kg/m3";
      pcs_secondary_function_timelevel[2] = 0;
   }
   //OK  LOPCalcSecondaryVariables_USER = MTM2CalcSecondaryVariables;  //SB:todo
   // 2 ELE values
   pcs_number_of_evals = 0;
   //	  pcs_eval_name[0] = "Darcy velocity";
#ifdef REACTION_ELEMENT
   pcs_number_of_evals = 1;
   pcs_eval_name[0] = pcs_primary_function_name[0];
   pcs_eval_unit[0] = "mol/kgH2O";
#endif
   // 3 ELE matrices
   // NUM
   pcs_num_name[0] = "CONCENTRATION0";
   /* SB: immer solver properties der ersten Komponente nehmen */
                                                  //SB ??
   pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_CONCENTRATION1";
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2003 OK Implementation
        WW Splitted for processes
last modified:
**************************************************************************/
void CRFProcess::ConfigHeatTransport()
{
   pcs_num_name[0] = "TEMPERATURE0";
   pcs_sol_name    = "LINEAR_SOLVER_PROPERTIES_TEMPERATURE1";
   // NOD
   if((int)continuum_vector.size() == 1)
   {
      pcs_number_of_primary_nvals = 1;
      pcs_primary_function_name[0] = "TEMPERATURE1";
      pcs_primary_function_unit[0] = "K";
      pcs_number_of_secondary_nvals = 0;
#ifdef REACTION_ELEMENT
      pcs_number_of_evals = 1;                    //MX
      pcs_eval_name[0] = "TEMPERATURE1";
      pcs_eval_unit[0] = "K";
#endif
   }
   if((int)continuum_vector.size() == 2)
   {
      pcs_number_of_primary_nvals = 2;
      pcs_primary_function_name[0] = "TEMPERATURE1";
      pcs_primary_function_unit[0] = "K";
      pcs_number_of_primary_nvals = 2;
      pcs_primary_function_name[1] = "TEMPERATURE2";
      pcs_primary_function_unit[1] = "K";
      pcs_number_of_secondary_nvals = 0;
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2003 WW Implementation
last modified:
**************************************************************************/
void CRFProcess::ConfigDeformation()
{
   // Generate high order nodes for all elements.
   m_msh->GenerateHighOrderNodes();               //WW
   type = 4;
   //	if (_pcs_type_name.find("DEFORMATION") != string::npos
   //			&& _pcs_type_name.find("FLOW") != string::npos) {
   if (getProcessType() == DEFORMATION_FLOW)
   {
      type = 41;
      cpl_type_name = "MONOLITHIC";
      pcs_deformation = 11;
   }
   pcs_coupling_iterations = 10;

   CNumerics* num = NULL;
   for (size_t i = 0; i < num_vector.size(); i++)
   {
      num = num_vector[i];
      if (num->pcs_type_name.find("DEFORMATION") != string::npos)
      {
         num->pcs_type_name = convertProcessTypeToString(this->getProcessType());
         if (num->nls_method == 1)                // Newton-Raphson
         {
            pcs_deformation = 101;
            if (type == 41)
               pcs_deformation = 110;
         }
         break;
      }
   }

   // Prepare for restart
   //RFConfigRenumber();

   // Geometry dimension
   problem_dimension_dm = m_msh->GetCoordinateFlag() / 10;
   problem_2d_plane_dm = 1;

   // NUM
   pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_DISPLACEMENT1";
   pcs_num_name[0] = "DISPLACEMENT0";
   pcs_num_name[1] = "PRESSURE0";
   if (pcs_type_name_vector[0].find("DYNAMIC") != string::npos)
      VariableDynamics();
   else
      VariableStaticProblem();

   // Coupling
   for (size_t i = 0; i < GetPrimaryVNumber(); i++)
      Shift[i] = i * m_msh->GetNodesNumber(true);
   // OBJ names are set to PCS name
   if (type == 41)
      SetOBJNames();                              //OK->WW please put to Config()
}


/**************************************************************************
FEMLib-Method: Static problems
Task:
Programing:
05/2005 WW Implementation
last modified:
**************************************************************************/
void CRFProcess::VariableStaticProblem()
{
   //----------------------------------------------------------------------
   // NOD Primary functions
   pcs_number_of_primary_nvals = 2;               //OK distinguish 2/3D problems, problem_dimension_dm;
   dm_number_of_primary_nvals = 2;
   pcs_primary_function_name[0] = "DISPLACEMENT_X1";
   pcs_primary_function_name[1] = "DISPLACEMENT_Y1";
   pcs_primary_function_unit[0] = "m";
   pcs_primary_function_unit[1] = "m";
   if(max_dim==2)
   {
      pcs_number_of_primary_nvals = 3;
      dm_number_of_primary_nvals = 3;
      pcs_primary_function_name[2] = "DISPLACEMENT_Z1";
      pcs_primary_function_unit[2] = "m";
   }
   if(type==41)                                   //Monolithic scheme
   {
      pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE1";
      pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa";
      pcs_number_of_primary_nvals++;
   }
   //----------------------------------------------------------------------
   // NOD Secondary functions
   pcs_number_of_secondary_nvals = 0;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XX";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XY";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YY";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_ZZ";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XX";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XY";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YY";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_ZZ";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_PLS";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   //  pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "POROPRESSURE0";
   //  pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   //  pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   //  pcs_number_of_secondary_nvals++;

   if(max_dim==2)                                 // 3D
   {
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XZ";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YZ";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XZ";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "--";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YZ";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "--";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
   }

}


/**************************************************************************
FEMLib-Method: Dynamic problems
Task:
Programing:
05/2005 WW/LD Implementation
last modified:
**************************************************************************/
void CRFProcess::VariableDynamics()
{
   //----------------------------------------------------------------------
   // NOD Primary functions
   pcs_number_of_primary_nvals = 2;
   dm_number_of_primary_nvals = 2;
   pcs_primary_function_name[0] = "ACCELERATION_X1";
   pcs_primary_function_name[1] = "ACCELERATION_Y1";
   pcs_primary_function_unit[0] = "m/s^2";
   pcs_primary_function_unit[1] = "m/s^2";
   if(max_dim==2)
   {
      pcs_number_of_primary_nvals = 3;
      dm_number_of_primary_nvals = 3;
      pcs_primary_function_name[2] = "ACCELERATION_Z1";
      pcs_primary_function_unit[2] = "m/s^2";
   }
   pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE_RATE1";
   pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa/s";
   pcs_number_of_primary_nvals++;

   //----------------------------------------------------------------------
   // NOD Secondary functions
   pcs_number_of_secondary_nvals = 0;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XX";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XY";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YY";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_ZZ";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   //  pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "POROPRESSURE0";
   //  pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   //  pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   //  pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XX";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XY";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YY";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_ZZ";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_PLS";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DISPLACEMENT_X1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DISPLACEMENT_Y1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_DM_X";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_DM_Y";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   // 3D
   if(max_dim==2)                                 // 3D
   {
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XZ";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YZ";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XZ";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YZ";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DISPLACEMENT_Z1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_DM_Z";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
02/2005 OK Implementation
02/2006 OK FLUX
**************************************************************************/
void CRFProcess::ConfigUnsaturatedFlow()
{
   if ((int) continuum_vector.size() == 1)
   {
      // 1.1 primary variables
      pcs_number_of_primary_nvals = 1;
      pcs_primary_function_name[0] = "PRESSURE1";
      pcs_primary_function_unit[0] = "Pa";
      // 1.2 secondary variables
      //OK LOPCalcSecondaryVariables_USER = MMPCalcSecondaryVariablesRichards; // p_c and S^l
      pcs_number_of_secondary_nvals = 0;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "SATURATION1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "SATURATION1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "PRESSURE_CAP1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
      pcs_number_of_secondary_nvals++;
                                                  //MB
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "FLUX";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
      pcs_number_of_secondary_nvals++;
                                                  //MB
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "FLUX";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      //TEST
      //#define DECOVALEX
#ifdef DECOVALEX
      // DECOVALEX Test
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_I";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
      pcs_number_of_secondary_nvals++;
#endif
      /* if(adaption) //WW, JOD removed
       {
       pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STORAGE_P";
       pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
       pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
       pcs_number_of_secondary_nvals++;
       }*/
      // Nodal velocity. WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "VELOCITY_X1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "VELOCITY_Y1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "VELOCITY_Z1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
                                                  //JOD
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "COUPLING";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
      pcs_number_of_secondary_nvals++;
                                                  //JOD
      pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "COUPLING";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;
   }
   else if ((int) continuum_vector.size() == 2)
   {
      dof = 2;                                    //WW
      // 1.1 primary variables
      pcs_number_of_primary_nvals = 2;            //YD
      pcs_primary_function_name[0] = "PRESSURE1";
      pcs_primary_function_unit[0] = "Pa";
      pcs_primary_function_name[1] = "PRESSURE2";
      pcs_primary_function_unit[1] = "Pa";
      // 1.2 secondary variables
      //OK LOPCalcSecondaryVariables_USER = MMPCalcSecondaryVariablesRichards; // p_c and S^l
      pcs_number_of_secondary_nvals = 0;
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "SATURATION1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "SATURATION1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "SATURATION2";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "SATURATION2";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "PRESSURE_CAP1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "PRESSURE_CAP2";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
      pcs_number_of_secondary_nvals++;            //WW
      // Nodal velocity. WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "VELOCITY_X1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "VELOCITY_Y1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "VELOCITY_Z1";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "VELOCITY_X2";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "VELOCITY_Y2";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;            //WW
      pcs_secondary_function_name[pcs_number_of_secondary_nvals]
         = "VELOCITY_Z2";
      pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
      pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
      pcs_number_of_secondary_nvals++;            //WW
                                                  // 03.03.2008. WW
      for (size_t i = 0; i < GetPrimaryVNumber(); i++)
         Shift[i] = i * m_msh->GetNodesNumber(true);
   }

   // Output material parameters
   //WW // TF
   configMaterialParameters();

   // 2 ELE values
   pcs_number_of_evals = 8;
   pcs_eval_name[0] = "VELOCITY1_X";
   pcs_eval_unit[0] = "m/s";
   pcs_eval_name[1] = "VELOCITY1_Y";
   pcs_eval_unit[1] = "m/s";
   pcs_eval_name[2] = "VELOCITY1_Z";
   pcs_eval_unit[2] = "m/s";
   pcs_eval_name[3] = "POROSITY";                 //MX 11.2005
   pcs_eval_unit[3] = "-";
   pcs_eval_name[4] = "POROSITY_IL";              //MX 11.2005
   pcs_eval_unit[4] = "-";
   pcs_eval_name[5] = "PERMEABILITY";             //MX 11.2005
   pcs_eval_unit[5] = "-";
   pcs_eval_name[6] = "n_sw";                     //MX 11.2005
   pcs_eval_unit[6] = "-";
   pcs_eval_name[7] = "n_sw_rate";                //MX 11.2005
   pcs_eval_unit[7] = "-";
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 05/2005 PCH Implementation
 last modified:
 **************************************************************************/
void CRFProcess::ConfigFluidMomentum()
{
   //pcs_num_name[0] = "VELOCITY1_X";
   // Nothing added in terms of matrix solver.
   // Just linear solver is good enough.
   pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
   // NOD values
   pcs_number_of_primary_nvals = 3;
   pcs_primary_function_name[0] = "VELOCITY1_X";
   pcs_primary_function_unit[0] = "m/s";
   pcs_primary_function_name[1] = "VELOCITY1_Y";
   pcs_primary_function_unit[1] = "m/s";
   pcs_primary_function_name[2] = "VELOCITY1_Z";
   pcs_primary_function_unit[2] = "m/s";

   // I'm adding this to initialize for Fluid Momentum process
   pcs_number_of_secondary_nvals = 0;
   pcs_number_of_evals = 3;

   pcs_eval_name[0] = "VELOCITY1_X";
   pcs_eval_unit[0] = "m/s";
   pcs_eval_name[1] = "VELOCITY1_Y";
   pcs_eval_unit[1] = "m/s";
   pcs_eval_name[2] = "VELOCITY1_Z";
   pcs_eval_unit[2] = "m/s";
}


/**************************************************************************/
void CRFProcess::ConfigRandomWalk()
{
   // Nothing added in terms of matrix solver.
   // Just linear solver is good enough.
   pcs_sol_name    = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";

   // NOD values
   pcs_number_of_primary_nvals = 0;
   pcs_number_of_secondary_nvals = 0;

   // 2 ELE values
   pcs_number_of_evals = 1;
   pcs_eval_name[0] = "CONCENTRATION0";
   pcs_eval_unit[0]  = "kg/m3";

}


/**************************************************************************
FEMLib-Method: For non-isothermal multi-phase flow
Task:
Programing:
02/2007 WW Implementation
**************************************************************************/
void CRFProcess::ConfigMultiPhaseFlow()
{
   dof = 2;
   // 1.1 primary variables
   pcs_number_of_primary_nvals = 2;
   pcs_primary_function_name[0] = "PRESSURE1";
   pcs_primary_function_unit[0] = "Pa";
   pcs_primary_function_name[1] = "PRESSURE2";
   pcs_primary_function_unit[1] = "Pa";
   // 1.2 secondary variables
   pcs_number_of_secondary_nvals = 0;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   // Nodal velocity.
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X2";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y2";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z2";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_W";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   //
   for(size_t i=0; i<GetPrimaryVNumber(); i++)    // 03.03.2008. WW
      Shift[i] = i*m_msh->GetNodesNumber(true);
}


/**************************************************************************
FEMLib-Method: For PS model for multiphase flow
Task:
Programing:
03/2009 PCH Implementation
**************************************************************************/
void CRFProcess::ConfigPS_Global()
{
   dof = 2;
   // 1.1 primary variables
   pcs_number_of_primary_nvals = 2;
   pcs_primary_function_name[0] = "PRESSURE1";
   pcs_primary_function_unit[0] = "Pa";
   pcs_primary_function_name[1] = "SATURATION2";
   pcs_primary_function_unit[1] = "m3/m3";
   // 1.2 secondary variables
   pcs_number_of_secondary_nvals = 0;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE2";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_CAP";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   // Nodal velocity.
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X2";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y2";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;
   pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z2";
   pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
   pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
   pcs_number_of_secondary_nvals++;

   //
   for(size_t i=0; i<GetPrimaryVNumber(); i++)    // 03.03.2008. WW
      Shift[i] = i*m_msh->GetNodesNumber(true);
}


//////////////////////////////////////////////////////////////////////////
// Configuration NOD
//////////////////////////////////////////////////////////////////////////
#ifndef NEW_EQS                                   //WW. 07.11.2008
/*************************************************************************
ROCKFLOW - Function:
Task: Config node values
Programming: 02/2003 OK Implementation
04/2004   WW   Modification for 3D problems
last modified:
**************************************************************************/
void CRFProcess::ConfigNODValues1(void)
{
   int i;
   int pcs_nval=0;
   const int DOF = GetPrimaryVNumber();
   anz_nval0 = anz_nval;
   number_of_nvals = 2*DOF + pcs_number_of_secondary_nvals;
   //NVAL
   pcs_nval_data = (PCS_NVAL_DATA *) \
      Malloc(number_of_nvals*sizeof(PCS_NVAL_DATA));
   /*----------------------------------------------------------------*/
   for(i=0;i<DOF;i++)
   {
      /* Primary variable - old time */
      //NVAL pcs_nval_data[pcs_nval] = (PCS_NVAL_DATA *) Malloc(sizeof(PCS_NVAL_DATA));
      strcpy(pcs_nval_data[pcs_nval].name,pcs_primary_function_name[pcs_nval-i]);
      //Change name for the previous time level
      /*
      char *ch = strchr(pcs_nval_data[pcs_nval].name, '1');
      if( ch != NULL )
      {
        int pos = ch-pcs_nval_data[pcs_nval].name;
        pcs_nval_data[pcs_nval].name[pos]='0';
      }
      else
        strcat(pcs_nval_data[pcs_nval].name, "0");   */
      //-------------------------------------------------------------------------------
      strcpy(pcs_nval_data[pcs_nval].einheit,pcs_primary_function_unit[pcs_nval-i]);
      pcs_nval_data[pcs_nval].timelevel = 0;
      pcs_nval_data[pcs_nval].speichern = 0;      //WW
      pcs_nval_data[pcs_nval].laden = 0;
      pcs_nval_data[pcs_nval].restart = 1;
      pcs_nval_data[pcs_nval].adapt_interpol = 1;
      pcs_nval_data[pcs_nval].vorgabe = 0.0;
#ifdef PCS_NOD
      pcs_nval_data[pcs_nval].nval_index = pcs_nval;
#else
      pcs_nval_data[pcs_nval].nval_index = anz_nval + pcs_nval;
#endif
      pcs_nval++;
      /* Primary variable - new time */
      //NVAL pcs_nval_data[pcs_nval] = (PCS_NVAL_DATA *) Malloc(sizeof(PCS_NVAL_DATA));
      strcpy(pcs_nval_data[pcs_nval].name,pcs_primary_function_name[pcs_nval-i-1]);
      strcpy(pcs_nval_data[pcs_nval].einheit,pcs_primary_function_unit[pcs_nval-i-1]);
      pcs_nval_data[pcs_nval].timelevel = 1;
      pcs_nval_data[pcs_nval].speichern = 1;

      pcs_nval_data[pcs_nval].laden = 0;
      pcs_nval_data[pcs_nval].restart = 1;
      pcs_nval_data[pcs_nval].adapt_interpol = 1;
      pcs_nval_data[pcs_nval].vorgabe = 0.0;
#ifdef PCS_NOD
      pcs_nval_data[pcs_nval].nval_index = pcs_nval;
#else
      pcs_nval_data[pcs_nval].nval_index = anz_nval + pcs_nval;
#endif
      pcs_nval++;
   }

   /*----------------------------------------------------------------*/
   /* Secondary variables */
   for(i=0;i<pcs_number_of_secondary_nvals;i++)
   {
      //NVAL pcs_nval_data[pcs_nval] = (PCS_NVAL_DATA *) Malloc(sizeof(PCS_NVAL_DATA));
      strcpy(pcs_nval_data[pcs_nval].name,pcs_secondary_function_name[i]);
      strcpy(pcs_nval_data[pcs_nval].einheit,pcs_secondary_function_unit[i]);
      //  pcs_nval_data[i+2]->timelevel = 1; // always at new time level
      pcs_nval_data[pcs_nval].timelevel = pcs_secondary_function_timelevel[i];
      if(pcs_nval_data[pcs_nval].timelevel==1)
         pcs_nval_data[pcs_nval].speichern = 1;
      else
         pcs_nval_data[pcs_nval].speichern = 0;
      pcs_nval_data[pcs_nval].laden = 0;
      pcs_nval_data[pcs_nval].restart = 1;
      pcs_nval_data[pcs_nval].adapt_interpol = 1;
      pcs_nval_data[pcs_nval].vorgabe = 0.0;
#ifdef PCS_NOD
      pcs_nval_data[pcs_nval].nval_index = pcs_nval;
#else
      pcs_nval_data[pcs_nval].nval_index = anz_nval + pcs_nval;
#endif
      if(type==4||type==41)
      {
         if(dm_number_of_primary_nvals==2||(dm_number_of_primary_nvals==3&&this->type==41))
         {
            //Block:
            //STRESS_ZX, STRESS_YZ, STRAIN_ZX, STRAIN_ZY and LUMPED_STRESS
            if(i==3||i==4||i==9||i==10||i==13)
               pcs_nval_data[pcs_nval].speichern = 0;
            if(!problem_2d_plane_dm)
            {
               //  Block STRESS_ZZ and STRAIN_ZZ
               if(i==5||i==11) pcs_nval_data[pcs_nval].speichern = 0;
            }
            //if(!pcs_plasticity)
            // {
            //   if(i==12) pcs_nval_data[pcs_nval].speichern = 0;  //STRAIN_PLS
            // }
         }
      }
      pcs_nval++;
   }
   pcs_nval_data = pcs_nval_data;
}


/*************************************************************************
ROCKFLOW - Function: CRFProcess::PCSConfigNODValues
Task: Config node values
Programming: 02/2003 OK Implementation
last modified:
**************************************************************************/
void CRFProcess::ConfigNODValues2(void)
{
   int i;

   number_of_nvals = 2*GetPrimaryVNumber() + pcs_number_of_secondary_nvals;
   for(i=0;i<number_of_nvals;i++)
   {
      ModelsAddNodeValInfoStructure(pcs_nval_data[i].name,\
         pcs_nval_data[i].einheit,\
         pcs_nval_data[i].speichern,\
         pcs_nval_data[i].laden,\
         pcs_nval_data[i].restart,\
         pcs_nval_data[i].adapt_interpol,\
         pcs_nval_data[i].vorgabe);
   }
}
#endif                                            //#ifndef NEW_EQS //WW. 07.11.2008
/**************************************************************************
FEMLib-Method:
Task:
Programing:
07/2004 OK Implementation
**************************************************************************/
int PCSGetNODValueIndex(const string &name,int timelevel)
{
   // PCS primary variables
   int pcs_vector_size =(int)pcs_vector.size();
   int i,j;
   CRFProcess *m_pcs = NULL;
   if(pcs_vector_size>0)
   {
      for(i=0;i<pcs_vector_size;i++)
      {
         m_pcs = pcs_vector[i];
         for(j=0;j<m_pcs->number_of_nvals;j++)
         {
            if((name.compare(m_pcs->pcs_nval_data[j].name)==0) && \
               (m_pcs->pcs_nval_data[j].timelevel==timelevel))
               return m_pcs->pcs_nval_data[j].nval_index;
         }
      }
   }
   cout << "Error in PCSGetNODValueIndex: " << name << endl;
   return -1;
}


//////////////////////////////////////////////////////////////////////////
// Configuration ELE
//////////////////////////////////////////////////////////////////////////

/*************************************************************************
ROCKFLOW - Function:
Task: Config element values
Programming: 02/2003 OK Implementation
last modified:
  06/2004  WW
**************************************************************************/
void CRFProcess::ConfigELEValues1(void)
{
   int i;
   if(pcs_number_of_evals)
      pcs_eval_data = (PCS_EVAL_DATA *) Malloc(pcs_number_of_evals*sizeof(PCS_EVAL_DATA));
   for(i=0;i<pcs_number_of_evals;i++)
   {
      //pcs_eval_data[i] = (PCS_EVAL_DATA *) Malloc(sizeof(PCS_EVAL_DATA));
      strcpy(pcs_eval_data[i].name,pcs_eval_name[i]);
      strcpy(pcs_eval_data[i].einheit,pcs_eval_unit[i]);
      pcs_eval_data[i].speichern = 1;
      pcs_eval_data[i].laden = 0;
      pcs_eval_data[i].restart = 1;
      pcs_eval_data[i].adapt_interpol = 1;
      pcs_eval_data[i].vorgabe = 0.0;
      pcs_eval_data[i].index = anz_eval + i;
      pcs_eval_data[i].eval_index = anz_eval + i; //SB
   }
}


/*************************************************************************
ROCKFLOW - Function:
Task: Config element values
Programming: 04/2003 OK Implementation
last modified:
**************************************************************************/
void CRFProcess::ConfigELEValues2(void)
{
   int i;
   for(i=0;i<pcs_number_of_evals;i++)
   {
      ModelsAddElementValInfoStructure(pcs_eval_data[i].name,\
         pcs_eval_data[i].einheit,\
         pcs_eval_data[i].speichern,\
         pcs_eval_data[i].laden,\
         pcs_eval_data[i].restart,\
         pcs_eval_data[i].adapt_interpol,\
         pcs_eval_data[i].vorgabe);
   }
}


/*************************************************************************
ROCKFLOW - Function: CRFProcess::PCSGetELEValueIndex
Task: Provide index for element values
Programming: 08/2003 SB Implementation
last modified:
**************************************************************************/
int PCSGetELEValueIndex(char *name)
{
   int i;
   CRFProcess *m_process = NULL;
   int j;
   int no_processes =(int)pcs_vector.size();
   for(j=0;j<no_processes;j++)
   {
      m_process = pcs_vector[j];
      for(i=0;i<m_process->pcs_number_of_evals;i++)
      {
         if(strcmp(m_process->pcs_eval_data[i].name,name)==0)
            return m_process->pcs_eval_data[i].eval_index;
      }
   }
   printf("PCSGetELEValueIndex Alert\n");
   printf("%s \n",name);
   return -1;
}


//////////////////////////////////////////////////////////////////////////
// Configuration ELE matrices
//////////////////////////////////////////////////////////////////////////

/**************************************************************************
FEMLib-Method:
Task:  Activate or deactivate elements specified in .pcs file
Programing:
05/2005 WW Implementation
**************************************************************************/
void CRFProcess::CheckMarkedElement()
{
   int i, j;
   long l;
   bool done;
   CElem *elem = NULL;
   CNode *node = NULL;
   for (l = 0; l < (long)m_msh->ele_vector.size(); l++)
   {
      elem = m_msh->ele_vector[l];
      done = false;
      for(i=0; i<NumDeactivated_SubDomains; i++)
      {
         if(elem->GetPatchIndex()== Deactivated_SubDomain[i])
         {
            elem->MarkingAll(false);
            done = true;
            break;
         }
      }
      if(done)
         continue;
      else
         elem->MarkingAll(true);

   }
   for (l = 0; l < (long)m_msh->nod_vector.size(); l++)
   {
      while(m_msh->nod_vector[l]->connected_elements.size())
         m_msh->nod_vector[l]->connected_elements.pop_back();
   }
   for (l = 0; l < (long)m_msh->ele_vector.size(); l++)
   {
      elem = m_msh->ele_vector[l];
      if(!elem->GetMark()) continue;
      for(i=0; i<elem->GetNodesNumber(m_msh->getOrder()); i++)
      {
         done = false;
         node = elem->GetNode(i);
         for(j=0; j<(int)node->connected_elements.size(); j++)
         {
            if(l==node->connected_elements[j])
            {
               done = true;
               break;
            }
         }
         if(!done)
            node->connected_elements.push_back(l);
      }

   }                                              //
}


//////////////////////////////////////////////////////////////////////////
// PCS Execution
//////////////////////////////////////////////////////////////////////////

/*************************************************************************
ROCKFLOW - Function: CRFProcess::
Task:
Programming:
02/2003 OK Implementation
04/2003 OK Storing time step results
08/2004 OK PCS2
01/2005 WW new ELE concept
03/2005 OK MultiMSH concept
06/2005 MB NEWTON error calculation
05/2007 WW DOF>1, unerrelaxation for Picard, and removing the old or spurious stuff
12/2007 WW Classes of sparse matrix (jagged diagonal storage) and linear solver
and parellelisation of them
last modified:
**************************************************************************/
double CRFProcess::Execute()
{
   int nidx1;
   double pcs_error;
   double relax;                                  //WW
   int i;
   long j, k = 0, g_nnodes = 0;                   //07.01.07 WW

   //----------------------------------------------------------------------
#ifdef USE_MPI                                 //WW
   long global_eqs_dim = pcs_number_of_primary_nvals*m_msh->GetNodesNumber(false);
   CPARDomain *dom = dom_vector[myrank];
   if(myrank==0)
#endif
      cout << "    ->Process " << pcs_number << ": " << convertProcessTypeToString (this->getProcessType()) << endl;
   //	if (! this->_pcs_type_name.compare("MASS_TRANSPORT")) {
   if (this->getProcessType() == MASS_TRANSPORT)
   {
      cout << "      for " << this->pcs_primary_function_name[0];
      cout << " pcs_component_number " << this->pcs_component_number;
      cout << endl;
   }
   //----------------------------------------------------------------------
   // 0 Initialize equations
#ifdef NEW_EQS                                 //WW
   if(!configured_in_nonlinearloop)
#if defined(USE_MPI)
   {
                                                  //02.2010.  WW
      dom->eqs->SetDOF(pcs_number_of_primary_nvals);
      dom->ConfigEQS(m_num, global_eqs_dim);
   }
#else
   // Also allocate temporary memory for linear solver. WW
   {
                                                  //_new 02/2010. WW
      eqs_new->SetDOF(pcs_number_of_primary_nvals);
      eqs_new->ConfigNumerics(m_num);
   }
   eqs_new->Initialize();
#endif
#else
   // System matrix
   SetLinearSolverType(eqs, m_num);               //WW
   SetZeroLinearSolver(eqs);
#endif
   /*
    //TEST_MPI
    string test = "rank";
    char stro[1028];
    sprintf(stro, "%d",myrank);
    string test1 = test+(string)stro+"Assemble.txt";
    ofstream Dum(test1.c_str(), ios::out); // WW
    dom->eqs->Write(Dum);   Dum.close();
    MPI_Finalize();
    exit(1);
    */
   //----------------------------------------------------------------------
   // Solution vector
   relax = 0.0;
   if (m_num->nls_relaxation < DBL_MIN)
      relax = 0.0;
   else
      relax = 1.0 - m_num->nls_relaxation;        //WW
   //NW. should mark active elements if any process uses deactivation
   //if(NumDeactivated_SubDomains>0)
   //TODO if it's nonlinear, CheckMarkedElement() has been already called
   if (hasAnyProcessDeactivatedSubdomains)
#ifdef NEW_EQS                                 //WW
      if(!configured_in_nonlinearloop)
#endif
         CheckMarkedElement();
   m_msh->SwitchOnQuadraticNodes(false);
   for (i = 0; i < pcs_number_of_primary_nvals; i++)
   {
                                                  //new time
      nidx1 = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;
      g_nnodes = m_msh->GetNodesNumber(false);    //WW
      for (j = 0; j < g_nnodes; j++)              //WW
#ifdef NEW_EQS
         // New EQS  WW
         eqs_new->x[j+i*g_nnodes] = GetNodeValue(m_msh->Eqs2Global_NodeIndex[j],nidx1);
#else
      eqs->x[j + i * g_nnodes] = GetNodeValue(
         m_msh->Eqs2Global_NodeIndex[j], nidx1);
#endif
   }
   /*---------------------------------------------------------------------*/
   /* 1 Calc element matrices */
   if ((aktueller_zeitschritt == 1) || (tim_type_name.compare("TRANSIENT")
      == 0))
   {
#ifdef USE_MPI                              //WW
      if(myrank==0)
#endif
         cout << "      Calculate element matrices" << endl;
   }
   /*---------------------------------------------------------------------*/
   /* 2 Assemble EQS */
#ifdef USE_MPI                                 //WW
   if(myrank==0)
#endif
      cout << "      Assemble equation system" << endl;
   //
   // Assembly
#ifdef USE_MPI                                 //WW
   clock_t cpu_time=0;                            //WW
   cpu_time = -clock();
#endif
   //
   GlobalAssembly();
#ifdef USE_MPI
   cpu_time += clock();
   cpu_time_assembly += cpu_time;
#endif
   //
   //----------------------------------------------------------------------
   /* 5 Solve EQS */
   //orig  cout << "      Solve equation system" << endl;
   // cout << "Solve equation system myrank =" <<myrank<<"of size="<<size<<endl;
#ifdef CHECK_EQS
   std::string eqs_name = convertProcessTypeToString(this->getProcessType()) + "_EQS.txt";
   MXDumpGLS((char*)eqs_name.c_str(),1,eqs->b,eqs->x);
#endif
   //.....................................................................
   // Execute linnear solver
   //
#ifdef NEW_EQS                                 //WW
#if defined(USE_MPI)
   iter_lin=dom->eqs->Solver(eqs_new->x, global_eqs_dim);  //21.12.2007
#else
#ifdef LIS
   iter_lin=eqs_new->Solver(this->m_num);                  //NW
#else
   iter_lin=eqs_new->Solver();
#endif
#endif
#else
   iter_lin=ExecuteLinearSolver();
#endif
   //
   //PCSDumpModelNodeValues();
   //----------------------------------------------------------------------
   // Error calculation
   //----------------------------------------------------------------------
   double val_n = 0.;                             //03.04.2007.  WW
   // NEWTON
   if (m_num->nls_method_name.find("NEWTON") != string::npos)
   {
      double Val;
      pcs_error = 0.0;

      // MSH data
      nidx1 = GetNodeValueIndex("HEAD") + 1;
      //for(long nn=0;nn<NodeListSize(); nn++)  {
      for (long nn = 0; nn < (long) m_msh->nod_vector.size(); nn++)
      {
#ifdef NEW_EQS
         Val = GetNodeValue(nn,nidx1)+ eqs_new->x[nn];
#else
         Val = GetNodeValue(nn, nidx1) + eqs->x[nn];
#endif
         pcs_error = max(pcs_error, fabs(Val - GetNodeValue(nn, nidx1)));
         SetNodeValue(nn, nidx1, Val);
      }
#ifdef USE_MPI                              //WW
      if(myrank==0)
#endif
                                                  //OK7
         cout << "      PCS error: " << pcs_error << endl;
   }                                              // NEWTON
   //----------------------------------------------------------------------
   else                                           //PICARD
   {
      //......................................................................
      pcs_error = CalcIterationNODError(1);       //OK4105//WW4117
#ifdef USE_MPI                              //WW
      if(myrank==0)
#endif
         cout << "      PCS error: " << pcs_error << endl;
      //--------------------------------------------------------------------
      // 7 Store solution vector in model node values table
      //....................................................................
      for (i = 0; i < pcs_number_of_primary_nvals; i++)
      {
         nidx1 = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;
         g_nnodes = m_msh->GetNodesNumber(false); //DOF>1, WW
         for (j = 0; j < g_nnodes; j++)
         {
            k = m_msh->Eqs2Global_NodeIndex[j];
            val_n = GetNodeValue(k, nidx1);       //03.04.2009. WW
#ifdef NEW_EQS
            SetNodeValue(k,nidx1,relax*val_n+(1.0-relax)*eqs_new->x[j+i*g_nnodes]);
            eqs_new->x[j+i*g_nnodes] = val_n;     // Used for time stepping. 03.04.2009. WW
#else
            SetNodeValue(k, nidx1, relax * val_n + (1.0 - relax) * eqs->x[j
               + i * g_nnodes]);
            eqs->x[j + i * g_nnodes] = val_n;     // Used for time stepping. 03.04.2009. WW
#endif
         }
      }
   }                                              // END PICARD
   //----------------------------------------------------------------------
#ifdef NEW_EQS                                 //WW
   if(!configured_in_nonlinearloop)
#if defined(USE_MPI)
      dom->eqs->Clean();
#else
   // Also allocate temporary memory for linear solver. WW
   eqs_new->Clean();
#endif
#endif
   return pcs_error;
}


/*************************************************************************
GEOSYS - Function:
Task:
Programming:
 08/2008 WW Implementation
 11/2008 WW Update
last modified:
**************************************************************************/
void CRFProcess::CopyU_n(double *temp_v)
{
   int i, nidx1;
   long g_nnodes, j, k;
   for(i=0; i<pcs_number_of_primary_nvals; i++)
   {
      nidx1 = GetNodeValueIndex(pcs_primary_function_name[i])+1;
      g_nnodes =m_msh->GetNodesNumber(false);     //DOF>1, WW
      for(j=0;j<g_nnodes;j++)
      {
         k = m_msh->Eqs2Global_NodeIndex[j];
         temp_v[j+i*g_nnodes] = GetNodeValue(k,nidx1);
      }
   }
}


/*************************************************************************
ROCKFLOW - Function:
Task: Initialize the equation system
Programming:
 05/2003 OK/WW Implementation
 08/2003 WW    Changes for momolithic coupled equations
last modified:
**************************************************************************/
#ifndef NEW_EQS                                   //WW 07.11.2008
void CRFProcess::InitEQS()
{
   int i;
   eqs->unknown_vector_indeces = \
      (int*) Malloc(pcs_number_of_primary_nvals*sizeof(int));
   eqs->unknown_node_numbers = \
      (long*) Malloc(pcs_number_of_primary_nvals*sizeof(long));
   eqs->unknown_update_methods = \
      (int*) Malloc(pcs_number_of_primary_nvals*sizeof(long));

   for(i=0;i<pcs_number_of_primary_nvals;i++)
   {
      eqs->unknown_vector_indeces[i] = GetNodeValueIndex(pcs_primary_function_name[i])+1;
      eqs->unknown_node_numbers[i] = m_msh->GetNodesNumber(true);
      eqs->unknown_update_methods[i] = 2;         //For Newton-Raphson method. WW
      if(type==41)                                // Monolithic scheme for p-u coupled equation
      {
         if(i==pcs_number_of_primary_nvals-1)
            eqs->unknown_node_numbers[i] = m_msh->GetNodesNumber(false);
         eqs->unknown_update_methods[i] = 0;
      }
   }
}
#endif
/*************************************************************************
ROCKFLOW - Function:
Task: Calculate element matrices
Programming:
05/2003 OK Implementation
09/2005 OK gas flow removed
last modified:
**************************************************************************/
void CRFProcess::CalculateElementMatrices(void)
{
   switch(this->type)
   {
      case 1:                                     //SM
         break;
      case 2:                                     //MTM2
         break;
      case 3:                                     //HTM
         break;
      case 5:                                     // Gas flow
         break;
      case 11:
         break;
      case 12:                                    //MMP
         break;
      case 13:                                    //MPC
         break;
      case 66:                                    //OF
         break;
      default:
         DisplayMsgLn("CalculateElementMatrices: no CalculateElementMatrices specified");
         abort();
   }
}


/*************************************************************************
GeoSys-Function:
Task: Assemble the global system equation
Programming:
01/2005 WW/OK Implementation
04/2005 OK MSH
07/2005 WW Change due to the geometry element objects applied
10/2005 OK DDC
11/2005 YD time step control
01/2006 OK/TK Tests
12/2007 WW Spase matrix class and condensation sequential and parallel loop
last modified:
10/2010 TF changed access to process type
**************************************************************************/
void CRFProcess::GlobalAssembly()
{
   // Tests
   if (!Tim)
      Tim = TIMGet(convertProcessTypeToString(this->getProcessType()));
   if (!Tim)
   {
      cout << "Error in CRFProcess::GlobalAssembly() - no TIM data" << endl;
      return;
   }

   if (!fem)
   {
      // Which process needs this?
      // Only one instance of CFiniteElementStd is required for each process
      // Use "new" in such way will cause memory problem.
      // Please move this declaration to pcs configuration. 	  WW
      if (m_msh)
         fem = new CFiniteElementStd(this, m_msh->GetCoordinateFlag());
   }

   CElem* elem = NULL;
   bool Check2D3D;
   Check2D3D = false;
   if (type == 66)                                //Overland flow
      Check2D3D = true;

   // DDC
   if (dom_vector.size() > 0)
   {
      cout << "      Domain Decomposition" << '\n';
      CPARDomain* m_dom = NULL;
      size_t j = 0;
#if defined(USE_MPI)                        //WW
      j = myrank;
#else
      for (j = 0; j < dom_vector.size(); j++)
      {
#endif
         m_dom = dom_vector[j];
#ifdef NEW_EQS
         m_dom->InitialEQS(this);
#else
         SetLinearSolver(m_dom->eqs);
         SetZeroLinearSolver(m_dom->eqs);
#endif
         for (size_t i = 0; i < m_dom->elements.size(); i++)
         {
            elem = m_msh->ele_vector[m_dom->elements[i]];
            if (elem->GetMark())
            {
               elem->SetOrder(false);
                                                  //WW
               fem->SetElementNodesDomain(m_dom->element_nodes_dom[i]);
               fem->ConfigElement(elem, Check2D3D);
               fem->m_dom = m_dom;                //OK
               fem->Assembly();
            }
         }
         // m_dom->WriteMatrix();
         //MXDumpGLS("rf_pcs.txt",1,m_dom->eqs->b,m_dom->eqs->x);
         //ofstream Dum("rf_pcs.txt", ios::out); // WW
         //m_dom->eqs->Write(Dum);
         //Dum.close();
         IncorporateSourceTerms(j);
         IncorporateBoundaryConditions(j);
         /*
          //TEST
          string test = "rank";
          char stro[64];
          sprintf(stro, "%d",j);
          string test1 = test+(string)stro+"Assemble.txt";
          ofstream Dum(test1.c_str(), ios::out);
          m_dom->eqs->Write(Dum);
          Dum.close();
          exit(1);
          */
#ifndef USE_MPI
      }
      // Assemble global system
      DDCAssembleGlobalMatrix();
#endif
   }
   else
   {                                              // STD
                                                  //YDTEST. Changed to DOF 15.02.2007 WW
      for (size_t ii = 0; ii < continuum_vector.size(); ii++)
      {
         continuum = ii;
         for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
         {
            elem = m_msh->ele_vector[i];
            if (elem->GetMark())                  // Marked for use
            {
               elem->SetOrder(false);
               fem->ConfigElement(elem, Check2D3D);
               fem->Assembly();
               // NEUMANN CONTROL---------
               if (Tim->time_control_name.compare("NEUMANN") == 0)
               {
                  Tim->time_step_length_neumann = MMin(
                     Tim->time_step_length_neumann, timebuffer);
                  Tim->time_step_length_neumann *= 0.5
                     * elem->GetVolume() * elem->GetVolume();
                  if (Tim->time_step_length_neumann < MKleinsteZahl)
                     Tim->time_step_length_neumann = 1.0e-5;
               }
               //------------------------------
            }
         }
      }
      //MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); //abort();
      //eqs_new->Write();
      IncorporateSourceTerms();
#ifdef GEM_REACT
      //		if ( _pcs_type_name.compare("MASS_TRANSPORT") == 0 && aktueller_zeitschritt > 1 && this->m_num->cpl_iterations > 1)
      if ( this->getProcessType() == MASS_TRANSPORT && aktueller_zeitschritt > 1 && this->m_num->cpl_iterations > 1)
         IncorporateSourceTerms_GEMS();
#endif
#ifndef NEW_EQS                             //WW. 07.11.2008
      SetCPL();                                   //OK
#endif
      IncorporateBoundaryConditions();
      // ofstream Dum("rf_pcs.txt", ios::out); // WW
      // eqs_new->Write(Dum);   Dum.close();
      //
      //
      // MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); //abort();
      //          MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); //abort();
   }
}


/*************************************************************************
GeoSys-Function:
Task: Integration
Programming:
05/2009 WW Implementation
**************************************************************************/
void CRFProcess::Integration(vector<double> &node_velue)
{
   //----------------------------------------------------------------------
   int k;
   long i;
   CElem* elem = NULL;
   bool Check2D3D;
   Check2D3D = false;
   double n_val[8];

   if(type == 66)                                 //Overland flow
      Check2D3D = true;

   vector<double> buffer((long)node_velue.size());
   for(i=0; i<(long)buffer.size(); i++)
      buffer[i] = 0.;

   for(i=0;i<(long)m_msh->ele_vector.size();i++)
   {
      elem = m_msh->ele_vector[i];
      if(!elem->GetMark()) continue;

      for(k=0; k<elem->GetNodesNumber(false); k++)
         n_val[k] = node_velue[elem->GetNodeIndex(k)];

      elem->SetOrder(false);
      fem->ConfigElement(elem,Check2D3D);
      fem->FaceIntegration(n_val);

      for(k=0; k<elem->GetNodesNumber(false); k++)
         buffer[elem->GetNodeIndex(k)] += n_val[k];

   }
   //----------------------------------------------------------------------
}


/*************************************************************************
GeoSys-Function:
Task: Calculate integration point velocity
Programming:
08/2005 WW Implementation
last modified:
**************************************************************************/
void CRFProcess::CalIntegrationPointValue()
{
   CElem* elem = NULL;
   cal_integration_point_value = false;
   continuum = 0;                                 //15.02.2007/
   //  cal_integration_point_value = true;
   // Currently, extropolation only valid for liquid and Richards flow.
   //	if (_pcs_type_name.find("LIQUID") != string::npos || _pcs_type_name.find(
   //			"RICHARD") != string::npos || _pcs_type_name.find(
   //			"MULTI_PHASE_FLOW") != string::npos || _pcs_type_name.find(
   //			"GROUNDWATER_FLOW") != string::npos || _pcs_type_name.find(
   //			"TWO_PHASE_FLOW") != string::npos
   //			|| _pcs_type_name.find("AIR_FLOW") != string::npos
   //			|| _pcs_type_name.find("PS_GLOBAL") != string::npos) //WW/CB
   if (getProcessType() == LIQUID_FLOW || getProcessType() == RICHARDS_FLOW
      || getProcessType() == MULTI_PHASE_FLOW
      || getProcessType() == GROUNDWATER_FLOW
      || getProcessType() == TWO_PHASE_FLOW
      || getProcessType() == AIR_FLOW
      || getProcessType() == PS_GLOBAL)           //WW/CB/TF
      cal_integration_point_value = true;
   if (!cal_integration_point_value)
      return;
   const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
   for (size_t i = 0; i < mesh_ele_vector_size; i++)
   {
      elem = m_msh->ele_vector[i];
      if (elem->GetMark())                        // Marked for use
      {
         fem->ConfigElement(elem);
         fem->Config();                           //OK4709
         // fem->m_dom = NULL; // To be used for parallization
         fem->Cal_Velocity();
      }
   }
   //	if (_pcs_type_name.find("TWO_PHASE_FLOW") != string::npos) //WW/CB
   if (getProcessType() == TWO_PHASE_FLOW)        //WW/CB
      cal_integration_point_value = false;
}


/*************************************************************************
GeoSys-Function: CalGPVelocitiesfromFluidMomentum
Task: Calculate gauss point velocities from fluid momentum solution
      extrapolate velocitiues from nodes to gauss points
Programming:
09/2009 SB Implementation

**************************************************************************/
void CRFProcess::CalGPVelocitiesfromFluidMomentum()
{
   long i;
   Mesh_Group::CElem* elem = NULL;
   int i_ind[3];

   cout << "      CalGPVelocitiesfromFluidMomentum()" << endl;

   // Get fluid_momentum process
   CRFProcess *m_pcs_fm =  PCSGet(FLUID_MOMENTUM);

   //  check all possibilities for grid orientation (ccord_flag)
   int coordinateflag = this->m_msh->GetCoordinateFlag();
                                                  // get index of velocity
   i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_X")+1;
   if(coordinateflag == 11)
                                                  // get index of velocity
      i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Y")+1;
   if(coordinateflag == 12)
                                                  // get index of velocity
      i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Z")+1;
                                                  // get index of velocity
   i_ind[1] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Y")+1;
   if(coordinateflag == 22)
                                                  // get index of velocity
      i_ind[1] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Z")+1;
                                                  // get index of velocity
   i_ind[2] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Z")+1;

   if((i_ind[0]<0)||(i_ind[1]<0)||(i_ind[2]<0))
      cout << " Error - wrong index in Cal_GP_Velocity_FM " << endl;

   // Loop over all elements
   for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
   {
      elem = m_msh->ele_vector[i];                // get element
      if (elem->GetMark())                        // Marked for use
      {
         fem->ConfigElement(elem);
         fem->Cal_GP_Velocity_FM(i_ind);
      }
   }                                              // end element loop

}


/*************************************************************************
ROCKFLOW - Function: AllocateLocalMatrixMemory
Task: As the function name
Programming:
01/2005 WW/OK Implementation
06/2005 OK MSH
last modified:
**************************************************************************/
void CRFProcess::AllocateLocalMatrixMemory()
{
   long i;
   //----------------------------------------------------------------------
   int up_type = 0;
   if(!M_Process) up_type=0;
   if(H_Process&&M_Process)
   {
      if(type!=4&&type!=41)
         up_type=1;
      else
      {
         if(type==4)
            up_type = 3;
         if(type==41)  up_type = 4;
      }
   }
   if(!H_Process) up_type = 2;
   if(MASS_TRANSPORT_Process)up_type = 5;         //SB for steady state element matrices in transport
   //----------------------------------------------------------------------
   ElementMatrix *eleMatrix = NULL;
   CElem* elem = NULL;
   for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
   {
      elem = m_msh->ele_vector[i];
      if (elem->GetMark())                        // Marked for use
      {
         eleMatrix = new ElementMatrix();
         eleMatrix->AllocateMemory(elem,up_type);
         Ele_Matrices.push_back(eleMatrix);
      }
   }
}


/*************************************************************************
FEMLib function
Task: Assemble global system matrix
Programming:
05/2003 OK Implementation
??/???? WW Moved from AssembleSystemMatrixNew
05/2006 WW Modified to enable dealing with the case of DOF>1
06/2006 WW Take the advantege of sparse matrix to enhance simulation
10/2007 WW Change for the new classes of sparse matrix and linear solver
**************************************************************************/
void CRFProcess::DDCAssembleGlobalMatrix()
{
   int ii,jj, dof;
   long i,j,j0,ig,jg, ncol;
   CPARDomain *m_dom = NULL;
   long *nodes2node = NULL;                       //WW
   double *rhs = NULL, *rhs_dom = NULL;
   double a_ij;
   double b_i=0.0;
   b_i = b_i;                                     //OK411
   int no_domains =(int)dom_vector.size();
   long no_dom_nodes;
   dof = pcs_number_of_primary_nvals;             //WW
   ncol = 0;                                      //WW
#ifndef USE_MPI
   int k;
   for(k=0;k<no_domains;k++)
   {
      m_dom = dom_vector[k];
#else
      m_dom = dom_vector[myrank];
#endif
      // RHS
#if defined(NEW_EQS)
      rhs = eqs_new->b;
      if(type==4)
         rhs_dom = m_dom->eqsH->b;
      else
         rhs_dom = m_dom->eqs->b;
#else
      rhs = eqs->b;
      rhs_dom = m_dom->eqs->b;
#endif

      no_dom_nodes = m_dom->nnodes_dom;           //WW
      if(type==4||type==41)
         no_dom_nodes = m_dom->nnodesHQ_dom;      //WW
      if(type==41) dof--;
      for(i=0;i<no_dom_nodes;i++)
      {
         //------------------------------------------
         // Use the feature of sparse matrix of FEM
         //WW
         ig = m_dom->nodes[i];
         ncol = m_dom->num_nodes2_node[i];
         nodes2node =  m_dom->node_conneted_nodes[i];
         for(j0=0;j0<ncol;j0++)
         {
            j=nodes2node[j0];
            if(j>=no_dom_nodes) continue;
            jg = m_dom->nodes[j];
            //------------------------------------------
            // DOF loop ---------------------------WW
            for(ii=0; ii<dof; ii++)
            {
               for(jj=0; jj<dof; jj++)
               {
                  // get domain system matrix
#ifdef  NEW_EQS                 //WW
                  if(type==4)
                     a_ij = (*m_dom->eqsH->A)(i+no_dom_nodes*ii,j+no_dom_nodes*jj);
                  else
                     a_ij = (*m_dom->eqs->A)(i+no_dom_nodes*ii,j+no_dom_nodes*jj);
                  (*eqs_new->A)(ig+Shift[ii],jg+Shift[jj]) += a_ij;
#else
                  SetLinearSolver(m_dom->eqs);
                  a_ij = MXGet(i+no_dom_nodes*ii,j+no_dom_nodes*jj);
                  // set global system matrix
                  SetLinearSolver(eqs);
                  MXInc(ig+Shift[ii],jg+Shift[jj],a_ij);
#endif
               }
            }
            // DOF loop ---------------------------WW
         }
         // set global RHS vector //OK
         for(ii=0; ii<dof; ii++)                  //WW
            rhs[ig+Shift[ii]] += rhs_dom[i+no_dom_nodes*ii];
      }

      // Mono HM------------------------------------WW
      if(type!=41)
#ifndef USE_MPI
         continue;
#else
      return;
#endif
      no_dom_nodes = m_dom->nnodes_dom;
      long  no_dom_nodesHQ = m_dom->nnodesHQ_dom;
      double a_ji = 0.0;
      for(i=0;i<no_dom_nodes;i++)
      {
         ig = m_dom->nodes[i];                    //WW
         ncol = m_dom->num_nodes2_node[i];
         nodes2node =  m_dom->node_conneted_nodes[i];
         for(j0=0;j0<ncol;j0++)
         {
            j=nodes2node[j0];
            jg = m_dom->nodes[j];
            for(ii=0; ii<dof; ii++)               //ww
            {
#if defined(NEW_EQS)
               // dom to global. WW
               a_ij = (*m_dom->eqsH->A)(i+no_dom_nodesHQ*dof,j+no_dom_nodesHQ*ii);
               a_ji = (*m_dom->eqsH->A)(j+no_dom_nodesHQ*ii, i+no_dom_nodesHQ*dof);
               (*eqs_new->A)(ig+Shift[ii],jg+Shift[problem_dimension_dm]) += a_ij;
               (*eqs_new->A)(jg+Shift[problem_dimension_dm],ig+Shift[ii]) += a_ji;
#else
               // get domain system matrix
               SetLinearSolver(m_dom->eqs);
               a_ij = MXGet(i+no_dom_nodesHQ*dof,j+no_dom_nodesHQ*ii);
               a_ji = MXGet(j+no_dom_nodesHQ*ii, i+no_dom_nodesHQ*dof);
               // set global system matrix
               SetLinearSolver(eqs);
               MXInc(ig+Shift[ii],jg+Shift[problem_dimension_dm],a_ij);
               MXInc(jg+Shift[problem_dimension_dm],ig+Shift[ii],a_ji);
#endif
            }
         }
      }
      for(i=0;i<no_dom_nodes;i++)
      {
         ig = m_dom->nodes[i];
         ncol = m_dom->num_nodes2_node[i];
         nodes2node =  m_dom->node_conneted_nodes[i];
         for(j0=0;j0<ncol;j0++)
         {
            j=nodes2node[j0];
            jg = m_dom->nodes[j];
            if(jg>=no_dom_nodes) continue;
            // get domain system matrix
#if defined(NEW_EQS)
            // dom to global. WW
            a_ij = (*m_dom->eqsH->A)(i+no_dom_nodesHQ*dof,j+no_dom_nodesHQ*dof);
            (*eqs_new->A)(ig+Shift[problem_dimension_dm],jg+Shift[problem_dimension_dm])
               += a_ij;
#else
            SetLinearSolver(m_dom->eqs);
            a_ij = MXGet(i+no_dom_nodesHQ*dof,j+no_dom_nodesHQ*dof);
            // set global system matrix
            SetLinearSolver(eqs);
            MXInc(ig+Shift[problem_dimension_dm],jg+Shift[problem_dimension_dm],a_ij);
#endif
         }
         //
         rhs[ig+Shift[problem_dimension_dm]] += rhs_dom[i+no_dom_nodesHQ*dof];
      }
      // Mono HM------------------------------------WW
#ifndef USE_MPI
   }
#endif

}


/*************************************************************************
ROCKFLOW - Function:
Task: Assemble system matrix
Programming: 05/2003 OK Implementation
ToDo: Prototyp function
last modified:
**************************************************************************/
void CRFProcess::AssembleSystemMatrixNew(void)
{
   switch(type)
   {
      case 1:
         //MakeGS_ASM_NEW(eqs->b,eqs->x,ddummy);
         // SMAssembleMatrix(eqs->b,eqs->x,ddummy,this);
         break;
      case 2:                                     //MTM2
         break;
      case 3:                                     //HTM
         break;
      case 5:                                     // Gas flow
         break;
      case 11:
         break;
      case 12:
         break;
      case 13:
         break;
      case 66:
         //MakeGS_ASM_NEW(eqs->b,eqs->x,ddummy);
         // SMAssembleMatrix(eqs->b,eqs->x,ddummy,this);
         break;
      default:
         DisplayMsgLn("CalculateElementMatrices: no CalculateElementMatrices specified");
         abort();
   }
#ifdef PARALLEL
   DDCAssembleGlobalMatrix();
#else
   IncorporateSourceTerms();
   IncorporateBoundaryConditions();
#endif
   //SetLinearSolver(eqs);
   //MXDumpGLS("global_matrix_dd.txt",1,eqs->b,eqs->x);
}


/**************************************************************************
FEMLib-Method: CRFProcess::IncorporateBoundaryConditions
Task: set PCS boundary conditions
Programing:
05/2006 WW Implementation
**************************************************************************/
void CRFProcess::SetBoundaryConditionSubDomain()
{
   int k;
   long i,j;
   CPARDomain *m_dom = NULL;
   CBoundaryConditionNode *m_bc_nv = NULL;
   CNodeValue *m_st_nv = NULL;
   //
   for(k=0;k<(int)dom_vector.size();k++)
   {
      m_dom = dom_vector[k];
      // BC
      for(i=0; i<(long)bc_node_value.size(); i++)
      {
         m_bc_nv = bc_node_value[i];
         for(j=0; j<(long)m_dom->nodes.size(); j++)
         {
            if(m_bc_nv->geo_node_number==m_dom->nodes[j])
            {
               bc_node_value_in_dom.push_back(i);
               bc_local_index_in_dom.push_back(j);
               break;
            }
         }
      }
      rank_bc_node_value_in_dom.push_back((long)bc_node_value_in_dom.size());
      // ST
      for(i=0; i<(long)st_node_value.size(); i++)
      {
         m_st_nv = st_node_value[i];
         for(j=0; j<(long)m_dom->nodes.size(); j++)
         {
            if(m_st_nv->geo_node_number==m_dom->nodes[j])
            {
               st_node_value_in_dom.push_back(i);
               st_local_index_in_dom.push_back(j);
               break;
            }
         }
      }
      rank_st_node_value_in_dom.push_back((long)st_node_value_in_dom.size());
   }
   long Size = (long)st_node_value.size();
   long l_index;
   for(i=0; i<Size; i++)
   {
      l_index = st_node_value[i]->geo_node_number;
      st_node_value[i]->node_value /= (double)node_connected_doms[l_index];
   }

}


/**************************************************************************
FEMLib-Method: CRFProcess::SetSTWaterGemSubDomain
Task: set source/sink terms for GEMS-flow coupling
Programing:
05/2006 WW Implementation
03/2010 KG44 modified to GEM
**************************************************************************/
void CRFProcess::SetSTWaterGemSubDomain(int myrank)
{
   int k;
   long i,j, dsize=0;
   CPARDomain *m_dom = NULL;
   long int m_stgem_nv = -1;
   //
   long Size = ( long ) Water_ST_vec.size();
   long l_index=-1;

   //	cout << "dom_vec_size: " << dom_vector.size() << endl;
   //	for ( k=0;k< ( int ) dom_vector.size();k++ )
   //	{
   k=myrank;                                      //do it for each domain only once!
   m_dom = dom_vector[k];
   dsize=(long) m_dom->nodes.size();
   // ST
   for ( i=0; i< Size; i++ )
   {
      m_stgem_nv = Water_ST_vec[i].index_node;
      for ( j=0; j< ( long ) m_dom->nodes.size(); j++ )
      {
         if ( m_stgem_nv==m_dom->nodes[j] )
         {
                                                  //index for Water_ST_vec
            stgem_node_value_in_dom.push_back ( i );
                                                  //index for RHS
            stgem_local_index_in_dom.push_back ( j );
            //	cout << "dom " << k <<  " i, j " << i << " " << j   << endl;
         }
      }
   }
                                                  // only one element per domain!
   rank_stgem_node_value_in_dom.push_back ( ( long ) stgem_node_value_in_dom.size() );
   //	cout << "dom " << k <<  " rank_stgem_node_value_in_dom " << (long) rank_stgem_node_value_in_dom[0]  << endl;

   //	}

   for ( i=0; i<Size; i++ )
   {
      l_index = Water_ST_vec[i].index_node;
      // cout << i << " " << node_connected_doms[l_index] << " " << endl;
                                                  //values for shared nodes are scaled
      Water_ST_vec[i].water_st_value /= ( double ) node_connected_doms[l_index];
   }

   return;
}


/**************************************************************************
FEMLib-Method: CRFProcess::IncorporateBoundaryConditions
Task: set PCS boundary conditions
Programing:
02/2004 OK Implementation
????    WW  and time curve
04/2005 OK MSH
05/2005 OK conditional BCs
04/2006 WW Changes due to the geometry object applied
04/2006 OK Conditions by PCS coupling OK
05/2006 WW Re-implement
05/2006 WW DDC
10/2007 WW Changes for the new classes of sparse matrix and linear solver
last modification:
**************************************************************************/
void CRFProcess::IncorporateBoundaryConditions(const int rank)
{
   static long i;
   static double bc_value, fac = 1.0, time_fac = 1.0;
   long bc_msh_node=-1;
   long bc_eqs_index, shift;
   int interp_method=0;
   int curve, valid=0;
   int ii, idx0, idx1;
   CBoundaryConditionNode* m_bc_node;             //WW
   CBoundaryCondition* m_bc;                      //WW
   CPARDomain *m_dom = NULL;
   CFunction* m_fct = NULL;                       //OK
   double *eqs_rhs = NULL;
   bool is_valid = false;                         //OK
#ifdef NEW_EQS
   Linear_EQS *eqs_p = NULL;
#endif
   //------------------------------------------------------------WW
   // WW
   double Scaling = 1.0;
   bool quadr = false;                            //15.4.2008. WW
   if(type==4||type==41)
   {
      fac = Scaling;
      quadr = true;
   }
   long begin = 0;
   long end = 0;
   long gindex=0;
   if(rank==-1)
   {
      begin = 0;
      end = (long)bc_node_value.size();
#ifdef NEW_EQS                              //WW
      eqs_p = eqs_new;
      eqs_rhs = eqs_new->b;                       //27.11.2007 WW
#else
      eqs_rhs = eqs->b;
#endif
   }
   else
   {
      m_dom = dom_vector[rank];
#ifdef NEW_EQS
      eqs_p = m_dom->eqs;
      if(type == 4 )                              //WW
      {
         eqs_p = m_dom->eqsH;
         eqs_rhs = m_dom->eqsH->b;
      }
      else
         eqs_rhs = m_dom->eqs->b;
#else
      eqs_rhs = m_dom->eqs->b;
#endif
      if(rank==0)
         begin = 0;
      else
         begin = rank_bc_node_value_in_dom[rank-1];
      end = rank_bc_node_value_in_dom[rank];
   }

   for(i=begin;i<end;i++)
   {
      gindex = i;
      if(rank>-1)
         gindex = bc_node_value_in_dom[i];
      m_bc_node = bc_node_value[gindex];
      m_bc = bc_node[gindex];
      shift = m_bc_node->msh_node_number-m_bc_node->geo_node_number;
      //
      if(rank>-1)
      {
         bc_msh_node = bc_local_index_in_dom[i];
         int dim_space = 0;
         if(shift==0)
            // 15.4.2008 WW if(m_msh->NodesNumber_Linear==m_msh->NodesNumber_Quadratic)
            dim_space = 0;
         else
         {
            // Following lines are changed since 15.4.2008 by WW
            if(quadr)
            {
               if(shift%m_msh->NodesNumber_Quadratic==0)
                  dim_space = shift/m_msh->NodesNumber_Quadratic;
               else
                  dim_space = m_msh->msh_max_dim;
            }
            else
            {
               if(shift%m_msh->NodesNumber_Linear==0)
                  dim_space = shift/m_msh->NodesNumber_Linear;
               else
                  dim_space = m_msh->msh_max_dim;
            }
         }
         shift = m_dom->shift[dim_space];
      }
      else
         bc_msh_node = m_bc_node->geo_node_number;
      //------------------------------------------------------------WW
      if(m_msh)                                   //OK
         //	    if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
         //          continue;
         time_fac = 1.0;
      if(bc_msh_node>=0)
      {
         //................................................................
         // Time dependencies - CURVE
         curve =  m_bc_node->CurveIndex;
         if(curve>0)
         {
            if(curve>10000000)                    /// 16.08.2010. WW
               time_fac = GetCurveValue(curve-10000000,interp_method,aktuelle_zeit,&valid);
            else
               time_fac = GetCurveValue(curve,interp_method,aktuelle_zeit,&valid);
            if(!valid) continue;
         }
         else
            time_fac = 1.0;
         //................................................................
         // Time dependencies - FCT
         if(m_bc_node->fct_name.length()>0)
         {
            m_fct = FCTGet(m_bc_node->fct_name);
            if(m_fct)
            {
               time_fac = m_fct->GetValue(aktuelle_zeit,&is_valid);
               //if(!valid) continue;
            }
            else
               cout << "Warning in CRFProcess::IncorporateBoundaryConditions - no FCT data" << endl;
         }
         //................................................................
         // Conditions
         if(m_bc_node->conditional)
         {
            int idx_1 = -1;                       //28.2.2007 WW
            for(ii=0;ii<dof;ii++)                 //28.2.2007 WW
            {
               if (convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()).find(pcs_primary_function_name[ii]) != string::npos)
               {
                  idx_1 = GetNodeValueIndex(pcs_primary_function_name[ii])+1;
                  break;
               }
            }
            bc_value = time_fac*fac* GetNodeValue(m_bc_node->msh_node_number_subst, idx_1);
            //WW  bc_value = time_fac*fac* GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1); // YD-----TEST---
         }
         else
                                                  // time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
            bc_value = time_fac*fac*m_bc_node->node_value;

         /// 16.08.2010. WW
         if(curve>10000000&&fabs(time_fac)>DBL_EPSILON)
            bc_value = bc_value/time_fac+time_fac;/// bc_value +time_fac;

         if(m_bc->isPeriodic())                   // JOD
            bc_value *= sin( 2 * 3.14159 * aktuelle_zeit / m_bc->getPeriodeTimeLength()  + m_bc->getPeriodePhaseShift());

         //----------------------------------------------------------------
         // MSH
         if(rank>-1)
            bc_eqs_index = bc_msh_node;
         else
                                                  //WW#
            bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
         //..............................................................
         // NEWTON WW
         if(m_num->nls_method_name.find("NEWTON")!=string::npos
            ||  type==4||type==41  )              //Solution is in the manner of increment !
         {
            idx0 = GetNodeValueIndex(convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()).c_str());
            if(type==4||type==41)
            {
                                                  //02.09.2010. WW
               if(m_bc_node->pcs_pv_name.find("PRESSURE")==string::npos)
               {
                  idx1 = idx0+1;
                  bc_value -=  GetNodeValue(m_bc_node->geo_node_number, idx0)
                     + GetNodeValue(m_bc_node->geo_node_number, idx1);
               }
            }
            else
               bc_value = bc_value - GetNodeValue(m_bc_node->geo_node_number, idx0);
         }
         //----------------------------------------------------------------
         bc_eqs_index += shift;
         /* // Make the follows as comment by WW. 04.03.2008
         //YD dual
         if(dof>1) //WW
         {
           for(ii=0;ii<dof;ii++)
           {
             if(m_bc->pcs_pv_name.find(pcs_primary_function_name[ii]) != string::npos)
             {
                //YD/WW
                //WW   bc_eqs_index += ii*(long)eqs->dim/dof;   //YD dual //DOF>1 WW
                bc_eqs_index += fem->NodeShift[ii];   //DOF>1 WW

         break;
         }
         }
         }
         */
#ifdef NEW_EQS                           //WW
         eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#else
         MXRandbed(bc_eqs_index,bc_value,eqs_rhs);
#endif
      }
   }
}


/**************************************************************************
FEMLib-Method: CRFProcess::IncorporateBoundaryConditions
Task: set PCS boundary conditions for FLUID_MOMENTUM depending on axis
Programing:
01/2007 PCH Implementation
last modification:
**************************************************************************/
void CRFProcess::IncorporateBoundaryConditions(const int rank, const int axis)
{
   static long i;
   static double bc_value, fac = 1.0, time_fac = 1.0;
   long bc_msh_node;
   long bc_eqs_index, shift;
   int interp_method=0;
   int curve, valid=0;
   int idx0, idx1;
   CBoundaryConditionNode* m_bc_node;             //WW
   CBoundaryCondition* m_bc;                      //WW
   CPARDomain *m_dom = NULL;
   CFunction* m_fct = NULL;                       //OK
   double *eqs_rhs = NULL;
   bool is_valid = false;                         //OK
#ifdef NEW_EQS
   Linear_EQS *eqs_p = NULL;
#endif
   //------------------------------------------------------------WW
   // WW
   double Scaling = 1.0;
   if(type==4||type==41) fac = Scaling;

   long begin = 0;
   long end = 0;
   long gindex=0;

   CBoundaryConditionsGroup *m_bc_group = NULL;
   //  m_bc_group = BCGetGroup(this->_pcs_type_name,this->pcs_primary_function_name[axis]);
                                                  // TF
   m_bc_group = BCGetGroup(convertProcessTypeToString(this->getProcessType()),this->pcs_primary_function_name[axis]);

   if(rank==-1)
   {
      begin = 0;
      end = (long)bc_node_value.size();
#ifdef NEW_EQS                              //WW
      eqs_p = eqs_new;
      eqs_rhs = eqs_new->b;                       //27.11.2007 WW
#else
      eqs_rhs = eqs->b;
#endif
   }
   else
   {
      m_dom = dom_vector[rank];
#ifdef NEW_EQS
      eqs_p = m_dom->eqs;
      if(type == 4 )                              //WW
      {
         eqs_p = m_dom->eqsH;
         eqs_rhs = m_dom->eqsH->b;
      } else
      eqs_rhs = m_dom->eqs->b;
#else
      eqs_rhs = m_dom->eqs->b;
#endif
      if(rank==0)
         begin = 0;
      else
         begin = rank_bc_node_value_in_dom[rank-1];
      end = rank_bc_node_value_in_dom[rank];
   }

   for(i=begin;i<end;i++)
   {
      gindex = i;
      if(rank>-1)
         gindex = bc_node_value_in_dom[i];
      m_bc_node = bc_node_value[gindex];

                                                  // PCH
      if(axis == 0 && m_bc_node->pcs_pv_name.find("VELOCITY1_X")!=string::npos)
      {
         m_bc = bc_node[gindex];
         shift = m_bc_node->msh_node_number-m_bc_node->geo_node_number;
         //
         if(rank>-1)
         {
            bc_msh_node = bc_local_index_in_dom[i];
            int dim_space = 0;
            if(m_msh->NodesNumber_Linear==m_msh->NodesNumber_Quadratic)
               dim_space = 0;
            else
            {
               if(shift%m_msh->NodesNumber_Quadratic==0)
                  dim_space = shift/m_msh->NodesNumber_Quadratic;
               else
                  dim_space = m_msh->msh_max_dim;
            }
            shift = m_dom->shift[dim_space];
         }
         else
            bc_msh_node = m_bc_node->geo_node_number;
         //------------------------------------------------------------WW
         if(m_msh)                                //OK
            //			if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
            //				continue;
            time_fac = 1.0;
         if(bc_msh_node>=0)
         {
            //................................................................
            // Time dependencies - CURVE
            curve =  m_bc_node->CurveIndex;
            if(curve>0)
            {
               time_fac = GetCurveValue(curve,interp_method,aktuelle_zeit,&valid);
               if(!valid) continue;
            }
            else
               time_fac = 1.0;
            //................................................................
            // Time dependencies - FCT
            if(m_bc_node->fct_name.length()>0)
            {
               m_fct = FCTGet(m_bc_node->fct_name);
               if(m_fct)
               {
                  time_fac = m_fct->GetValue(aktuelle_zeit,&is_valid);
                  //if(!valid) continue;
               }
               else
               {
                  cout << "Warning in CRFProcess::IncorporateBoundaryConditions - no FCT data" << endl;
               }
            }
            //................................................................
            // Conditions
            if(m_bc_node->conditional)
            {
               bc_value = time_fac*fac* GetNodeValue(m_bc_node->msh_node_number_subst,
                                                  //WW  bc_value = time_fac*fac* GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1); // YD-----TEST---
                  GetNodeValueIndex(pcs_primary_function_name[0])+1);
            }
            else
                                                  // time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
               bc_value = time_fac*fac*m_bc_node->node_value;
            //----------------------------------------------------------------
            // MSH
            if(rank>-1)
               bc_eqs_index = bc_msh_node;
            else
                                                  //WW#
               bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
            //..............................................................
            // NEWTON WW
            if(m_num->nls_method_name.find("NEWTON")!=string::npos
               || type==4||type==41  )            //Solution is in the manner of increment !
            {
               idx0 = GetNodeValueIndex(convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()));
               if(type==4||type==41)
               {
                                                  //30.08.2010. WW
                  if(m_bc_node->pcs_pv_name.find("PRESSURE")==string::npos)
                  {
                     idx1 = idx0+1;
                     bc_value -=  GetNodeValue(m_bc_node->geo_node_number, idx0)
                        + GetNodeValue(m_bc_node->geo_node_number, idx1);
                  }
               }
               else
                  bc_value = bc_value - GetNodeValue(m_bc_node->geo_node_number, idx0);
            }
            //----------------------------------------------------------------
            bc_eqs_index += shift;
            if((int)continuum_vector.size() > 1)
            {
               //YD/WW
               if(m_bc_node->pcs_pv_name.find(pcs_primary_function_name[continuum]) == string::npos)
                  continue;
            }
#ifdef NEW_EQS                        //WW
            eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#else
            MXRandbed(bc_eqs_index,bc_value,eqs_rhs);
#endif
         }
      }
                                                  // PCH
      else if(axis == 1 && m_bc_node->pcs_pv_name.find("VELOCITY1_Y")!=string::npos)
      {
         m_bc = bc_node[gindex];
         shift = m_bc_node->msh_node_number-m_bc_node->geo_node_number;
         //
         if(rank>-1)
         {
            bc_msh_node = bc_local_index_in_dom[i];
            int dim_space = 0;
            if(m_msh->NodesNumber_Linear==m_msh->NodesNumber_Quadratic)
               dim_space = 0;
            else
            {
               if(shift%m_msh->NodesNumber_Quadratic==0)
                  dim_space = shift/m_msh->NodesNumber_Quadratic;
               else
                  dim_space = m_msh->msh_max_dim;
            }
            shift = m_dom->shift[dim_space];
         }
         else
            bc_msh_node = m_bc_node->geo_node_number;
         //------------------------------------------------------------WW
         if(m_msh)                                //OK
            //			if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
            //				continue;
            time_fac = 1.0;
         if(bc_msh_node>=0)
         {
            //................................................................
            // Time dependencies - CURVE
            curve =  m_bc_node->CurveIndex;
            if(curve>0)
            {
               time_fac = GetCurveValue(curve,interp_method,aktuelle_zeit,&valid);
               if(!valid) continue;
            }
            else
               time_fac = 1.0;
            //................................................................
            // Time dependencies - FCT
            if(m_bc_node->fct_name.length()>0)
            {
               m_fct = FCTGet(m_bc_node->fct_name);
               if(m_fct)
               {
                  time_fac = m_fct->GetValue(aktuelle_zeit,&is_valid);
                  //if(!valid) continue;
               }
               else
               {
                  cout << "Warning in CRFProcess::IncorporateBoundaryConditions - no FCT data" << endl;
               }
            }
            //................................................................
            // Conditions
            if(m_bc_node->conditional)
            {
               bc_value = time_fac*fac* GetNodeValue(m_bc_node->msh_node_number_subst,
                                                  //WW  bc_value = time_fac*fac* GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1); // YD-----TEST---
                  GetNodeValueIndex(pcs_primary_function_name[0])+1);
            }
            else
                                                  // time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
               bc_value = time_fac*fac*m_bc_node->node_value;
            //----------------------------------------------------------------
            // MSH
            if(rank>-1)
               bc_eqs_index = bc_msh_node;
            else
                                                  //WW#
               bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
            //..............................................................
            // NEWTON WW
            if(m_num->nls_method_name.find("NEWTON")!=string::npos
               || type==4||type==41  )            //Solution is in the manner of increment !
            {
               idx0 = GetNodeValueIndex(convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()).c_str());
               if(type==4||type==41)
               {
                  idx1 = idx0+1;
                  bc_value -=  GetNodeValue(m_bc_node->geo_node_number, idx0)
                     + GetNodeValue(m_bc_node->geo_node_number, idx1);
               }
               else
                  bc_value = bc_value - GetNodeValue(m_bc_node->geo_node_number, idx0);
            }
            //----------------------------------------------------------------
            bc_eqs_index += shift;
            if((int)continuum_vector.size() > 1)
            {
               //YD/WW
               if(m_bc_node->pcs_pv_name.find(pcs_primary_function_name[continuum]) == string::npos)
                  continue;
            }
#ifdef NEW_EQS                        //WW
            eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#else
            MXRandbed(bc_eqs_index,bc_value,eqs_rhs);
#endif
         }
      }
                                                  // PCH
      else if(axis == 2 && m_bc_node->pcs_pv_name.find("VELOCITY1_Z")!=string::npos)
      {
         m_bc = bc_node[gindex];
         shift = m_bc_node->msh_node_number-m_bc_node->geo_node_number;
         //
         if(rank>-1)
         {
            bc_msh_node = bc_local_index_in_dom[i];
            int dim_space = 0;
            if(m_msh->NodesNumber_Linear==m_msh->NodesNumber_Quadratic)
               dim_space = 0;
            else
            {
               if(shift%m_msh->NodesNumber_Quadratic==0)
                  dim_space = shift/m_msh->NodesNumber_Quadratic;
               else
                  dim_space = m_msh->msh_max_dim;
            }
            shift = m_dom->shift[dim_space];
         }
         else
            bc_msh_node = m_bc_node->geo_node_number;
         //------------------------------------------------------------WW
         if(m_msh)                                //OK
            //			if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
            //				continue;
            time_fac = 1.0;
         if(bc_msh_node>=0)
         {
            //................................................................
            // Time dependencies - CURVE
            curve =  m_bc_node->CurveIndex;
            if(curve>0)
            {
               time_fac = GetCurveValue(curve,interp_method,aktuelle_zeit,&valid);
               if(!valid) continue;
            }
            else
               time_fac = 1.0;
            //................................................................
            // Time dependencies - FCT
            if(m_bc_node->fct_name.length()>0)
            {
               m_fct = FCTGet(m_bc_node->fct_name);
               if(m_fct)
               {
                  time_fac = m_fct->GetValue(aktuelle_zeit,&is_valid);
                  //if(!valid) continue;
               }
               else
               {
                  cout << "Warning in CRFProcess::IncorporateBoundaryConditions - no FCT data" << endl;
               }
            }
            //................................................................
            // Conditions
            if(m_bc_node->conditional)
            {
               bc_value = time_fac*fac* GetNodeValue(m_bc_node->msh_node_number_subst,
                                                  //WW  bc_value = time_fac*fac* GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1); // YD-----TEST---
                  GetNodeValueIndex(pcs_primary_function_name[0])+1);
            }
            else
                                                  // time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
               bc_value = time_fac*fac*m_bc_node->node_value;
            //----------------------------------------------------------------
            // MSH
            if(rank>-1)
               bc_eqs_index = bc_msh_node;
            else
                                                  //WW#
               bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
            //..............................................................
            // NEWTON WW
            if(m_num->nls_method_name.find("NEWTON")!=std::string::npos
               || type==4||type==41  )            //Solution is in the manner of increment !
            {
               idx0 = GetNodeValueIndex(convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()).c_str());
               if(type==4||type==41)
               {
                  idx1 = idx0+1;
                  bc_value -=  GetNodeValue(m_bc_node->geo_node_number, idx0)
                     + GetNodeValue(m_bc_node->geo_node_number, idx1);
               }
               else
                  bc_value = bc_value - GetNodeValue(m_bc_node->geo_node_number, idx0);
            }
            //----------------------------------------------------------------
            bc_eqs_index += shift;
            if((int)continuum_vector.size() > 1)
            {
               //YD/WW
               if(m_bc_node->pcs_pv_name.find(pcs_primary_function_name[continuum]) == std::string::npos)
                  continue;
            }
#ifdef NEW_EQS                        //WW
            eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#else
            MXRandbed(bc_eqs_index,bc_value,eqs_rhs);
#endif
         }
      }
   }

   //-----------------------------------------------------------------------
   /* irreg. Zeilen/Spalten regularisieren */
   /*
     else if (GetNodeState(NodeNumber[i]) == -2 || GetNodeState(NodeNumber[i]) == -4) { // irreg.Knoten
       if (GetRFControlGridAdapt())
         if (AdaptGetMethodIrrNodes() == 1) {
           MXSet(i, i, MKleinsteZahl);
           rechts[i] = 0.0;
         }
       }
     }
   */

}


/**************************************************************************
FEMLib-Method:
Task: PCS source terms into EQS
Programing:
04/2004 OK Implementation
08/2004 WW Extension for monolithic PCS and time curve
last modification:
02/2005 MB River Condition and CriticalDepth
05/2005 WW Dynamic problems
07/2005 WW Changes due to the geometry object applied
03/2006 WW Re-arrange
04/2006 OK CPL
05/2006 WW DDC
08/2006 YD FCT use
**************************************************************************/
void CRFProcess::IncorporateSourceTerms(const int rank)
{
   double value = 0, fac = 1.0, time_fac;
   int interp_method = 0;
   int curve, valid = 0;
   long msh_node, shift;
   long bc_eqs_index = -1;
   MshElemType::type EleType;                     //ii
   double q_face = 0.0;
   CElem* elem = NULL;
   CElem* face = NULL;
   ElementValue* gp_ele = NULL;
   CPARDomain *m_dom = NULL;
   double *eqs_rhs = NULL;
   double vel[3];
   bool is_valid;                                 //YD
   CFunction* m_fct = NULL;                       //YD
   long i;                                        //, group_vector_length;

   double Scaling = 1.0;
   if (type == 4)
      fac = Scaling;

   CNodeValue *cnodev = NULL;
   CSourceTerm *m_st = NULL;
   //
   long begin = 0;
   long end = 0;
   long gindex = 0;
   int dim_space = 0;                             //kg44 better define here and not in a loop!

   if (rank == -1)
   {
      begin = 0;
      end = (long) st_node_value.size();
#ifdef NEW_EQS                              //WW
      eqs_rhs = eqs_new->b;                       //27.11.2007 WW
#else
      eqs_rhs = eqs->b;
#endif
   }
   else
   {
      m_dom = dom_vector[rank];
#ifdef NEW_EQS
      if(type == 4 )
         eqs_rhs = m_dom->eqsH->b;
      else
         eqs_rhs = m_dom->eqs->b;
#else
      eqs_rhs = m_dom->eqs->b;
#endif
      if (rank == 0)
         begin = 0;
      else
         begin = rank_st_node_value_in_dom[rank - 1];
      end = rank_st_node_value_in_dom[rank];
   }

   for (i = begin; i < end; i++)
   {
      gindex = i;
      if (rank > -1)
         gindex = st_node_value_in_dom[i];

      cnodev = st_node_value[gindex];
      shift = cnodev->msh_node_number - cnodev->geo_node_number;
      if (rank > -1)
      {
         msh_node = st_local_index_in_dom[i];
         dim_space = 0;
         if (m_msh->NodesNumber_Linear == m_msh->NodesNumber_Quadratic)
            dim_space = 0;
         else
         {
            if (shift % m_msh->NodesNumber_Quadratic == 0)
               dim_space = shift / m_msh->NodesNumber_Quadratic;
            else
               dim_space = m_msh->msh_max_dim;
         }
         shift = m_dom->shift[dim_space];
      }
      else
      {
         msh_node = cnodev->msh_node_number;
         msh_node -= shift;
      }
      value = cnodev->node_value;
      //--------------------------------------------------------------------
      // Tests
      if (msh_node < 0)
         continue;
      m_st = NULL;
      if (st_node.size() > 0 && (long) st_node.size() > i)
      {
         m_st = st_node[gindex];
         //--------------------------------------------------------------------
         // CPL
         //if(m_st->_pcs_type_name_cond.size()>0) continue; // this is a CPL source term, JOD removed
         //--------------------------------------------------------------------
         // system dependent YD
         if (cnodev->getProcessDistributionType() == FiniteElement::SYSTEM_DEPENDENT)
         {
            long no_st_ele = (long) m_st->element_st_vector.size();
            for (long i_st = 0; i_st < no_st_ele; i_st++)
            {
               long ele_index = m_st->element_st_vector[i_st];
               elem = m_msh->ele_vector[ele_index];
               if (elem->GetMark())
               {
                  fem->ConfigElement(elem);
                  fem->Cal_Velocity();
               }
               gp_ele = ele_gp_value[ele_index];
               gp_ele->GetEleVelocity(vel);
               EleType = elem->GetElementType();
               if (EleType == MshElemType::LINE)  //Line
                  cnodev->node_value += vel[0];
                                                  //Traingle & Qua
               if (EleType == MshElemType::TRIANGLE || EleType == MshElemType::QUAD)
               {
                  for (long i_face = 0; i_face
                     < (long) m_msh->face_vector.size(); i_face++)
                  {
                     face = m_msh->face_vector[i_face];
                     if ((size_t) m_st->element_st_vector[i_st]
                        == face->GetOwner()->GetIndex())
                        q_face = PointProduction(vel,
                           m_msh->face_normal[i_face])
                           * face->GetVolume();   //
                     //for(i_node)
                  }
                  cnodev->node_value = +q_face / 2;
               }
               // cout<<"  value  "<<m_st_group->group_vector[i]->node_value<<endl;
            }
         }
         //--------------------------------------------------------------------
         // MB
         //if(m_st->conditional && !m_st->river)
         //{

         GetNODValue(value, cnodev, m_st);

      }                                           // st_node.size()>0&&(long)st_node.size()>i
      //----------------------------------------------------------------------------------------
      //--------------------------------------------------------------------
      // Please do not move the this section
      curve = cnodev->CurveIndex;
      if (curve > 0)
      {
         time_fac = GetCurveValue(curve, interp_method, aktuelle_zeit,
            &valid);
         if (!valid)
         {
            cout << "\n!!! Time dependent curve is not found. Results are not guaranteed "
               << endl;
            cout << " in void CRFProcess::IncorporateSourceTerms(const double Scaling)"
               << endl;
            time_fac = 1.0;
         }
      } else
      time_fac = 1.0;

      // Time dependencies - FCT    //YD
      if (m_st)                                   //WW
      {
                                                  //WW/YD //OK
         if (m_msh && m_msh->geo_name.find("LOCAL") != string::npos)
         {
            if (m_st->getFunctionName().length() > 0)
            {
               m_fct = FCTGet(pcs_number);
               if (m_fct)
                  time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid,
                     m_st->getFunctionMethod());  //fct_method. WW
               else
                  cout << "Warning in CRFProcess::IncorporateSourceTerms - no FCT data"
                     << endl;
            }
         }
         else
         {
            if (m_st->getFunctionName().length() > 0)
            {
               m_fct = FCTGet(m_st->getFunctionName());
               if (m_fct)
               {
                  time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
               }
               else
               {
                  cout << "Warning in CRFProcess::IncorporateSourceTerms - no FCT data"
                     << endl;
               }
            }
         }
      }
      //----------------------------------------------------------------------------------------
      value *= time_fac * fac;
      //------------------------------------------------------------------
      // EQS->RHS
      if (rank > -1)
         bc_eqs_index = msh_node + shift;
      else
         bc_eqs_index = m_msh->nod_vector[msh_node]->GetEquationIndex() + shift;
      eqs_rhs[bc_eqs_index] += value;
   }
   //====================================================================

   // if coupling to GEMS exist----------------------------------------------------
   // HS, added 11.2008
   // KG44 03/03/2010 modified to hopefully soon work with parallel solvers
   long gem_node_index = -1, glocalindex = -1;
   if (flag_couple_GEMS == 1 && aktueller_zeitschritt > 1)
   {
      begin = 0;
      if (rank == -1)                             // serial version
      {
         end = (long ) Water_ST_vec.size();
      } else                                      // parallel version
      {
         end = 0;
         if (rank_stgem_node_value_in_dom.size() > 0)
            end = rank_stgem_node_value_in_dom[0];
      }
      // only when switch is on and not in the first time step
      // loop over the Water_ST_vec vector,
      // add the excess water to the right-hand-side of the equation
      for (i = begin; i < end; i++)
      {
         if (rank > -1)                           // parallel version: stgem_node_value_in_dom and stgem_local_index_in_dom contain only values for the corresponding domain  == rank
         {
            //				cout << "rank " << rank ;
            gindex = stgem_node_value_in_dom[i];  //contains indexes to water-st_vec
            //				cout << " gindex " << gindex << " i " << i << endl ;
                                                  // contains index to node
            glocalindex = stgem_local_index_in_dom[i];
            //				cout << " gem_node_index " << gem_node_index << endl;
            eqs_rhs[glocalindex] += Water_ST_vec[gindex].water_st_value;

         } else                                   // serial version
         {
            gem_node_index = Water_ST_vec[i].index_node;
            eqs_rhs[gem_node_index] += Water_ST_vec[i].water_st_value;
         }
      }
      // after finished adding to RHS, clear the vector
      Water_ST_vec.clear();
      if (rank > -1)
      {
         stgem_node_value_in_dom.clear();
         stgem_local_index_in_dom.clear();
         rank_stgem_node_value_in_dom.clear();
      }
   }
}


#ifndef NEW_EQS                                   //WW
/**************************************************************************
FEMLib-Method:
Task:
Programing:
02/2004 OK Implementation
11/2004 OK NUM
07/2006 WW Parallel BiCGStab
last modification:
**************************************************************************/
int CRFProcess::ExecuteLinearSolver(void)
{
   long iter_count = 0;
   long iter_sum = 0;
   //WW  int found = 0;
   //-----------------------------------------------------------------------
   // Set EQS
   //cout << "Before SetLinearSolver(eqs) myrank = "<< myrank<< '\n';
   SetLinearSolver(eqs);
   //--------------------------------------------------------------------
   // NUM
   //WW  found = 1;
   cg_maxiter        = m_num->ls_max_iterations;  //OK lsp->maxiter;
   cg_eps            = m_num->ls_error_tolerance; //OK lsp->eps;
   //cg_repeat = lsp->repeat;
   vorkond           = m_num->ls_precond;         //OK lsp->precond;
   linear_error_type = m_num->ls_error_method;    //OK lsp->criterium;

   //  cout << "Before eqs->LinearSolver(eqs) myrank = "<< myrank<< '\n';

#ifdef USE_MPI
   //WW
   long dim_eqs = 0;
   if(type==41||type==4)                          // DOF >1
   {
      dom_vector[myrank]->quadratic = true;
      if(type==4)
         dim_eqs = pcs_number_of_primary_nvals*m_msh->GetNodesNumber(true);
      else if(type==41)
         dim_eqs = pcs_number_of_primary_nvals*m_msh->GetNodesNumber(true)
               +m_msh->GetNodesNumber(false);
   }
   else
   {
      dom_vector[myrank]->quadratic = false;
      dim_eqs = m_msh->GetNodesNumber(false);
   }
   iter_count =  SpBICGSTAB_Parallel(dom_vector[myrank], eqs->x, dim_eqs);
#else
   iter_count = eqs->LinearSolver(eqs->b,eqs->x,eqs->dim);
#endif

   eqs->master_iter = iter_count;
   if(iter_count>=cg_maxiter)
   {
      cout << "Warning in CRFProcess::ExecuteLinearSolver() - Maximum iteration number reached" << endl;
      return -1;
   }
   iter_sum += iter_count;
   //-----------------------------------------------------------------------
   // Clean results ?
   /*
     for (i=0;i<eqs->dim;i++)
       if (fabs(eqs->x[i])<MKleinsteZahl)
         eqs->x[i] = 0.0;
   */
   //-----------------------------------------------------------------------
   return iter_sum;
}
#endif
/**************************************************************************
FEMLib-Method:
Task:
Programing:
06/2005 PCH Overriding
last modification:
**************************************************************************/
#ifndef NEW_EQS                                   //WW 07.11.2008
int CRFProcess::ExecuteLinearSolver(LINEAR_SOLVER *eqs)
{
   long iter_count;
   long iter_sum = 0;
   int found = 0;
   //-----------------------------------------------------------------------
   // Set EQS
   SetLinearSolver(eqs);
   //--------------------------------------------------------------------
   // NUM
   found = 1;
   cg_maxiter        = m_num->ls_max_iterations;  //OK lsp->maxiter;
   cg_eps            = m_num->ls_error_tolerance; //OK lsp->eps;
   //cg_repeat = lsp->repeat;
   vorkond           = m_num->ls_precond;         //OK lsp->precond;
   linear_error_type = m_num->ls_error_method;    //OK lsp->criterium;

   iter_count = eqs->LinearSolver(eqs->b,eqs->x,eqs->dim);
   eqs->master_iter = iter_count;
   if(iter_count>=cg_maxiter)
   {
      cout << "Warning in CRFProcess::ExecuteLinearSolver() - Maximum iteration number reached" << endl;
      return -1;
   }
   iter_sum += iter_count;
   //-----------------------------------------------------------------------
   return iter_sum;
}
#endif
//WW
int CRFProcess::GetNODValueIndex(const string &name,int timelevel)
{
   for(int j=0;j<number_of_nvals;j++)
   {
      if((name.compare(pcs_nval_data[j].name)==0) && \
         (pcs_nval_data[j].timelevel==timelevel))
         return pcs_nval_data[j].nval_index;
   }
   cout << "Error in PCSGetNODValueIndex: " << name << endl;
   return -1;
}


///////////////////////////////////////////////////////////////////////////
// Specials
///////////////////////////////////////////////////////////////////////////

/*-------------------------------------------------------------------------
ROCKFLOW - Function: PCSRestart
Task: Insert process to list
Programming:
06/2003 OK Implementation
11/2004 OK file_name_base
last modified:
-------------------------------------------------------------------------*/
void PCSRestart()
{
   /*OK411
     int j;
     CRFProcess *m_pcs = NULL;
     int nidx0,nidx1;
     int i;
     int no_processes =(int)pcs_vector.size();
     if(no_processes==0)
       return; //OK41
     int ok = 0;
     //----------------------------------------------------------------------
     string file_name_base = pcs_vector[0]->file_name_base;
   //OK  ok = ReadRFRRestartData(file_name_base);
   if(ok==0){
   cout << "RFR: no restart data" << endl;
   return;
   }
   //----------------------------------------------------------------------
   for(i=0;i<no_processes;i++){
   m_pcs = pcs_vector[i];
   for(j=0;j<m_pcs->GetPrimaryVNumber();j++) {
   // timelevel=0;
   nidx0 = m_pcs->GetNodeValueIndex(m_pcs->GetPrimaryVName(j));
   // timelevel= 1;
   nidx1 = nidx0+1;
   CopyNodeVals(nidx1,nidx0);
   }
   }
   */
}


/**************************************************************************
FEMLib-Method:
Task: Relocate Deformation process
Programing:
09/2004 WW Implementation
10/2010 TF changes due to conversion from std::string to enum for process type
**************************************************************************/
void RelocateDeformationProcess(CRFProcess *m_pcs)
{
   //   string pcs_name_dm = m_pcs->_pcs_type_name;
   ProcessType pcs_name_dm (m_pcs->getProcessType());

   string num_type_name_dm;
   // Numerics
   if(m_pcs->num_type_name.compare("STRONG_DISCONTINUITY")==0)
   {
      num_type_name_dm = m_pcs->num_type_name;
      enhanced_strain_dm=1;
   }

   delete m_pcs;
   m_pcs = dynamic_cast<CRFProcess *> (new CRFProcessDeformation());
   m_pcs->setProcessType (pcs_name_dm);

   if(enhanced_strain_dm==1) m_pcs->num_type_name = num_type_name_dm;
   pcs_deformation = 1;
}


/*************************************************************************
ROCKFLOW - Function: CRFProcess::PCSMoveNOD
Task:
Programming:
08/2004 MB/OK Implementation
last modified:
**************************************************************************/
void CRFProcess::PCSMoveNOD(void)
{
   switch(this->type)
   {
      case 1:
         MSHMoveNODUcFlow(this);
         break;
      default:
         DisplayMsgLn("PCSMoveNOD: no valid process");
         abort();
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
09/2004 OK Implementation
10/2004 OK 2nd version
**************************************************************************/
std::string PCSProblemType()
{
   std::string pcs_problem_type;
   size_t no_processes(pcs_vector.size());

   for (size_t i = 0; i < no_processes; i++)
   {
      switch (pcs_vector[i]->getProcessType())
      {
         case LIQUID_FLOW:
            pcs_problem_type = "LIQUID_FLOW";
            break;
         case OVERLAND_FLOW:
            pcs_problem_type = "OVERLAND_FLOW";
            break;
         case GROUNDWATER_FLOW:
            pcs_problem_type = "GROUNDWATER_FLOW";
            break;
         case TWO_PHASE_FLOW:
            pcs_problem_type = "TWO_PHASE_FLOW";
            break;
         case RICHARDS_FLOW:                      //MX test 04.2005
            pcs_problem_type = "RICHARDS_FLOW";
            break;
         case DEFORMATION:
            if (pcs_problem_type.empty())
               pcs_problem_type = "DEFORMATION";
            else
               pcs_problem_type += "+DEFORMATION";
            break;
         case DEFORMATION_FLOW:
            if (pcs_problem_type.empty())
               pcs_problem_type = "DEFORMATION";
            else
               pcs_problem_type += "+DEFORMATION";
            break;
         case HEAT_TRANSPORT:
            if (pcs_problem_type.empty())
               pcs_problem_type = "HEAT_TRANSPORT";
            else
               pcs_problem_type += "+HEAT_TRANSPORT";
            break;
         case MASS_TRANSPORT:
            if (pcs_problem_type.empty())
               pcs_problem_type = "MASS_TRANSPORT";
            else
               pcs_problem_type += "+MASS_TRANSPORT";
            break;
         case FLUID_MOMENTUM:
            if (pcs_problem_type.empty())
               pcs_problem_type = "FLUID_MOMENTUM";
            else
               pcs_problem_type += "+FLUID_MOMENTUM";
            break;
         case RANDOM_WALK:
            if (pcs_problem_type.empty())
               pcs_problem_type = "RANDOM_WALK";
            else
               pcs_problem_type += "+RANDOM_WALK";
            break;
         default:
            pcs_problem_type = "";
      }
   }

   //	//----------------------------------------------------------------------
   //	CRFProcess* m_pcs = NULL;
   //	// H process
   //	for (i = 0; i < no_processes; i++) {
   //		m_pcs = pcs_vector[i];
   //		switch (m_pcs->pcs_type_name[0]) {
   //		case 'L':
   //			pcs_problem_type = "LIQUID_FLOW";
   //			break;
   //			// case 'U':
   //			//  pcs_problem_type = "UNCONFINED_FLOW";
   //			//  break;
   //		case 'O':
   //			pcs_problem_type = "OVERLAND_FLOW";
   //			break;
   //		case 'G':
   //			pcs_problem_type = "GROUNDWATER_FLOW";
   //			break;
   //		case 'T':
   //			pcs_problem_type = "TWO_PHASE_FLOW";
   //			break;
   //		case 'C':
   //			pcs_problem_type = "COMPONENTAL_FLOW";
   //			break;
   //		case 'R': //MX test 04.2005
   //			pcs_problem_type = "RICHARDS_FLOW";
   //			break;
   //		}
   //	}
   //	//----------------------------------------------------------------------
   //	// M process
   //	for (i = 0; i < no_processes; i++) {
   //		m_pcs = pcs_vector[i];
   //		switch (m_pcs->pcs_type_name[0]) {
   //		case 'D':
   //			if (pcs_problem_type.empty())
   //				pcs_problem_type = "DEFORMATION";
   //			else
   //				pcs_problem_type += "+DEFORMATION";
   //			break;
   //		}
   //	}
   //	//----------------------------------------------------------------------
   //	// T process
   //	for (i = 0; i < no_processes; i++) {
   //		m_pcs = pcs_vector[i];
   //		switch (m_pcs->pcs_type_name[0]) {
   //		case 'H':
   //			if (pcs_problem_type.empty())
   //				pcs_problem_type = "HEAT_TRANSPORT";
   //			else
   //				pcs_problem_type += "+HEAT_TRANSPORT";
   //			break;
   //		}
   //	}
   //	//----------------------------------------------------------------------
   //	// CB process
   //	for (i = 0; i < no_processes; i++) {
   //		m_pcs = pcs_vector[i];
   //		switch (m_pcs->pcs_type_name[0]) {
   //		case 'M':
   //			if (pcs_problem_type.empty())
   //				pcs_problem_type = "MASS_TRANSPORT";
   //			else
   //				pcs_problem_type += "+MASS_TRANSPORT";
   //			break;
   //		}
   //	}
   //	//----------------------------------------------------------------------
   //	//----------------------------------------------------------------------
   //	// FM process
   //	for (i = 0; i < no_processes; i++) {
   //		m_pcs = pcs_vector[i];
   //		switch (m_pcs->pcs_type_name[0]) {
   //		case 'F':
   //			if (pcs_problem_type.empty())
   //				pcs_problem_type = "FLUID_MOMENTUM";
   //			else
   //				pcs_problem_type += "+FLUID_MOMENTUM";
   //			break;
   //		}
   //	}
   //	//----------------------------------------------------------------------
   //	for (i = 0; i < no_processes; i++) {
   //		m_pcs = pcs_vector[i];
   //		switch (m_pcs->pcs_type_name[7]) { // _pcs_type_name[7] should be 'W' because 'R' is reserved for Richard Flow.
   //		case 'W':
   //			if (pcs_problem_type.empty())
   //				pcs_problem_type = "RANDOM_WALK";
   //			else
   //				pcs_problem_type += "+RANDOM_WALK";
   //			break;
   //		}
   //	}

   return pcs_problem_type;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
11/2004 OK Implementation
**************************************************************************/
void CRFProcess::CalcELEMassFluxes(void)
{
   /*OK411
     int i;
     double e_value = -1.0;
     string e_value_name;
     double geo_factor, density;
     double velocity = 0.0;
     int e_idx;
     int phase = 0;
     long e;
     CMediumProperties* m_mmp = NULL;
     CFluidProperties* m_mfp = NULL;
   m_mfp = mfp_vector[phase]; //OK ToDo
   //======================================================================
   for(e=0;e<ElementListLength;e++){
   m_mmp = mmp_vector[ElGetElementGroupNumber(e)];
   geo_factor = m_mmp->geo_area;
   density = m_mfp->Density();
   for(i=0;i<pcs_number_of_evals;i++){
   e_value_name = pcs_eval_data[i].name;
   e_idx = PCSGetELEValueIndex(pcs_eval_data[i].name);
   if(e_value_name.find("MASS_FLUX1_X")!=string::npos){
   velocity = ElGetElementVal(e,PCSGetELEValueIndex("VELOCITY1_X"));
   e_value = geo_factor * density * velocity;
   ElSetElementVal(e,e_idx,e_value);
   }
   if(e_value_name.find("MASS_FLUX1_Y")!=string::npos){
   velocity = ElGetElementVal(e,PCSGetELEValueIndex("VELOCITY1_Y"));
   e_value = geo_factor * density * velocity;
   ElSetElementVal(e,e_idx,e_value);
   }
   if(e_value_name.find("MASS_FLUX1_Z")!=string::npos){
   velocity = ElGetElementVal(e,PCSGetELEValueIndex("VELOCITY1_Z"));
   e_value = geo_factor * density * velocity;
   ElSetElementVal(e,e_idx,e_value);
   }
   }
   }
   //======================================================================
   */
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
02/2005 OK Implementation
02/2006 YD Dual Richards
02/2007 WW General function for all unsaturated flow
**************************************************************************/
void CRFProcess::CalcSecondaryVariables(bool initial)
{
   //  char pcsT;
   //  pcsT = _pcs_type_name[0];
   //  if(type==1212) pcsT = 'V'; //WW
   //  switch(pcsT){
   //    case 'L':
   //      break;
   //    case 'U':
   //      break;
   //    case 'G':
   //      break;
   //    case 'T':
   //      break;
   //    case 'C':
   //      break;
   //    case 'R': // Richards flow
   //	  if(_pcs_type_name[1] == 'I')	// PCH To make a distinction with RANDOM WALK.
   //        CalcSecondaryVariablesUnsaturatedFlow(initial); // WW
   //      break;
   //    case 'D':
   //      break;
   //    case 'V':
   //      CalcSecondaryVariablesUnsaturatedFlow(initial); //WW
   //      break;
   //	case 'P':
   //      CalcSecondaryVariablesPSGLOBAL(); //WW
   //      break;
   //  }

   switch (getProcessType())
   {
      case LIQUID_FLOW:
         break;
      case GROUNDWATER_FLOW:
         break;
      case TWO_PHASE_FLOW:
         break;
      case RICHARDS_FLOW:                         // Richards flow
                                                  // WW
         CalcSecondaryVariablesUnsaturatedFlow(initial);
         break;
      case DEFORMATION || DEFORMATION_FLOW || DEFORMATION_DYNAMIC:
         break;
      case PS_GLOBAL:
         CalcSecondaryVariablesPSGLOBAL();        //WW
         break;
      default:
         if (type == 1212)
                                                  //WW
            CalcSecondaryVariablesUnsaturatedFlow(initial);
         break;
   }
}


//////////////////////////////////////////////////////////////////////////
// ReMove site
//////////////////////////////////////////////////////////////////////////

/*************************************************************************
ROCKFLOW - Function: GetCompNamehelp
Task: Namepatch, until primary function names are finally sorted out //SB:todo
Programming:	08/2003 SB Implementation
last modified:
superseded by GetPFNamebyCPName() but left here, as not all files are already in the new concept
**************************************************************************/
/* SB: namepatch
   Repariert kurzfristig die Ausgabe
   input: datafield_n[j].name
   wenn dar Name "CONCENTRATIONx" ist, wird er durch den enstprechenden Komponentennamen ersetzt, sonst bleibts */

char *GetCompNamehelp(char *inname)
{
   int comp;                                      //WW, phase;
   char * outname, help[MAX_ZEILE];
   CRFProcess *m_pcs=NULL;
   outname = inname;
   //WW phase = 0;
   for(comp=0;comp<GetRFProcessNumComponents();comp++)
   {
      sprintf(help,"%s%d","CONCENTRATION",comp+1);
      /*  help has to be a part of inname (strstr) and also have the same length (strcmp) */
      if(strstr(inname, help) && (strcmp(inname,help) == 0))
      {
         m_pcs = m_pcs->GetProcessByFunctionName(help);
         if(m_pcs == NULL) break;
         //		outname = GetTracerCompName(phase,m_pcs->GetProcessComponentNumber()-1);
         //		outname = cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname;
         outname = (char *) cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname.data();
         return outname;
      };
   };
   return outname;
}                                                 //SB:namepatch


/*************************************************************************
ROCKFLOW - Function: GetCPNamebyPFName(string )
Task: Replaces CP Name by Primary function name for output input
Input:	component property name
Output: primary function name
Programming:	10/2004 SB Implementation
10/2010 TF changed access to process type
**************************************************************************/
string GetPFNamebyCPName(string inname)
{
   int i, j, k;
   int pcs_vector_size = (int) pcs_vector.size();
   string outname;
   char  help[MAX_ZEILE];
   CRFProcess *m_pcs=NULL;
   outname = "dummy";
   if(pcs_vector_size>0)
      for(i=0;i<pcs_vector_size;i++)
   {
      m_pcs = pcs_vector[i];
      //	if(m_pcs->_pcs_type_name.compare("MASS_TRANSPORT") == 0){ // if this is mass transport // TF
      if(m_pcs->getProcessType() == MASS_TRANSPORT)
      {
         j = m_pcs->GetProcessComponentNumber();
         k = cp_vec[j]->transport_phase;
         outname = cp_vec[m_pcs->GetProcessComponentNumber()]->compname;
         if(outname == inname)                    //right process found
         {
            sprintf(help,"%s%d","CONCENTRATION",j);
            outname = help;
            return outname;
         }
      }
   }
   // for(i=0;comp<GetRFProcessNumComponents();i++) {
   //	sprintf(help,"%s%d","CONCENTRATION",i);
   /*  help has to be a part of inname (strstr) and also have the same length (strcmp) */
   //	if(strstr(inname, help) && (strcmp(inname,help) == 0)){
   //		m_pcs = m_pcs->GetProcessByFunctionName(help);
   //		if(m_pcs == NULL) break;
   //		outname = cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname;
   //		outname = (char *) cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname.data();
   //		if(outname.compare(inname) == 0)
   //			return outname;
   //	};
   // };
   // Inname is not from a mass transport process, therefore return inname
   return inname;
}                                                 //SB:namepatch


//========================================================================
//OK former model functions
int GetRFControlGridAdapt(void)
{
   //OK  return (get_rfcp_adaptive_mesh_refinement_flag(rfcp));
   if(show_onces_adp)
      cout << "GetRFControlGridAdapt - to be removed" << endl;
   show_onces_adp = false;
   return 0;
}


int GetRFControlModel(void)
{
   if(show_onces_mod)
      cout << "GetRFControlModel - to be removed" << endl;
   show_onces_mod = false;
   return -1;
}


int GetRFProcessChemicalModel(void)
{
   cout << "GetRFProcessChemicalModel - to be removed" << endl;
   return 0;
}


int GetRFProcessFlowModel(void)
{
   if(show_onces_mod_flow)
      cout << "GetRFProcessFlowModel - to be removed" << endl;
   show_onces_mod_flow = false;
   return 0;
}


int GetRFProcessHeatReactModel(void)
{
   cout << "GetRFProcessHeatReactModel - to be removed" << endl;
   return 0;
}


int GetRFProcessNumPhases(void)
{
   //DisplayMsgLn("GetRFProcessNumPhases - to be removed");
   int no_phases = (int)mfp_vector.size();
   return no_phases;
}


int GetRFProcessProcessing(char* rfpp_type)
{
   bool pcs_flow = false;
   bool pcs_deformation = false;
   CRFProcess* m_pcs = NULL;
   size_t no_processes = pcs_vector.size();
   for (size_t i = 0; i < no_processes; i++)
   {
      m_pcs = pcs_vector[i];
      //		if (m_pcs->_pcs_type_name.find("DEFORMATION") != string::npos)
      if (isDeformationProcess (m_pcs->getProcessType()))
         pcs_deformation = true;
      //		if (m_pcs->_pcs_type_name.find("FLOW") != string::npos)
      if (isFlowProcess (m_pcs->getProcessType()))
         pcs_flow = true;
   }

   if (strcmp(rfpp_type, "SD") == 0)
   {
      if (pcs_flow && pcs_deformation)
         return 1;
   } else
   cout << "GetRFProcessProcessing - to be removed" << endl;
   return 0;
}


int GetRFProcessProcessingAndActivation(const char*)
{
   cout << "GetRFProcessProcessingAndActivation - to be removed" << endl;
   return 0;
}


long GetRFProcessNumComponents(void)
{
   //DisplayMsgLn("GetRFProcessNumComponents - to be removed");
   int no_components = (int)cp_vec.size();
   return no_components;
}


int GetRFControlModex(void)
{
   cout << "GetRFControlModex - to be removed" << endl;
   return 0;
}


int GetRFProcessDensityFlow(void)
{
   if(show_onces_density)
      cout << "GetRFProcessDensityFlow - to be removed" << endl;
   show_onces_density = false;
   return 0;
}


int GetRFProcessNumContinua(void)
{
   cout << "GetRFProcessNumContinua - to be removed" << endl;
   return 0;
}


int GetRFProcessNumElectricFields(void)
{
   cout << "GetRFProcessNumElectricFields - to be removed" << endl;
   return 0;
}


int GetRFProcessNumTemperatures(void)
{
   cout << "GetRFProcessNumTemperatures - to be removed" << endl;
   return -1;
}


int GetRFProcessSimulation(void)
{
   cout << "GetRFProcessSimulation - to be removed" << endl;
   return -1;
}


/**************************************************************************
   ROCKFLOW - Funktion: ModelsAddNodeValInfoStructure

   Aufgabe:
   Fuellt die Knotendaten-Infostruktur mit den zugehoerigen Modelldaten.
   Wird vom Modell der Reihe nach fuer jede Knotengroesse aufgerufen.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E:char *name         :Name der Knotengroesse fuer Ergebnisdatei
   E:char *einheit      :Name der phys. Einheit fuer Ergebnisdatei
   E:int speichern      :Werte sollen gespeichert werden (0/1)
E:int laden          :Werte sollen geladen werden falls vorhanden (0/1)
E:int restart        :Werte sollen bei Restart geladen werden (0/1)
E:int adapt_interpol :Werte sollen beim verfeinern auf Kinder interpoliert (0/1)
E:double vorgabe     :Vorgabe falls keine Restartdaten oder Anfangsbedingungen vorhanden sind

Ergebnis:
Knotenindex der gerade vergeben wurde

Programmaenderungen:
09/2000   CT    Erste Version

**************************************************************************/
int ModelsAddNodeValInfoStructure(char *name, char *einheit, int speichern, int laden, int restart, int adapt_interpol, double vorgabe)
{
   anz_nval++;
   nval_data = (NvalInfo *) Realloc(nval_data, anz_nval * sizeof(NvalInfo));

   nval_data[anz_nval - 1].name = NULL;
   nval_data[anz_nval - 1].einheit = NULL;

   if (name)
   {
      nval_data[anz_nval - 1].name = (char *) Malloc(((int)strlen(name) + 1) * sizeof(char));
      strcpy(nval_data[anz_nval - 1].name, name);
   }
   if (einheit)
   {
      nval_data[anz_nval - 1].einheit = (char *) Malloc(((int)strlen(einheit) + 1) * sizeof(char));
      strcpy(nval_data[anz_nval - 1].einheit, einheit);
   }

   nval_data[anz_nval - 1].speichern = speichern;
   nval_data[anz_nval - 1].laden = laden;
   nval_data[anz_nval - 1].restart = restart;
   nval_data[anz_nval - 1].adapt_interpol = adapt_interpol;
   nval_data[anz_nval - 1].vorgabe = vorgabe;

   return anz_nval - 1;
}


/**************************************************************************
   ROCKFLOW - Funktion: ModelsAddElementValInfoStructure

   Aufgabe:
   Fuellt die Elementdaten-Infostruktur mit den zugehoerigen Modelldaten.
   Wird vom Modell der Reihe nach fuer jede Elementgroesse aufgerufen.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E:char *name         :Name der Elementgroesse fuer Ergebnisdatei
   E:char *einheit      :Name der phys. Einheit fuer Ergebnisdatei
   E:int speichern      :Werte sollen gespeichert werden (0/1)
E:int laden          :Werte sollen geladen werden falls vorhanden (0/1)
E:int restart        :Werte sollen bei Restart geladen werden (0/1)
E:int adapt_interpol :Werte sollen beim verfeinern auf Kinder interpoliert (0/1)
E:double vorgabe     :Vorgabe falls keine Restartdaten oder Anfangsbedingungen vorhanden sind

Ergebnis:
Elementindex der gerade vergeben wurde

Programmaenderungen:
09/2000   CT    Erste Version

**************************************************************************/
int ModelsAddElementValInfoStructure(char *name, char *einheit, int speichern, int laden, int restart, int adapt_interpol, double vorgabe)
{
   anz_eval++;
   eval_data = (EvalInfo *) Realloc(eval_data, anz_eval * sizeof(EvalInfo));

   eval_data[anz_eval - 1].name = NULL;
   eval_data[anz_eval - 1].einheit = NULL;

   if (name)
   {
      eval_data[anz_eval - 1].name = (char *) Malloc(((int)strlen(name) + 1) * sizeof(char));
      strcpy(eval_data[anz_eval - 1].name, name);
   }
   if (einheit)
   {
      eval_data[anz_eval - 1].einheit = (char *) Malloc(((int)strlen(einheit) + 1) * sizeof(char));
      strcpy(eval_data[anz_eval - 1].einheit, einheit);
   }

   eval_data[anz_eval - 1].speichern = speichern;
   eval_data[anz_eval - 1].laden = laden;
   eval_data[anz_eval - 1].restart = restart;
   eval_data[anz_eval - 1].adapt_interpol = adapt_interpol;
   eval_data[anz_eval - 1].vorgabe = vorgabe;

   return anz_eval - 1;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
last modified:
**************************************************************************/
void PCSDelete()
{
   for(int i=0;i<(int)pcs_vector.size();i++)
   {
      delete pcs_vector[i];
   }
   pcs_vector.clear();
   pcs_no_components = 0;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2005 OK Implementation
last modified:
**************************************************************************/
void CRFProcess::SetNodeValue(long n,int nidx,double value)
{
#ifdef gDEBUG
   if(nidx<0)
   {
      cout<<" Fatal error in  CRFProcess::SetNodeValue() "<<endl;
      abort();
   }
#endif
   nod_val_vector[n][nidx] = value;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
09/2005 PCH Implementation
last modified:
**************************************************************************/
void CRFProcess::SetElementValue(long n,int nidx,double value)
{
#ifdef gDEBUG
   if(nidx<0)
   {
      cout<<" Fatal error in  CRFProcess::SetElementValue() "<<endl;
      abort();
   }
#endif
   ele_val_vector[n][nidx] = value;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2005 OK Implementation
last modified:
**************************************************************************/
double CRFProcess::GetNodeValue(long n,int nidx)
{
   double value;
#ifdef gDEBUG
   if(nidx<0)
   {
      cout<<" Fatal error in  CRFProcess::GetNodeValue() "<<endl;
      abort();
   }
#endif
   value = nod_val_vector[n][nidx];
   return value;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
09/2005 PCH Implementation
last modified:
**************************************************************************/
double CRFProcess::GetElementValue(long n,int nidx)
{
   double value;
#ifdef gDEBUG
   if(nidx<0)
   {
      cout<<" Fatal error in CRFProcess::GetElementValue() "<<endl;
      abort();
   }
#endif
   value = ele_val_vector[n][nidx];
   return value;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2005 OK Implementation
01/2006 OK Test
**************************************************************************/
int CRFProcess::GetNodeValueIndex(const std::string &var_name)
{
   int nidx = -2;

   for (size_t i = 0; i < nod_val_name_vector.size(); i++)
   {
      if (nod_val_name_vector[i].compare(var_name) == 0)
      {
         nidx = i;
#ifdef gDEBUG
         if(nidx<0)
         {
            cout<<" Fatal error in CRFProcess::GetNodeValueIndex() "<<endl;
            abort();
         }
#endif
         return nidx;
      }
   }

   // Suppress the following error message when Fluid Momentum process is on.
   CRFProcess* m_pcs = PCSGet(FLUID_MOMENTUM);
   if (m_pcs)
      ;                                           // Don't print any error message.
   else
      cout << "Error in CRFProcess::GetNodeValueIndex - " << convertProcessTypeToString (getProcessType())
         << ", " << var_name << ", NIDX = " << nidx << endl;
   //abort();
   return nidx;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
09/2005 PCH Implementation
last modified:
**************************************************************************/
int CRFProcess::GetElementValueIndex(const string &var_name)
{
   int nidx = -1;
   for (size_t i = 0; i < this->ele_val_name_vector.size(); i++)
   {
      if (this->ele_val_name_vector[i].compare(var_name) == 0)
      {
         nidx = i;
         return nidx;
      }
   }
   return nidx;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2005 OK Implementation
05/2005 OK pcs_pv_name,
12/2005 OK RESTART
07/2006 OK/MX MSH
**************************************************************************/
void CRFProcess::SetIC()
{
   CInitialCondition* m_ic = NULL;
   // HS, for MASS_TRANSPORT PCS,
   // it is not necessary to use PrimaryVarible as second check.
   // nidx will give the proper IC pointer.
   if ( this->getProcessType() == MASS_TRANSPORT )
   {
       for (int i = 0; i < pcs_number_of_primary_nvals; i++)
       {
          int nidx = GetNodeValueIndex(pcs_primary_function_name[i]);

          // PrimaryVariable pv_i (convertPrimaryVariable(pcs_primary_function_name[i]));
          for (size_t j = 0; j < ic_vector.size(); j++)
          {
             m_ic = ic_vector[j];
             m_ic->m_msh = m_msh;                     //OK/MX

             if (m_ic->getProcess() == this )
             {
                m_ic->Set(nidx);
                m_ic->Set(nidx + 1);
             }
          }
       }
   }
   else // otherwise PrimaryVariable check is still performed.
   {
          for (int i = 0; i < pcs_number_of_primary_nvals; i++)
          {
              int nidx = GetNodeValueIndex(pcs_primary_function_name[i]);
              PrimaryVariable pv_i (convertPrimaryVariable(pcs_primary_function_name[i]));
              for (size_t j = 0; j < ic_vector.size(); j++)
              {
                  m_ic = ic_vector[j];
                  m_ic->m_msh = m_msh;                     //OK/MX

                  if (m_ic->getProcessType() != this->getProcessType())
                      continue;

                  m_ic->setProcess(this);
                  if (m_ic->getProcessPrimaryVariable() == pv_i)
                  {
                      m_ic->Set(nidx);
                      m_ic->Set(nidx + 1);
                  } // end of if
              } // end of for j
          } // end of for i

   } // end of if-else

}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2005 OK Implementation
last modified:
**************************************************************************/
#ifndef NEW_EQS                                   //WW. 07.11.2008
void CRFProcess::SetNODValues()
{
   for(long i=0;i<(long)m_msh->nod_vector.size();i++)
   {
      //    SetNODValue(i,GetNODValueIndex(_pcs_type_name),eqs->x[i]);
      //    SetNODValue(i,GetNODValueIndex(_pcs_type_name)+1,eqs->x[i]);
                                                  //WW
      SetNodeValue(m_msh->Eqs2Global_NodeIndex[i],0,eqs->x[i]);
                                                  //WW
      SetNodeValue(m_msh->Eqs2Global_NodeIndex[i],1,eqs->x[i]);
   }
}
#endif

/**************************************************************************
FEMLib-Method:
Task: Ermittelt den Fehler bei Iterationen
 new_iteration     : Vektor des neuen Iterationsschritts
 old_iteration_ndx : Knotenindex fuer Werte des alten Iterationsschritts
 reference_ndx     : Knotenindex fuer Werte des alten Zeitschritts (als Referenz)
 method            : Methode der Fehlerermittlung
Programing:
01/2005 OK NUM implementation
05/2005 OK MSH
08/2005 WW Re-implememtation based on NUMCalcIterationError
01/2007 WW For DOF>1
11/2007 WW Changes for the new classes of sparse and linear solver
last modification:
**************************************************************************/
double CRFProcess::CalcIterationNODError(int method)
{

   static long i, k, g_nnodes;
   static double error, change, max_c, min_c;
   //-----------------------------------------------------WW
   int nidx1;
   int ii;
   g_nnodes = m_msh->GetNodesNumber(false);       //WW
   //-----------------------------------------------------WW
   double *eqs_x = NULL;                          // 11.2007. WW
#ifdef NEW_EQS
   eqs_x = eqs_new->x;
#else
   eqs_x = eqs->x;
#endif
   error = 0.;
   change = 0.;

   max_c = 0.;
   min_c = 1.e99;

   switch (method)
   {
      default:
      case 0:
         return 0.;
         // Maximum error
      case 1:
                                                  //WW
         for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
         {
                                                  //new time
            nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
            for (i = 0l; i < g_nnodes; i++)       //WW
            {
               k = m_msh->Eqs2Global_NodeIndex[i];
               error = max(error, fabs(GetNodeValue(k, nidx1) - eqs_x[i+ii*g_nnodes]));
            }
         }
         return error;
      case 2:
                                                  //WW
         for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
         {
                                                  //new time
            nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
            for (i = 0l; i < g_nnodes; i++)       //WW
            {
               k = m_msh->Eqs2Global_NodeIndex[i];
               error = max(error, fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k, nidx1))
                  /(fabs(eqs_x[i+ii*g_nnodes]) + fabs(GetNodeValue(k, nidx1))+ MKleinsteZahl) );
            }
         }
         return error;
      case 3:
                                                  //WW
         for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
         {
                                                  //new time
            nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
            for (i = 0l; i < g_nnodes; i++)       //WW
            {
               k = m_msh->Eqs2Global_NodeIndex[i];
               error = max(error, fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k, nidx1)));
               max_c = max(max(max_c, fabs(fabs(eqs_x[i+ii*g_nnodes]))),fabs(GetNodeValue(k, nidx1)));
            }
         }
         return error / (max_c + MKleinsteZahl);
      case 4:
                                                  //WW
         for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
         {
                                                  //new time
            nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
            for (i = 0l; i < g_nnodes; i++)       //WW
            {
               k = m_msh->Eqs2Global_NodeIndex[i];
               error = max(error, fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k, nidx1)));
               min_c = min(min_c, fabs(eqs_x[i+ii*g_nnodes]));
               max_c = max(max_c, fabs(eqs_x[i+ii*g_nnodes]));
            }
         }
         return error / (max_c - min_c + MKleinsteZahl);
      case 5:
                                                  //WW
         for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
         {
                                                  //new time
            nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
            for (i = 0l; i < g_nnodes; i++)       //WW
            {
               k = m_msh->Eqs2Global_NodeIndex[i];
               error = max(error, fabs(eqs_x[i+ii*g_nnodes] -  GetNodeValue(k, nidx1))
                  / (fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k, nidx1-1)) + MKleinsteZahl));
            }
         }
         return error;
      case 6:
                                                  //WW
         for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
         {
                                                  //new time
            nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
            for (i = 0l; i < g_nnodes; i++)       //WW
            {
               k = m_msh->Eqs2Global_NodeIndex[i];
               error = max(error, fabs(eqs_x[i+ii*g_nnodes] -  GetNodeValue(k, nidx1)));
               change = max(change, fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k, nidx1-1)));
            }
         }
         return error / (change + MKleinsteZahl);
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 OK Implementation
04/2006 YD Add contiuum coupling OK???Why ere
04/2007 WW Remove the spurious stuff
08/2008 WW Time step size control (First)
12/2008 WW Time step size control (Update)
**************************************************************************/
double CRFProcess::ExecuteNonLinear()
{
   double nonlinear_iteration_error=0.0;
#ifdef NEW_EQS
   configured_in_nonlinearloop = true;
   // Also allocate temporary memory for linear solver. WW
   //
#if defined(USE_MPI)
   CPARDomain *dom = dom_vector[myrank];
   dom->eqs->SetDOF(pcs_number_of_primary_nvals); //_new 02/2010 WW
   dom->ConfigEQS(m_num, pcs_number_of_primary_nvals
      *m_msh->GetNodesNumber(false));
#else
   eqs_new->SetDOF(pcs_number_of_primary_nvals);  //_new 02/2010. WW
   eqs_new->ConfigNumerics(m_num);
#endif
   //
#endif
   //..................................................................
   // PI time step size control. 29.08.2008. WW
   if(Tim->GetTimeStepCrtlType()>0)
      CopyU_n(aproblem->GetBufferArray());
   if (hasAnyProcessDeactivatedSubdomains)
      this->CheckMarkedElement();                 //NW
   for(iter_nlin=0;iter_nlin<pcs_nonlinear_iterations;iter_nlin++)
   {
      cout << "    PCS non-linear iteration: " << iter_nlin << "/"
         << pcs_nonlinear_iterations << endl;
      nonlinear_iteration_error = Execute();
      if(mobile_nodes_flag ==1)
         PCSMoveNOD();
      if(!Tim)                                    //OK
         continue;
      if(nonlinear_iteration_error < pcs_nonlinear_iteration_tolerance)
      {
         Tim->repeat = false;                     //OK/YD
                                                  //OK/YD
         Tim->nonlinear_iteration_error = nonlinear_iteration_error;
         break;
      }
      else
         Tim->repeat = true;                      //OK/YD
   }
   // PI time step size control. 27.08.2008. WW
   if(Tim->GetTimeStepCrtlType()>0)
      PI_TimeStepSize(aproblem->GetBufferArray());
   // 8 Calculate secondary variables
   if(accepted)                                   // 27.08.2008. WW
      CalcSecondaryVariables();                   // Moved here from Execute() WW
   // Release temporary memory of linear solver. WW
#ifdef NEW_EQS                                 //WW
#if defined(USE_MPI)
   dom->eqs->Clean();
#else
   eqs_new->Clean();                              // Release buffer momery WW
#endif
   configured_in_nonlinearloop = false;
#endif
   return nonlinear_iteration_error;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
02/2007 WW Implementation
**************************************************************************/
void CRFProcess::Extropolation_GaussValue()
{
   int k, NS;
   long i = 0;
   int idx[3];
   // const long LowOrderNodes= m_msh->GetNodesNumber(false);
   Mesh_Group::CElem* elem = NULL;

   //
   NS =m_msh->GetCoordinateFlag()/10;
   idx[0] = GetNodeValueIndex("VELOCITY_X1");
   idx[1] = GetNodeValueIndex("VELOCITY_Y1");
   idx[2] = GetNodeValueIndex("VELOCITY_Z1");
   for (i = 0; i <(long)m_msh->GetNodesNumber(false); i++)
   {
      for(k=0; k<NS; k++)
         SetNodeValue(i, idx[k], 0.0);
   }
   if(type==1212 || type==1313)                   // Multi-phase flow
   {
      idx[0] = GetNodeValueIndex("VELOCITY_X2");
      idx[1] = GetNodeValueIndex("VELOCITY_Y2");
      idx[2] = GetNodeValueIndex("VELOCITY_Z2");
      for (i = 0; i <(long)m_msh->GetNodesNumber(false); i++)
      {
         for(k=0; k<NS; k++)
            SetNodeValue(i, idx[k], 0.0);
      }
   }
   //
   for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
   {
      elem = m_msh->ele_vector[i];
      if (elem->GetMark())                        // Marked for use
      {
         fem->ConfigElement(elem);
         for(k=0; k<NS; k++)
            fem->ExtropolateGauss(this, k);
      }
   }
}


/**************************************************************************
FEMLib-Method:
Task:    Calculate the material values at Gauss points and extropolate them
         to element node
Programing:
04/2007 WW Implementation
**************************************************************************/
void CRFProcess::Extropolation_MatValue()
{
   //	if (_pcs_type_name.find("FLOW") == string::npos)
   if (! isFlowProcess(this->getProcessType()))
      return;
   if (additioanl2ndvar_print < 0)
      return;

   //
   int NS = m_msh->GetCoordinateFlag() / 10;
   //
   if ((additioanl2ndvar_print > 0) && (additioanl2ndvar_print < 3))
   {
      int idx[3];
      idx[0] = GetNodeValueIndex("PERMEABILITY_X1");
      idx[1] = GetNodeValueIndex("PERMEABILITY_Y1");
      if (NS > 2)
         idx[2] = GetNodeValueIndex("PERMEABILITY_Z1");
      for (long i = 0; i < m_msh->GetNodesNumber(false); i++)
      {
         for (int k = 0; k < NS; k++)
            SetNodeValue(i, idx[k], 0.0);
      }
   }
   if (additioanl2ndvar_print > 1)
   {
      int idxp = GetNodeValueIndex("POROSITY");
      for (long i = 0; i < m_msh->GetNodesNumber(false); i++)
         SetNodeValue(i, idxp, 0.0);
   }
   //
   continuum = 0;
   if (continuum_vector.size() == 2)
      continuum = 1;

   Mesh_Group::CElem* elem = NULL;
   for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
   {
      elem = m_msh->ele_vector[i];
      if (elem->GetMark())                        // Marked for use
      {
         fem->ConfigElement(elem);
         fem->CalcNodeMatParatemer();
      }
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
04/2005 OK Implementation
**************************************************************************/
void PCSDelete(const std::string &m_pcs_type_name)
{
   ProcessType pcs_type (convertProcessType (m_pcs_type_name));
   CRFProcess* m_pcs = NULL;
   for (size_t i = 0; i < pcs_vector.size(); i++)
   {
      m_pcs = pcs_vector[i];
      //		if (m_pcs->_pcs_type_name.compare(m_pcs_type_name) == 0) { TF
      if (m_pcs->getProcessType () == pcs_type)
      {
         delete m_pcs;
         pcs_vector.erase(pcs_vector.begin() + i);
      }
   }
}


/**************************************************************************
   GeoSys - Function: Reallocation

   Aufgabe:
        Reallocte memory by new operator
   09/2005   WW    Erste Version

**************************************************************************/
template <class T> T *resize(T *array, size_t old_size, size_t new_size)
{
   T *temp = new T[new_size];
   for(size_t i=0; i<old_size; i++) temp[i] = array[i];
   for(size_t i=old_size; i<new_size; i++) temp[i] = 0;
   delete [] array;
   array = temp;
   return temp;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
11/2005 MB Implementation
**************************************************************************/
#ifndef NEW_EQS                                   //WW. 07.11.2008
void CRFProcess::CalcFluxesForCoupling(void)
{
   int i, j;
   double flux;
   long n_index;
   long NodeIndex_GW;
   long NodeIndex_OLF;
   int no_nodes = m_msh->getNumberOfMeshLayers() + 1;
   long no_richards_problems = (long) (m_msh->ele_vector.size()
      / m_msh->getNumberOfMeshLayers());
   long IndexBottomNode;
   long IndexTopNode;
   int NoOfGWNodes = 0;
   double AverageZ_GW = 0.0;
   double AverageZ_OLF = 0.0;
   double AverageH_GW = 0.0;
   double AverageH_OLF = 0.0;
   double dh;
   int idxFLUX;
   int idxHead_GW;
   int idxHead_OLF;
   Mesh_Group::CElem* m_ele_GW = NULL;
   Mesh_Group::CElem* m_ele_OLF = NULL;

   // Get processes
   CRFProcess* m_pcs_GW(PCSGet(GROUNDWATER_FLOW));
   if (!m_pcs_GW)                                 //OK
   {
      cout << "Fatal error: no GROUNDWATER_FLOW process" << endl;
      return;
   }
   CRFProcess* m_pcs_OLF(PCSGet(OVERLAND_FLOW));
   if (!m_pcs_OLF)                                //OK
   {
      cout << "Fatal error: no OVERLAND_FLOW process" << endl;
      return;
   }

   // Get meshes
   CFEMesh* m_msh_GW = m_pcs_GW->m_msh;
   CFEMesh* m_msh_OLF = m_pcs_OLF->m_msh;

   //Get indeces
   idxHead_GW = m_pcs_GW->GetNodeValueIndex("HEAD") + 1;
   idxHead_OLF = m_pcs_OLF->GetNodeValueIndex("HEAD") + 1;
   idxFLUX = GetNodeValueIndex("FLUX") + 1;

   for (i = 0; i < no_richards_problems; i++)
   {
      IndexBottomNode = ((i + 1) * no_nodes) - 1;

      // ToDo safe somewhere else so that this has to be done only once
      //-----------------------------------------------------------------
      // Get Nearest GW and OLF Element
      GEOLIB::Point pnt;
      pnt[0] = m_msh->nod_vector[IndexBottomNode]->X();
      pnt[1] = m_msh->nod_vector[IndexBottomNode]->Y();
      pnt[2] = m_msh->nod_vector[IndexBottomNode]->Z();

      long EleNumber = m_msh_GW->GetNearestELEOnPNT(&pnt);

      //GW and OLF use the same Numbering !!!
      m_ele_GW = m_msh_GW->ele_vector[EleNumber];
      m_ele_OLF = m_msh_OLF->ele_vector[EleNumber];

      //-----------------------------------------------------------------
      // Get Average values for element //ToDo encapsulate //WW: CElement::elemnt_averag??e
      NoOfGWNodes = m_ele_OLF->GetNodesNumber(m_msh_GW->getOrder());
      for (j = 0; j < NoOfGWNodes; j++)
      {
         NodeIndex_GW = m_ele_GW->GetNodeIndex(j);
         NodeIndex_OLF = m_ele_OLF->GetNodeIndex(j);

         AverageZ_GW += m_pcs_GW->GetNodeValue(NodeIndex_GW, idxHead_GW);
         AverageZ_OLF += m_msh_OLF->nod_vector[NodeIndex_OLF]->Z();
      }
      AverageZ_GW = AverageZ_GW / NoOfGWNodes;
      AverageZ_OLF = AverageZ_OLF / NoOfGWNodes;

      //-----------------------------------------------------------------
      // UsatZone exists -> Flux from this
      if (AverageZ_GW < AverageZ_OLF)
      {
         n_index = m_msh->Eqs2Global_NodeIndex[IndexBottomNode];
         if (m_msh->nod_vector[IndexBottomNode]->GetMark())
         {
            flux = eqs->b[IndexBottomNode];
            //FLUXES IN NEW VERSION WITH VELOCITIES !!!!!
            //WAIT FOR SEBASTIANS MASS TRANSPORT IN USAT ZONE !!!!!
            //TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
            flux = 0.00001;
            //flux = 1;
            SetNodeValue(n_index, idxFLUX, flux);
         }
      }

      //-----------------------------------------------------------------
      // No UsatZone -> Calculate Flux from leakage terms
      if (AverageZ_GW >= AverageZ_OLF)
      {
         //SetRichardsNodesToFullySaturated??
         IndexTopNode = i * no_nodes;             // Top Node of Richard Column -> Flux for OLF
                                                  // Bottom Node of Richards Column -> Flux for GW
         IndexBottomNode = ((i + 1) * no_nodes) - 1;
         //-----------------------------------------------------------------
         // Get Average values for element //ToDo encapsulate
         for (j = 0; j < NoOfGWNodes; j++)
         {
            NodeIndex_GW = m_ele_GW->GetNodeIndex(j);
            NodeIndex_OLF = m_ele_OLF->GetNodeIndex(j);
            AverageH_GW += m_pcs_GW->GetNodeValue(NodeIndex_GW, idxHead_GW);
            AverageH_OLF += m_pcs_OLF->GetNodeValue(NodeIndex_OLF,
               idxHead_OLF);
         }
         AverageH_GW = AverageH_GW / NoOfGWNodes;
         AverageH_OLF = AverageH_OLF / NoOfGWNodes;
         //Calculate the vertical leakage
         dh = AverageH_GW - AverageH_OLF;
         // get kf fully saturated of uppermost element ?
         // or user defined value: entry resistance / leakage factor ?
         //flux = dh * 0.001;
         flux = dh * 1.;

         //1. Add reacharge value to GW flow -> Add to flux off IndexBottomNode
         //Achtung nur zum Testen Source f\A8\B9r GW flow durchgehend !!!!!!
         //SetNodeValue(IndexBottomNode, idxFLUX, -flux);  //H_OLF  > H_GW -> + flux_GW
         SetNodeValue(IndexBottomNode, idxFLUX, 0.00001);

         //2. Add reacharge value to OLF -> Add to flux off IndexTopNode
                                                  //H_OLF  > H_GW -> - flux_OLF
         SetNodeValue(IndexTopNode, idxFLUX, flux);
         //3. Set flag to set reacharge to Usat to zero ???

      }
   }
}


/**************************************************************************
FEMLib-Method:
Task: Ermittelt den Fehler bei Kopplungs Iterationen
Programing:
11/2005 MB implementation
**************************************************************************/
double CRFProcess::CalcCouplingNODError()
{
   static long i, k;
   static double error, change, max_c, min_c;
   int ndx0 = GetNodeValueIndex(m_num->cpl_variable);
   int ndx1 = GetNodeValueIndex(m_num->cpl_variable)+1;

   error = 0.;
   change = 0.;

   max_c = 0.;
   min_c = 1.e99;

   for (i = 0l; i < eqs->dim; i++)
   {
      k = m_msh->Eqs2Global_NodeIndex[i];
      if(m_msh->nod_vector[k]->GetMark())
         error = max(error, fabs(GetNodeValue(k, ndx1) - GetNodeValue(k, ndx0)));
   }
   cout <<  "  Coupling error: " << error << endl;

   return error;

}
#endif
/**************************************************************************
FEMLib-Method:
Task:
Programing:
11/2005 MB implementation
**************************************************************************/
void CRFProcess::CopyCouplingNODValues()
{
   //Carefull if cpl_variable = primary variable -> need extra column in NodeValueTable !
   int nidx0 = GetNodeValueIndex(m_num->cpl_variable);
   for (long l = 0; l < m_msh->GetNodesNumber(false); l++)
   {
      SetNodeValue(l, nidx0, GetNodeValue(l, nidx0+1));
   }
   //	if (_pcs_type_name.find("RICHARDS") != string::npos) { //WW
   if (this->getProcessType() == RICHARDS_FLOW)   //WW
   {
      nidx0 = GetNodeValueIndex("SATURATION1");
      for (long l = 0; l < m_msh->GetNodesNumber(false); l++)
         SetNodeValue(l, nidx0, GetNodeValue(l, nidx0 + 1));
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
11/2005 MB implementation
02/2006 WW Modified for the cases of high order element and saturation
08/2008 WW Make it twofold copy: forward and backward
**************************************************************************/
void CRFProcess::CopyTimestepNODValues(bool forward)
{
   bool Quadr = false;                            //WW
   if (type == 4 || type == 41)
      Quadr = true;

   for (int j = 0; j < pcs_number_of_primary_nvals; j++)
   {
      int nidx0 = GetNodeValueIndex(pcs_primary_function_name[j]);
      int nidx1 = nidx0 + 1;
      if (!forward)                               //08.2008. WW
      {
         nidx0++;
         nidx1--;
      }
      for (long l = 0; l < m_msh->GetNodesNumber(Quadr); l++)
         SetNodeValue(l, nidx0, GetNodeValue(l, nidx1));
      //WW
      //		if (_pcs_type_name.find("RICHARDS") != string::npos || type == 1212) { //Multiphase. WW
                                                  //Multiphase. WW
      if (this->getProcessType() == RICHARDS_FLOW || type == 1212)
      {
         if (j == 1 && type == 1212)              // Multiphase. WW
            continue;
         if (j == 0)
            nidx0 = GetNodeValueIndex("SATURATION1");
         else if (j == 1)
            nidx0 = GetNodeValueIndex("SATURATION2");
         nidx1 = nidx0 + 1;
         if (!forward)                            //27.08.2008. WW
         {
            nidx0++;
            nidx1--;
         }
         //
         for (long l = 0; l < m_msh->GetNodesNumber(false); l++)
            SetNodeValue(l, nidx0, GetNodeValue(l, nidx1));
      }
   }
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
03/2006 SB Implementation
last modified:
**************************************************************************/
int PCSGetPCSIndex(const string &pcs_type_name, const string &comp_name)
{
   ProcessType pcs_type (convertProcessType (pcs_type_name));

   CRFProcess *m_pcs = NULL;
   int i, pcs_no;
   int no_processes = (int) pcs_vector.size();
   string testname;
   pcs_no = -1;
   for (i = 0; i < no_processes; i++)
   {
      m_pcs = pcs_vector[i];
      //		if (m_pcs->pcs_type_name.compare(pcs_type_name) == 0) {
      if (m_pcs->getProcessType () == pcs_type)
      {
         testname = m_pcs->pcs_primary_function_name[0];
         if (testname.compare(comp_name) == 0)
         {
            //        cout << " Found in PCSGetbyTypeandCompName for PCSType/Compname " << pcs_type_name << ", " << comp_name;
            //        cout << " Process number " << m_pcs->pcs_number << ", compnumber " << m_pcs->pcs_component_number << endl;
            pcs_no = i;
            return pcs_no;
         }
      }
   }
   return pcs_no;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
12/2005 SB Implementation
last modified:
10/2010 TF restructured function
**************************************************************************/
CRFProcess* PCSGet(const std::string &pcs_type_name, const std::string &comp_name)
{
   ProcessType pcs_type (convertProcessType (pcs_type_name));
   size_t no_processes (pcs_vector.size());
   for (size_t i = 0; i < no_processes; i++)
   {
      //		if (pcs_vector[i]->pcs_type_name.compare(pcs_type_name) == 0) { // TF
      if (pcs_vector[i]->getProcessType () == pcs_type)
      {
         if (comp_name.compare (pcs_vector[i]->pcs_primary_function_name[0]) == 0)
         {
            return pcs_vector[i];
         }
      }
   }
   return NULL;
}


CRFProcess* PCSGet(ProcessType pcs_type, const std::string &comp_name)
{
   size_t no_processes (pcs_vector.size());
   for (size_t i = 0; i < no_processes; i++)
   {
      if (pcs_vector[i]->getProcessType () == pcs_type)
      {
         if (comp_name.compare (pcs_vector[i]->pcs_primary_function_name[0]) == 0)
         {
            return pcs_vector[i];
         }
      }
   }
   return NULL;
}


/**************************************************************************
PCSLib-Method:
12/2005 OK Implementation
**************************************************************************/
CRFProcess* PCSGet(const string &var_name,bool bdummy)
{
   string pcs_var_name;
   CRFProcess *m_pcs = NULL;
   bdummy = bdummy;                               //WW
   for(size_t i=0;i<pcs_vector.size();i++)
   {
      m_pcs = pcs_vector[i];
      for(size_t j=0;j<m_pcs->GetPrimaryVNumber();j++)
      {
         pcs_var_name = m_pcs->pcs_primary_function_name[j];
         if(pcs_var_name.compare(var_name)==0)
         {
            return m_pcs;
         }
      }
      for(size_t j=0;j<m_pcs->GetSecondaryVNumber();j++)
      {
         pcs_var_name = m_pcs->pcs_secondary_function_name[j];
         if(pcs_var_name.compare(var_name)==0)
         {
            return m_pcs;
         }
      }
   }
   return NULL;
}


/**************************************************************************
PCSLib-Method:
05/2006 CMCD Implementation
**************************************************************************/
CRFProcess* PCSGetFluxProcess()
{
   CRFProcess *m_pcs = NULL;
   bool found = false;
   const size_t no_processes (pcs_vector.size());

   for (size_t i = 0; i < no_processes; i++)
   {
      //		if (pcs_vector[i]->_pcs_type_name == "LIQUID_FLOW") { // TF
      if (pcs_vector[i]->getProcessType () == LIQUID_FLOW)
      {
         m_pcs = pcs_vector[i];
         found = true;
      }
      //		if (pcs_vector[i]->_pcs_type_name == "GROUNDWATER_FLOW") {
      if (pcs_vector[i]->getProcessType() == GROUNDWATER_FLOW)
      {
         m_pcs = pcs_vector[i];
         found = true;
      }
      if (found)
         return m_pcs;
   }
   return NULL;
}


/*************************************************************************
GeoSys-FEM Function:
Task:
For Richards, dual-Richards and multi-phase flow.
Programming:
02/2007 WW Implementation
**************************************************************************/
void CRFProcess::CalcSecondaryVariablesUnsaturatedFlow(bool initial)
{
   int ii;
   long i;
   double p_cap;
   int idxp,idxcp,idxS;
   //----------------------------------------------
   vector<string> secondSNames;
   vector<string> secondCPNames;
   secondSNames.push_back("SATURATION1");
   if(type==1212)                                 // Multiphase flow
      secondCPNames.push_back("PRESSURE1");
   else
      secondCPNames.push_back("PRESSURE_CAP1");
   if(continuum_vector.size()==2)
   {
      secondSNames.push_back("SATURATION2");
      secondCPNames.push_back("PRESSURE_CAP2");
   }

   //
   CElem* elem =NULL;
   CFiniteElementStd* fem = GetAssembler();
   //----------------------------------------------------------------------
   for(ii=0;ii<(int)secondSNames.size();ii++)
   {
      idxS = GetNodeValueIndex(secondSNames[ii].c_str())+1;
      if(type==1212)                              // Multiphase flow
      {
         idxcp = GetNodeValueIndex(pcs_primary_function_name[ii])+1;
         idxp = GetNodeValueIndex("PRESSURE_W");
      }
      else
      {
         idxp = GetNodeValueIndex(pcs_primary_function_name[ii])+1;
         idxcp = GetNodeValueIndex(secondCPNames[ii].c_str());
      }
      //----------------------------------------------------------------------
      // Capillary pressure
      for(i=0;i<(long)m_msh->GetNodesNumber(false);i++)
      {
         // Copy the new saturation to the old level
         //
         if(type==1212)                           // Multiphase flow
         {
            // p_w = p_g-p_c
            p_cap = GetNodeValue(i,idxcp+2)-GetNodeValue(i,idxcp);
            SetNodeValue(i,idxp,p_cap);
         }
         else
         {
            p_cap = -GetNodeValue(i,idxp);
            SetNodeValue(i,idxcp,p_cap);
         }
         SetNodeValue(i,idxS, 0.0);
      }
   }
   // Cal S
   for(ii=0;ii<(int)secondSNames.size();ii++)
   {
      continuum = ii;
      for(i=0;i<(long)m_msh->ele_vector.size();i++)
      {
         elem = m_msh->ele_vector[i];
         if (elem->GetMark())                     // Marked for use
         {
            elem->SetOrder(false);
            fem->ConfigElement(elem,false);
            fem->CalcSatution();
         }
      }
   }
   //
   if(!initial)
      return;
   //----------
   for(ii=0;ii<(int)secondSNames.size();ii++)
   {
      idxS = GetNodeValueIndex(secondSNames[ii].c_str())+1;
      for(i=0;i<(long)m_msh->GetNodesNumber(false);i++)
         SetNodeValue(i,idxS-1,  GetNodeValue(i, idxS));
   }
}


/*************************************************************************
GeoSys-FEM Function:
Task: Updating secondary variables for Multiphase flow in PS_GLOBAL

Programming:
03/2009 PCH Implementation
**************************************************************************/
void CRFProcess::CalcSecondaryVariablesPSGLOBAL()
{
   long i;
   int ndx_pressure1, ndx_p_cap, ndx_pressure2, ndx_s_wetting, ndx_s_nonwetting;

   // The primary variables
   ndx_pressure1 = GetNodeValueIndex("PRESSURE1");
   ndx_s_nonwetting = GetNodeValueIndex("SATURATION2");

   // The secondary variables
   ndx_pressure2 = GetNodeValueIndex("PRESSURE2");
   ndx_p_cap = GetNodeValueIndex("PRESSURE_CAP");
   ndx_s_wetting = GetNodeValueIndex("SATURATION1");

   double pressure1, pressure2, p_cap, s_wetting;
   for(i=0;i<(long)m_msh->GetNodesNumber(false);i++)
   {
      pressure1 = GetNodeValue(i,ndx_pressure1+1);// New
      pressure2 = GetNodeValue(i,ndx_pressure1);  // Old

      // Let's get capillary pressure before updating pressure2
      // by accessing the primary variable of the saturation equation
      // not the secondary variable of it.
      int ndx_sat2 = GetNodeValueIndex("SATURATION2");
      double sat2 = GetNodeValue(i,ndx_sat2+1);
      // Due to the iterative solution scheme in solving Snw with no
      // explicit boundary condition for non-zero flux condition,
      // Snw may become negative particularly the density difference
      // between two fluids is big. To prevent negative Snw, the
      // saturation restriction added.
      CMediumProperties* mmp = NULL;
      if(mmp_vector.size() > 1)
      {
         double sum = 0.0;
         CNode* thisNode = m_msh->nod_vector[i];
         int NumOfNeighborElements = (int)thisNode->connected_elements.size();
         // Harmonic mean
         for(int i=0; i< NumOfNeighborElements; ++i)
         {
            // Mount neighboring elemenets and get the corresponding material group one by one.
            int eleIdx = thisNode->connected_elements[i];
            CElem* thisEle = m_msh->ele_vector[eleIdx];
            int matgrp = thisEle->GetPatchIndex();
            mmp = mmp_vector[matgrp];
            mmp->mode = 2;
            sum += 1.0/sat2;
            //			sum += 1.0/MRange(mmp->saturation_res[1],sat2,1.0-mmp->saturation_res[0]);
         }
         sat2 = (double)NumOfNeighborElements/sum;
      }
      else
      {
         mmp = mmp_vector[0];
         //		sat2 = MRange(mmp->saturation_res[1],sat2,1.0-mmp->saturation_res[0]);
      }
      s_wetting = 1.0 - sat2;
      // Assigning the secondary variable, Sw
      SetNodeValue(i,ndx_s_wetting,s_wetting);
      // Assigning the primary variable Snw here one more time
      // to completely bound the range of saturation
      //	SetNodeValue(i,ndx_s_nonwetting,sat2);
      //	SetNodeValue(i,ndx_s_nonwetting+1,sat2);

      // Assigning the secondary variable, Pc
      if(mmp_vector.size() > 1)
         p_cap = GetCapillaryPressureOnNodeByNeighobringElementPatches(i, 2, 1.0-sat2);
      else
         p_cap = mmp->CapillaryPressureFunction(i,NULL,1.0,0,1.0-sat2);

      SetNodeValue(i,ndx_p_cap,p_cap);

      pressure2 = pressure1 + p_cap;
      // Assigning the secondary variables, Pnw
      SetNodeValue(i,ndx_pressure2,pressure2);
   }
}


/**************************************************************************
FEMLib-Method:
Task: Calculate saturation on node by averaging the patches of the
         neighboring elements in three means

  0: Arithmatic mean
  1: Geomtric mean
  2: Harmonic mean

Programing:
03/2009 PCH Implementation
last modification:
*************************************************************************/
double CRFProcess::GetCapillaryPressureOnNodeByNeighobringElementPatches(int nodeIdx, int meanOption, double Sw)
{
   double p_cap = 0.0, sum = 0.0;

   CNode* thisNode = m_msh->nod_vector[nodeIdx];
   int NumOfNeighborElements = (int)thisNode->connected_elements.size();

   switch (meanOption)
   {
      case 0:
         break;
      case 1:
         break;
      case 2:                                     // Harmonic mean
         for(int i=0; i< NumOfNeighborElements; ++i)
         {
            // Mount neighboring elemenets and get the corresponding material group one by one.
            int eleIdx = thisNode->connected_elements[i];
            CElem* thisEle = m_msh->ele_vector[eleIdx];
            int matgrp = thisEle->GetPatchIndex();
            CMediumProperties* mmp = mmp_vector[matgrp];
            sum += 1.0/mmp->CapillaryPressureFunction(i,NULL,1.0,0,Sw);
         }
         p_cap = (double)NumOfNeighborElements/sum;
         break;

      default:
         cout<<"Please define the option for various means!"<<endl;
         cout<<"The code stops at GetCapillaryPressureOnNodeByNeighobringElementPatches function!"<<endl;
         abort();
         break;
   }

   return p_cap;
}


/*************************************************************************
GeoSys-FEM Function:
Task: Calculates saturation for richards flow,
      alternative to CRFProcess::CalcSecondaryVariablesRichards
      uses pressure of only one node
     evoked by switch case in pcs-file: SATURATION_SWITCH = true
Programming:
06/2007 JOD Implementation
**************************************************************************/
void CRFProcess::CalcSaturationRichards(int timelevel, bool update)
{
   double p_cap, saturation, volume_sum;
   int idxp, idxcp, idxS, idx_tS = -1;
   const size_t number_continuum (continuum_vector.size());
   size_t i_s, i_e;

   CMediumProperties* m_mmp = NULL;
   CElem* elem = NULL;
   //WW  CFiniteElementStd* fem = GetAssembler();

   if (continuum_ic)                              //Create IC: for both continua
   {
      i_s = 0;
      i_e = number_continuum;
   }
   else
   {
      i_s = continuum;
      i_e = continuum + 1;
   }

   for (size_t i_pv = i_s; i_pv < i_e; i_pv++)
   {
      idxp = GetNodeValueIndex(pcs_primary_function_name[i_pv]) + timelevel;
      idxS = GetNodeValueIndex(pcs_secondary_function_name[i_pv
         * number_continuum]) + timelevel;
      idxcp = GetNodeValueIndex(pcs_secondary_function_name[i_pv
         * number_continuum + number_continuum * 2]) + timelevel;
      if (continuum_vector.size() > 1)
         idx_tS = GetNodeValueIndex("TOTAL_SATURATION") + timelevel;

      for (long i = 0; i < (long) m_msh->GetNodesNumber(false); i++)
      {
         // Capillary pressure
         p_cap = -GetNodeValue(i, idxp);
         if (timelevel == 1 && update)
            SetNodeValue(i, idxcp - 1, GetNodeValue(i, idxcp));
         SetNodeValue(i, idxcp, p_cap);
         if (timelevel == 1 && update)
            SetNodeValue(i, idxS - 1, GetNodeValue(i, idxS));

         // Liquid saturation
         if (continuum_vector.size() > 1)
            SetNodeValue(i, idx_tS, 0.0);
         //

         saturation = 0., volume_sum = 0.;
         size_t elemsCnode = m_msh->nod_vector[i]->connected_elements.size();

         for (size_t j = 0; j < elemsCnode; j++)
         {
            elem = m_msh->ele_vector[m_msh->nod_vector[i]->connected_elements[j]];
            m_mmp = mmp_vector[elem->GetPatchIndex()];
            volume_sum += elem->volume;
            saturation += m_mmp->SaturationCapillaryPressureFunction(
               p_cap,
               (int) mfp_vector.size() - 1) * elem->volume;
         }
         saturation /= volume_sum;
         SetNodeValue(i, idxS, saturation);

      }

   }

   if (continuum > 0)
   {
      for (long i = 0; i < (long) m_msh->GetNodesNumber(false); i++)
      {
         double total_S = 0;
         for (size_t j = 0; j < continuum_vector.size(); j++)
         {
            idxS = GetNodeValueIndex(pcs_secondary_function_name[j
               * number_continuum]) + timelevel;
            total_S += GetNodeValue(i, idxS) * continuum_vector[j];
         }
         SetNodeValue(i, idx_tS, total_S);
      }
   }
}


/**************************************************************************
   GeoSys - Function: Get mean element value for element index from secondary node values
                  of process pcs_name and for variable var_name; old and new timelevel
    01/2006   SB    Implementation
    02/2008   CB    generalization
**************************************************************************/
double PCSGetEleMeanNodeSecondary_2(long index, int pcsT, const string &var_name, int timelevel)
{

   double val = 1.0;                              // As this returns saturation, default is fully saturated = 1.0;
   int idx=0, j;                                  //OK411
   long enode;
   CRFProcess *m_pcs = NULL;
   CRFProcess *cplpcs = NULL;
   CElem* elem =NULL;

   // Get index of secondary node value
   switch(pcsT)
   {
      case 0:                                     //Liquid_Flow
         break;
      case 1:                                     //Groundwater Flow
         break;
      case 66:                                    //Overland Flow
         break;
      case 5:                                     //Air Flow
         break;
      case 11:                                    //Componental Flow
         break;
      case 1212:                                  //Multiphase Flow
         break;
      case 12:                                    // Two_phase_Flow
         m_pcs = PCSGet(TWO_PHASE_FLOW);
         if(m_pcs)
         {
            if(m_pcs->pcs_type_number==0)
               cplpcs = pcs_vector[m_pcs->pcs_number+1];
            else if(m_pcs->pcs_type_number==1)
               cplpcs = pcs_vector[m_pcs->pcs_number-1];
            idx = cplpcs->GetNodeValueIndex(var_name)+timelevel;
         }
         break;
      case 22:                                    // Richards flow
         m_pcs = PCSGet(RICHARDS_FLOW);
         if(m_pcs)
         {
            idx = m_pcs->GetNodeValueIndex(var_name)+timelevel;
            cplpcs = m_pcs;
         }
         break;
      default:
         break;
   }

   if(m_pcs)
   {
      // Get element with index index
      elem = m_pcs->m_msh->ele_vector[index];
      val = 0.0;
      for(j=0; j<elem->GetVertexNumber(); j++)    // average all adjoining nodes
      {
         enode = elem->GetNodeIndex(j);
         val +=  cplpcs->GetNodeValue(enode,idx);
      }
      val = val/((double)elem->GetVertexNumber());
   }
   return val;
}


/**************************************************************************
   GeoSys - Function: Get mean element value for element index from secondary node values
                  of process pcs_name and for variable var_name; old and new timelevel
    01/2006   SB    Implementation
**************************************************************************/
double PCSGetEleMeanNodeSecondary(long index, const string &pcs_name,
const string &var_name, int timelevel)
{
   double val = 1.0;                              // As this returns saturation, default is fully saturated = 1.0;
   int idx, j;
   long enode;
   CRFProcess *m_pcs = NULL;

   // Get process by process name
   ProcessType pcs_type (convertProcessType (pcs_name));
   m_pcs = PCSGet(pcs_type);
   if (m_pcs)
   {
      // Get index of secondary node value
      idx = m_pcs->GetNodeValueIndex(var_name) + timelevel;
      // Get element with index index
      CElem* elem = NULL;
      elem = m_pcs->m_msh->ele_vector[index];
      val = 0.0;
                                                  // average all adjoining nodes
      for (j = 0; j < elem->GetVertexNumber(); j++)
      {
         enode = elem->GetNodeIndex(j);
         val += m_pcs->GetNodeValue(enode, idx);
      }
      val = val / ((double) elem->GetVertexNumber());
   }
   return val;
}


/*************************************************************************
GeoSys-FEM Function:
01/2006 OK Implementation
**************************************************************************/
#ifndef NEW_EQS                                   //WW. 07.11.2008
void CRFProcess::SetNODFlux()
{
   long i;
   //----------------------------------------------------------------------
   int nidx;
   nidx = GetNodeValueIndex("FLUX");
   if(nidx<0)
      return;
   double m_val;
   for(i=0;i<(long)m_msh->nod_vector.size();i++)
   {
      m_val = eqs->b[i];                          //? m_nod->eqs_index
      SetNodeValue(i,nidx,m_val);
   }
   //----------------------------------------------------------------------
}
#endif

/*************************************************************************
GeoSys-FEM Function:
01/2006 OK Implementation
**************************************************************************/
#ifndef NEW_EQS                                   //WW. 07.11.2008
void CRFProcess::AssembleParabolicEquationRHSVector()
{
   //OK  long i;
   //----------------------------------------------------------------------
   // Init
   /*	PCH & WW
     for(i=0;i<eqs->dim;i++)
     {
       eqs->b[i] = 0.0;
     }
     //----------------------------------------------------------------------
     CElem* m_ele = NULL;
     for(i=0;i<(long)m_msh->ele_vector.size();i++)
     {
       m_ele = m_msh->ele_vector[i];
       if(m_ele->GetMark()) // Marked for use
   {
   fem->ConfigElement(m_ele,false);
   fem->AssembleParabolicEquationRHSVector();
   //fem->AssembleParabolicEquationLHSMatrix();
   }
   }
   */
   //----------------------------------------------------------------------
}
#endif                                            //#ifndef NEW_EQS //WW. 07.11.2008
/*************************************************************************
GeoSys-FEM Function:
06/2006 YD Implementation
02/2008 JOD removed
03/2008 HS/KG activated for adaptive time step
Reload primary variable
**************************************************************************/
void CRFProcess::PrimaryVariableReload()
{
   //  char pcsT;
   //  pcsT = _pcs_type_name[0];
   //  switch(pcsT){
   //    case 'L':
   //      break;
   //    case 'U':
   //      break;
   //    case 'G':
   //      break;
   //    case 'T':
   //      break;
   //    case 'C':
   //      break;
   //    case 'M':
   //      PrimaryVariableReloadTransport();
   //      break;
   //    case 'R': // Richards flow
   //      PrimaryVariableReloadRichards();
   //      break;
   //  }

   switch (this->getProcessType())
   {
      case MASS_TRANSPORT:
         PrimaryVariableReloadTransport();
         break;
      case RICHARDS_FLOW:
         PrimaryVariableReloadRichards();
         break;
      default:
         break;
   }
}


/*************************************************************************
GeoSys-FEM Function:
06/2006 YD Implementation
02/2008 JOD removed
Reload primary variable of Richards Flow
**************************************************************************/
void CRFProcess::PrimaryVariableReloadRichards()
{
   int i;
   int idxp,idx_storage;
   double storage_p;

   idxp = GetNodeValueIndex(pcs_primary_function_name[0]);
   idx_storage = GetNodeValueIndex("STORAGE_P");
   for(i=0;i<(long)m_msh->GetNodesNumber(false);i++)
   {
      storage_p = GetNodeValue(i,idx_storage);
      SetNodeValue(i,idxp,storage_p);
      SetNodeValue(i,idxp+1,storage_p);
   }
   CalcSecondaryVariables(0);
   CalcSecondaryVariables(1);
}


/*************************************************************************
GeoSys-FEM Function:
12/2007 kg44 Implementation
Reload primary variable for Transport
**************************************************************************/
void CRFProcess::PrimaryVariableReloadTransport()
{
   long i;
   int idxp,idx_storage;
   double conc_back;
   char* mcomp_name;

                                                  //kg44 test
   idxp = GetNodeValueIndex(pcs_primary_function_name[0]);

   mcomp_name = new char[80];
   sprintf(mcomp_name, "%s%i","CONC_BACK_",pcs_component_number);

   idx_storage = GetNodeValueIndex(mcomp_name);

   for(i=0;i<(long)m_msh->GetNodesNumber(false);i++)
   {
      conc_back = GetNodeValue(i,idx_storage);
      SetNodeValue(i,idxp,conc_back);
      SetNodeValue(i,idxp+1,conc_back);

   }

   CalcSecondaryVariables(0);
   CalcSecondaryVariables(1);
}


/*************************************************************************
GeoSys-FEM Function:
12/2007 kg44 Implementation
Reload primary variable for Transport
**************************************************************************/
void CRFProcess::PrimaryVariableStorageTransport()
{
   long i;
   int idxp,idx_storage;
   double concentration;
   char* mcomp_name;

   idxp = GetNodeValueIndex(pcs_primary_function_name[0]);
   mcomp_name= new char[80];
   sprintf(mcomp_name, "%s%i","CONC_BACK_",pcs_component_number);
   //   cout << "mcomp_name"<< mcomp_name<<endl;
   idx_storage = GetNodeValueIndex(mcomp_name);

   for(i=0;i<(long)m_msh->GetNodesNumber(false);i++)
   {
      concentration = GetNodeValue(i,idxp);
      SetNodeValue(i,idx_storage,concentration);
      //    SetNodeValue(i,idx_storage+1,concentration);
   }
}


/*************************************************************************
GeoSys-FEM Function:
06/2006 YD Implementation
Reload primary variable of Richards Flow
**************************************************************************/
void CRFProcess::PrimaryVariableStorageRichards()
{
   int i;
   int idxp,idx_storage;
   double pressure;

   idxp = GetNodeValueIndex(pcs_primary_function_name[0])+1;
   idx_storage = GetNodeValueIndex("STORAGE_P");
   for(i=0;i<(long)m_msh->GetNodesNumber(false);i++)
   {
      pressure = GetNodeValue(i,idxp);
      SetNodeValue(i,idx_storage,pressure);
      SetNodeValue(i,idx_storage+1,pressure);
   }
}


//*************************************************************************
//GeoSys-FEM Function:
//12/2007 kg44 Implementation
//check change of concentration and set new time step factor
//**************************************************************************/
#ifdef kg44                                       // WW
double CRFProcess::GetNewTimeStepSizeTransport(double mchange)
{
   long i, mnode=-1;
   int comp;
   int idxn,idxo;
   double conc_new, conc_old,/*time_coeff,*/ max_change=1.0e-10, tchange=1.0;
   char* mcomp_name;

   idxo = GetNodeValueIndex(pcs_primary_function_name[0]);
   comp = pcs_component_number;                   // get component number
   //   cout << "comp number "<<comp<<endl;
   mcomp_name= new char[80];
   sprintf(mcomp_name, "%s%li","CONC_BACK_",comp);
   //   cout << "mcomp_name"<< mcomp_name<<endl;
   idxn = GetNodeValueIndex(mcomp_name);
   for(i=0;i<(long)m_msh->GetNodesNumber(false);i++)
   {
      conc_old = abs(GetNodeValue(i,idxo));
      conc_new = abs(GetNodeValue(i,idxn));
      if (((conc_old) > MKleinsteZahl)&&((conc_new) > MKleinsteZahl))
      {
         max_change=MMax(max_change,conc_new/conc_old); mnode=i;
      }
   }
   tchange=mchange/max_change;
   if (tchange > 2.0) tchange=2.0;
   cout << "Transport: max change of "<< max_change << " at node " << mnode << " factor " << tchange<< endl;
   return tchange;
}
#endif

/**************************************************************************
FEMLib-Method:
11/2005 MB Implementation
03/2006 OK 2nd version (binary coupling)
**************************************************************************/
#ifndef NEW_EQS                                   //WW. 07.11.2008
void CRFProcess::SetCPL()
{
   int i;
   double value = 0.0;
   //----------------------------------------------------------------------
   // Nothing to do
   if ((cpl_type_name.size() == 0) || (cpl_type_name.compare("PARTITIONED")
      == 0))
      return;
   //----------------------------------------------------------------------
   // PCS CPL
   CRFProcess* m_pcs_cpl = PCSGet(cpl_type_name);
   if (!m_pcs_cpl)
   {
      cout << "Fatal error in CRFProcess::SetCPL: no PCS data" << endl;
      return;
   }
   //----------------------------------------------------------------------
   // MSH data for PCS CPL
   CFEMesh* m_msh_cpl = m_pcs_cpl->m_msh;
   if (!m_msh_cpl)
   {
      cout << "Fatal error in CRFProcess::SetCPL: no MSH data" << endl;
      return;
   }
   //----------------------------------------------------------------------
   // GEO data for PCS CPL
   Surface* m_sfc = GEOGetSFCByName(m_msh_cpl->geo_name);
   if (!m_sfc)
   {
      cout << "Fatal error in CRFProcess::SetCPL: no GEO data" << endl;
      return;
   }
   //----------------------------------------------------------------------
   //......................................................................
   // MSH nodes of PCS CPL
   cout << "      ->CPL: " << cpl_type_name << ": ";
   vector<long> cpl_msh_nodes_vector;
   m_msh_cpl->GetNODOnSFC(m_sfc, cpl_msh_nodes_vector);
   if ((int) cpl_msh_nodes_vector.size() == 0)
   {
      cout << "Warning in CRFProcess::SetCPL: no MSH nodes found" << endl;
   }
   cout << "CPL nodes = " << (int) cpl_msh_nodes_vector.size() << endl;
   //.....................................................................-
   // MSH nodes of PCS
   cout << "      ->CPL: " << convertProcessTypeToString(this->getProcessType()) << ": ";
   vector<long> msh_nodes_vector;
   m_msh->GetNODOnSFC(m_sfc, msh_nodes_vector);
   if ((int) msh_nodes_vector.size() == 0)
   {
      cout << "Warning in CRFProcess::SetCPL: no MSH nodes found" << endl;
   }
   cout << "CPL nodes = " << (int) msh_nodes_vector.size() << endl;
   //----------------------------------------------------------------------
   if (m_msh_cpl->pcs_name.compare("RICHARDS_FLOW") == 0)
   {
      m_msh->SetNODPatchAreas();
      int nidx = GetNodeValueIndex("WDEPTH");
      long st_node_number;
      double st_node_value = 0.0;
      CNode* m_nod = NULL;
      for (i = 0; i < (int) msh_nodes_vector.size(); i++)
      {
         value = -2.314e-02;
         st_node_number = msh_nodes_vector[i];
         m_nod = m_msh->nod_vector[st_node_number];
         st_node_value = GetNodeValue(st_node_number, nidx);
         st_node_value /= m_nod->patch_area;
         value *= st_node_value;
         //cout << "CPL value = " << value << endl;
         eqs->b[st_node_number] += value;
      }
   }

   //	if (_pcs_type_name.compare("RICHARDS_FLOW") == 0
   //				&& m_msh_cpl->pcs_name.compare("OVERLAND_FLOW") == 0) { // ToDo
   if (this->getProcessType() == RICHARDS_FLOW
                                                  // ToDo
      && m_msh_cpl->pcs_name.compare("OVERLAND_FLOW") == 0)
   {
      long msh_node_number;
      long cpl_msh_nod_number;
      long cpl_msh_ele_number;
      value = 0.0;
      cout << "CPL value = " << value << endl;
      // PCS-CON
      CRFProcess* m_pcs_cond = PCSGet(cpl_type_name);
      //int nidx = m_pcs_cond->GetNodeValueIndex(m_pcs_cond->pcs_primary_function_name[0]);
      int cpl_nidx = m_pcs_cond->GetNodeValueIndex("WDEPTH");
      //----------------------------------------------------------------------
      // ELE of PCS_CPL related to NOD of PCS
      //  CFEMesh* m_msh_this = MSHGet("RICHARDS_FLOW_LOCAL");
      CElem* m_ele_cnd = NULL;
      //----------------------------------------------------------------------
      //  CSourceTermGroup *m_st_group = NULL;
      //  CSourceTerm *m_st = NULL;
      //  m_st_group = STGetGroup(_pcs_type_name,pcs_primary_function_name[0]);
      //----------------------------------------------------------------------
      double cpl_ele_val = 0.0;
      int j;
      CNodeValue *cnodev = NULL;
      //  for(i=0;i<(int)m_st_group->group_vector.size();i++)
      //  ofstream st_out_file("st_out_file.txt",ios::app);
      for (i = 0; i < (int) st_node_value.size(); i++)
      {
         cnodev = st_node_value[i];
         // MSH-PCS
         //m_nod = m_msh_this->nod_vector[m_st_group->group_vector[i]->msh_node_number];
                                                  //m_st_group->group_vector[i]->msh_node_number; //0
         msh_node_number = cnodev->msh_node_number;
         // MSH-PCS-CPL
         cpl_msh_nod_number = msh_node_number;
         cpl_msh_ele_number = pcs_number;         //OK:TODO
         m_ele_cnd = m_pcs_cond->m_msh->ele_vector[cpl_msh_ele_number];
         for (j = 0; j < m_ele_cnd->GetNodesNumber(false); j++)
         {
            cpl_ele_val += m_pcs_cond->GetNodeValue(
               m_ele_cnd->nodes_index[j], cpl_nidx);
         }
         cpl_ele_val /= m_ele_cnd->GetNodesNumber(false);
         // VAL-CON
         value = 2.314e-02 * cpl_ele_val * 1e-2;
         //    st_out_file << value << endl;
         // EQS-RHS
         eqs->b[msh_node_number] += value;
      }
      //----------------------------------------------------------------------
      /*
       CNodeValue* m_node_value = NULL;
       m_node_value = new CNodeValue();
       m_node_value->msh_node_number = msh_node_number;
       m_node_value->geo_node_number = m_pnt->id;
       m_node_value->node_value = value;
       CSourceTermGroup *m_st_group = NULL;
       m_st_group = STGetGroup(_pcs_type_name,pcs_primary_function_name[0]);
       m_st_group->group_vector.push_back(m_node_value);
       m_st_group->st_group_vector.push_back(m_st); //OK
       */
   }

   //	if (_pcs_type_name.compare("GROUNDWATER_FLOW") == 0
   //				&& m_msh_cpl->pcs_name.compare("OVERLAND_FLOW") == 0) { // ToDo
   if (this->getProcessType() == GROUNDWATER_FLOW
                                                  // ToDo
      && m_msh_cpl->pcs_name.compare("OVERLAND_FLOW") == 0)
   {
      long ie = (long) msh_nodes_vector.size()
         / (m_msh->getNumberOfMeshLayers() + 1);
      long of_node_number, gf_node_number;
      double of_node_value, gf_node_value;
      //  CNode* m_nod = NULL;
      int of_nidx = GetNodeValueIndex("WDEPTH");

      for (i = 0; i < ie; i++)
      {
         of_node_number = msh_nodes_vector[i];    //ToDo
         of_node_value = m_pcs_cpl->GetNodeValue(of_node_number, of_nidx);
         //m_nod = m_msh->nod_vector[gf_node_number];
         //st_node_value /= m_nod->patch_area;
         gf_node_value = of_node_value * 2e-11;
         if (gf_node_value > 1e-13)
            cout << "CPL value = " << gf_node_value << endl;
         gf_node_number = msh_nodes_vector[i];
         eqs->b[gf_node_number] += gf_node_value;
      }
   }

   int idx = GetNodeValueIndex("FLUX") + 1;
   //for(i=0;i<(int)msh_nodes_vector.size();i++)
   for (i = 0; i < 1; i++)
   {
      SetNodeValue(msh_nodes_vector[i], idx, value);
   }
}
#endif

/**************************************************************************
PCSLib-Method:
04/2006 OK Implementation
**************************************************************************/
void CRFProcess::CreateBCGroup()
{
   cout << "->Create BC" << '\n';
   CBoundaryConditionsGroup *m_bc_group = NULL;

   std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));

   for (size_t i = 0; i < GetPrimaryVNumber(); i++)
   {
      BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
      m_bc_group = new CBoundaryConditionsGroup();
                                                  //OK
      m_bc_group->setProcessTypeName(pcs_type_name);
                                                  //OK
      m_bc_group->setProcessPrimaryVariableName(pcs_primary_function_name[i]);
      m_bc_group->Set(this, Shift[i]);
      bc_group_list.push_back(m_bc_group);
   }
}


/**************************************************************************
PCSLib-Method:
04/2006 OK Implementation
**************************************************************************/
void CRFProcess::CreateSTGroup()
{
   cout << "->Create ST" << '\n';
   CSourceTermGroup *m_st_group = NULL;
   // WW
   std::ifstream *iSourceNBC_RHS_file = NULL;
   std::ofstream *oSourceNBC_RHS_file = NULL;

   std::string pcs_type_name (convertProcessTypeToString (this->getProcessType()));

   if (WriteSourceNBC_RHS == 1)
   {
      string m_file_name = FileName + "_" + pcs_type_name
         + "_source_Neumann_RHS.bin";
      oSourceNBC_RHS_file = new ofstream(m_file_name.c_str(), ios::binary);
   }
   else if (WriteSourceNBC_RHS == 2)
   {
      string m_file_name = FileName + "_" + pcs_type_name
         + "_source_Neumann_RHS.bin";
      iSourceNBC_RHS_file = new ifstream(m_file_name.c_str(), ios::binary);
      if (!iSourceNBC_RHS_file->good())
         cout << "_source_Neumann_RHS file is not found" << endl;
   }

   for (size_t i = 0; i < GetPrimaryVNumber(); i++)
   {
      //OK m_st_group = m_st_group->Get(pcs_primary_function_name[i]);
      m_st_group = STGetGroup(pcs_type_name, pcs_primary_function_name[i]);
      if (!m_st_group)
      {
         m_st_group = new CSourceTermGroup();
                                                  //OK
         m_st_group->pcs_type_name = pcs_type_name;
                                                  //OK
         m_st_group->pcs_pv_name = pcs_primary_function_name[i];
         //      if(iSourceNBC_RHS_file)  // Read from data. WW
         //        m_st_group->Read(*iSourceNBC_RHS_file);
         //	  else
         m_st_group->Set(this, Shift[i]);
         st_group_list.push_back(m_st_group);
      }
   }
   if (oSourceNBC_RHS_file)                       // WW
      //    WriteRHS_of_ST_NeumannBC(*oSourceNBC_RHS_file);

      if (iSourceNBC_RHS_file)                    //WW
   {
      iSourceNBC_RHS_file->close();
      delete iSourceNBC_RHS_file;
      iSourceNBC_RHS_file = NULL;
   }
   if (oSourceNBC_RHS_file)                       //WW
   {
      oSourceNBC_RHS_file->close();
      delete oSourceNBC_RHS_file;
      oSourceNBC_RHS_file = NULL;
   }
}

/**************************************************************************
    PCSLib-Method:
    08/2006 OK Implementation
 **************************************************************************/
double CRFProcess::CalcELEFluxes(const GEOLIB::Polyline* const ply)
{
   CRFProcess* m_pcs_flow = NULL;
   //	if (_pcs_type_name.find("FLOW") != string::npos) {
   if (isFlowProcess (this->getProcessType()))
   {
      m_pcs_flow = this;
   }
   else
   {
      m_pcs_flow = PCSGet(GROUNDWATER_FLOW);
   }

   int v_eidx[3];
   v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
   v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
   v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");

   for (size_t i = 0; i < 3; i++)
   {
      if (v_eidx[i] < 0)
      {
         std::cout << i << " " << v_eidx[i] << "Velocity output is not specified"
            << std::endl;
         return 0.0;
      }
   }

   // Get elements at GEO
   std::vector<size_t> ele_vector_at_geo;
   m_msh->GetELEOnPLY(ply, ele_vector_at_geo);

   CElem* m_ele = NULL;
   CEdge* m_edg = NULL;
   vec<CEdge*> ele_edges_vector(15);
   double edg_normal_vector[3];
   double vn;
   double edg_length = 0.0;
   double vn_vec[3];
   double edge_vector[3];
   double f_n_sum = 0.0;
   double C_ele = 0.0;
   vec<long> element_nodes(20);
   double f[3], v[3];

   for (size_t i = 0; i < ele_vector_at_geo.size(); i++)
   {
      m_ele = m_msh->ele_vector[ele_vector_at_geo[i]];
      m_ele->SetNormalVector();
      m_ele->GetNodeIndeces(element_nodes);

      // velocity vector
      for (size_t j = 0; j < 3; j++)
      {
         v[j] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[j]);
      }

      // edge projection // edge marked
      m_ele->GetEdges(ele_edges_vector);
      for (size_t j = 0; j < static_cast<size_t>(m_ele->GetEdgesNumber()); j++)
      {
         m_edg = ele_edges_vector[j];
         if (m_edg->GetMark())
         {
            m_edg->SetNormalVector(m_ele->normal_vector, edg_normal_vector);
            edg_length = m_edg->Length();
            m_edg->GetEdgeVector(edge_vector);
         }
      }
      vn = MSkalarprodukt(v, edg_normal_vector, 3);
      for (size_t j = 0; j < 3; j++)
      {
         vn_vec[j] = vn * edg_normal_vector[j];
      }

      switch (this->getProcessType())
      {
         case GROUNDWATER_FLOW:
            // Volume flux = v_n * l^e * z^e
            for (size_t j = 0; j < 3; j++)
            {
               f[j] = vn_vec[j] * edg_length * m_ele->GetFluxArea();
            }
            break;
         case MASS_TRANSPORT:
            // Mass flux = v_n * l^e * z^e * C^e
            C_ele = 0.0;
            for (int j = 0; j < m_ele->GetNodesNumber(false); j++)
            {
               size_t nidx = GetNodeValueIndex(pcs_primary_function_name[0]) + 1;
               C_ele = +GetNodeValue(element_nodes[j], nidx);
            }
            C_ele /= (double) m_ele->GetNodesNumber(false);
            for (size_t j = 0; j < 3; j++)
            {
               f[j] = vn_vec[j] * edg_length * m_ele->GetFluxArea() * C_ele;
            }
            break;
         default:
            break;
      }

      // set
      f_n_sum += MBtrgVec(f, 3);
   }

   ele_vector_at_geo.clear();

   return f_n_sum;
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2004 OK Implementation
 08/2006 OK new
 **************************************************************************/
void CRFProcess::CalcELEVelocities(void)
{
   // If not FLUID_MOMENTUM,
   //	if (_pcs_type_name.compare("RANDOM_WALK") != 0) {
   if (this->getProcessType() != RANDOM_WALK)
   {
      int eidx[3];
      eidx[0] = GetElementValueIndex("VELOCITY1_X");
      eidx[1] = GetElementValueIndex("VELOCITY1_Y");
      eidx[2] = GetElementValueIndex("VELOCITY1_Z");
      for (size_t i = 0; i < 3; i++)
      {
         if (eidx[i] < 0)
         {
            cout << "Fatal error in CRFProcess::CalcELEVelocities - abort"
               << endl;                           //abort();	// PCH commented abort() out for FM.
         }
      }

      FiniteElement::ElementValue* gp_ele = NULL;
      double vx, vy, vz;
      const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
      for (size_t i = 0; i < mesh_ele_vector_size; i++)
      {
         gp_ele = ele_gp_value[i];
         vx = gp_ele->Velocity(0, 0);
         SetElementValue(i, eidx[0], vx);
         SetElementValue(i, eidx[0] + 1, vx);
         vy = gp_ele->Velocity(1, 0);
         SetElementValue(i, eidx[1], vy);
         SetElementValue(i, eidx[1] + 1, vy);
         vz = gp_ele->Velocity(2, 0);
         SetElementValue(i, eidx[2], vz);
         SetElementValue(i, eidx[2] + 1, vz);
      }
   }
}


/*************************************************************************
GeoSys-FEM Function:
08/2006 OK Implementation
**************************************************************************/
#ifndef NEW_EQS                                   //WW. 07.11.2008
                                                  //(vector<long>&ele_number_vector)
void CRFProcess::AssembleParabolicEquationRHSVector(CNode*m_nod)
{
   //cout << "CRFProcess::AssembleParabolicEquationRHSVector" << endl;
   int i;
   long ldummy;
   double ddummy;
   //----------------------------------------------------------------------
   // Init
   for(i=0;i<(int)m_nod->connected_elements.size();i++)
   {
      eqs->b[m_nod->connected_elements[i]] = 0.0;
   }
   //----------------------------------------------------------------------
   CElem* m_ele = NULL;
   CEdge* m_edg = NULL;
   double edg_normal_vector[3];
   double edge_mid_point[3];
   vec<CEdge*>ele_edges_vector(15);
   int j;
   double aux_vector[3];
   double* gravity_center = NULL;
   double check_sign;
   //----------------------------------------------------------------------
   // Element velocity
   int v_eidx[3];
   //kg44  v_eidx[0] = GetElementValueIndex("VELOCITY1_X");
   //kg44  v_eidx[1] = GetElementValueIndex("VELOCITY1_Y");
   //kg44  v_eidx[2] = GetElementValueIndex("VELOCITY1_Z");
   CRFProcess* m_pcs_flow = NULL;
   //  if(_pcs_type_name.find("FLOW")!=string::npos) {
   if(isFlowProcess(this->getProcessType()))
   {
      m_pcs_flow = this;
   }
   else
   {
      m_pcs_flow = PCSGet(GROUNDWATER_FLOW);
   }
   v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
   v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
   v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");
   for(i=0;i<3;i++)
   {
      if(v_eidx[i]<0)
      {
         cout << v_eidx[i] << i << " Warning in CRFProcess::AssembleParabolicEquationRHSVector - no PCS-VEL data" << endl;
         return;
      }
   }
   double v[3];
   //======================================================================
   // Topology
   for(i=0;i<(int)m_nod->connected_elements.size();i++)
   {
      m_ele = m_msh->ele_vector[m_nod->connected_elements[i]];
      m_ele->SetNormalVector();                   //OK_BUGFIX
      v[0] = m_pcs_flow->GetElementValue(m_ele->GetIndex(),v_eidx[0]);
      v[1] = m_pcs_flow->GetElementValue(m_ele->GetIndex(),v_eidx[1]);
      v[2] = m_pcs_flow->GetElementValue(m_ele->GetIndex(),v_eidx[2]);
      m_ele->SetMark(false);
      switch(m_ele->GetElementType())
      {
         //------------------------------------------------------------------
         // line elements
         case MshElemType::LINE:
            v[1] = GetElementValue(m_ele->GetIndex(),v_eidx[0]);
            v[0] = GetElementValue(m_ele->GetIndex(),v_eidx[1]);
            if(m_nod->connected_elements.size()==1)
            {
               m_ele->SetMark(true);
               break;
            }
            gravity_center = m_ele->GetGravityCenter();
            aux_vector[0] = gravity_center[0] - m_nod->X();
            aux_vector[1] = gravity_center[1] - m_nod->Y();
            aux_vector[2] = gravity_center[2] - m_nod->Z();
            check_sign = MSkalarprodukt(v,aux_vector,3);
            if(check_sign<0.0)
               m_ele->SetMark(true);
            break;
            //------------------------------------------------------------------
            // tri elements
         case MshElemType::TRIANGLE:
            m_ele->GetEdges(ele_edges_vector);
            for(j=0;j<(int)m_ele->GetEdgesNumber();j++)
            {
               m_edg = ele_edges_vector[j];
               if(m_edg->GetMark())
               {
                  m_edg->SetNormalVector(m_ele->normal_vector,edg_normal_vector);
                  break;
                  /*
                              m_edg->GetEdgeMidPoint(edge_mid_point);
                             gravity_center = m_ele->GetGravityCenter();
                             aux_vector[0] = gravity_center[0] - edge_mid_point[0];
                             aux_vector[1] = gravity_center[1] - edge_mid_point[1];
                             aux_vector[2] = gravity_center[2] - edge_mid_point[2];
                             check_sign = MSkalarprodukt(edg_normal_vector,aux_vector,3);
                             if(check_sign<0.0) break;
                  */
               }
            }
            if(m_edg->GetMark()) break;
            //----------------------------------------------------------------
            // ToDo
         default:
            cout << "Warning in CRFProcess::AssembleParabolicEquationRHSVector - not implemented for this element type" << endl;
            break;
      }                                           // switch
   }
   //======================================================================
   for(i=0;i<(int)m_nod->connected_elements.size();i++)
   {
      m_ele = m_msh->ele_vector[m_nod->connected_elements[i]];
      switch(m_ele->GetElementType())
      {
         //------------------------------------------------------------------
         // line elements
         case MshElemType::LINE:
            if(m_ele->GetMark())
            {
               cout << m_ele->GetIndex() << endl;
               ldummy = m_nod->GetIndex();
               ddummy = eqs->b[m_nod->GetIndex()];
               fem->ConfigElement(m_ele,false);
               fem->AssembleParabolicEquationRHSVector();
               ddummy = eqs->b[m_nod->GetIndex()];
            }
            break;
            //------------------------------------------------------------------
            // tri elements
         case MshElemType::TRIANGLE:
            m_edg->GetEdgeMidPoint(edge_mid_point);
            gravity_center = m_ele->GetGravityCenter();
            aux_vector[0] = gravity_center[0] - edge_mid_point[0];
            aux_vector[1] = gravity_center[1] - edge_mid_point[1];
            aux_vector[2] = gravity_center[2] - edge_mid_point[2];
            check_sign = MSkalarprodukt(edg_normal_vector,aux_vector,3);
            if(check_sign<0.0) continue;
            {
               //cout << m_ele->GetIndex() << endl;
               fem->ConfigElement(m_ele,false);
               fem->AssembleParabolicEquationRHSVector();
            }
            break;
            //----------------------------------------------------------------
            // ToDo
         default:
            cout << "Warning in CRFProcess::AssembleParabolicEquationRHSVector - not implemented for this element type" << endl;
            break;
      }                                           // switch
   }
   //======================================================================
}
#endif
/**************************************************************************
PCSLib-Method:
08/2006 OK Implementation
compare with CMCDs PCSGetFluxProcess
**************************************************************************/
CRFProcess* PCSGetFlow()
{
   CRFProcess *m_pcs = NULL;
   for(int i=0;i<(int)pcs_vector.size();i++)
   {
      m_pcs = pcs_vector[i];
      //    if(m_pcs->_pcs_type_name.find("FLOW")!=string::npos)
      if(isFlowProcess (m_pcs->getProcessType()))
      {
         return m_pcs;
      }
   }
   return NULL;
}


/**************************************************************************
PCSLib-Method:
based on MMPCalcSecondaryVariables
01/2007 OK Implementation
**************************************************************************/

/**************************************************************************
PCSLib-Method:
based on MMPCalcSecondaryVariables
01/2007 OK Implementation
**************************************************************************/
void CRFProcess::SetBC()
{
   CBoundaryCondition *m_bc = NULL;
   CBoundaryConditionNode *m_node = NULL;
   int nidx = GetNodeValueIndex(pcs_primary_function_name[0]);
   for(long i=0;i<(long)bc_node_value.size();i++)
   {
      m_node = bc_node_value[i];
      m_bc = bc_node[i];
                                                  // old time
      SetNodeValue(m_node->msh_node_number,nidx,m_node->node_value);
                                                  // new time
      SetNodeValue(m_node->msh_node_number,nidx+1,m_node->node_value);
   }
}


/**************************************************************************
Task: Postprocessing function calculates the NAPL saturation after
      Flow, Transport and kkinetic NAPL dissolution

Programing:
   01/2008   CB   Implementation                                          */
/**************************************************************************/
void CalcNewNAPLSat(CRFProcess*m_pcs)
{

   long i, j, k, l;
   int idx0, idx1, idx2, idxC;
   long nnodes, nNAPLcomps;
   double conc, rho_N_new, rho_N_old, rho_N_fluid;
   double mass, volume;
   double satu_N_new, satu_N_old ;

   vector<int>pcs_napl_comps_vector;
   vector<double>molar_weights_vector;
   vector<double>molar_densities_vector;

   string var_name;
   int no_processes;

   CRFProcess *n_pcs = NULL;
   CFluidProperties* m_mfp = NULL;

   nnodes = (long) fem_msh_vector[0]->nod_vector.size();

   // make sure m_pcs is NAPL-phase
   if(m_pcs->pcs_type_number==0)
      m_pcs = pcs_vector[m_pcs->pcs_number+1];    // this is the NAPL phase

   m_mfp = mfp_vector[m_pcs->pcs_type_number];
   //m_mfp->mode = 1; // CB ??

   // Get indices of node value variables for phase 2
   idx0 = m_pcs->GetNodeValueIndex("SATURATION2");// old timelevel
   idx1 = m_pcs->GetNodeValueIndex("DENSITY2");
   idx2 = m_pcs->GetNodeValueIndex("SATURATION1");// old timelevel

   i = j = k = l = 0;
   no_processes =(int)pcs_vector.size();

   for(i=0;i<no_processes;i++)
   {
      n_pcs = pcs_vector[i];
      //     if(n_pcs->_pcs_type_name.compare("MASS_TRANSPORT")==0){ // TF
      if(n_pcs->getProcessType () == MASS_TRANSPORT)
      {
         j = n_pcs->GetProcessComponentNumber();
         if(cp_vec[j]->transport_phase==3)        // is in napl
         {
            pcs_napl_comps_vector.push_back(i);   // store processes that are NAPL dissolution components
                                                  // get the corresponding molar weight
            molar_weights_vector.push_back(cp_vec[j]->molar_weight);
                                                  // get the corresponding densities
            molar_densities_vector.push_back(cp_vec[j]->molar_density);
            l++;
         }
      }
   }
   nNAPLcomps = l;

   for(i = 0; i < nnodes; i++)
   {
      conc = rho_N_new = rho_N_old = rho_N_fluid = satu_N_new = satu_N_old = mass = volume = 0;

      // deteremine the old NAPL DENSITY after flow / transport step
      for(j = 0; j < nNAPLcomps; j++)
      {
         l=pcs_napl_comps_vector[j];
         idxC = pcs_vector[l]->GetNodeValueIndex(pcs_vector[l]->pcs_primary_function_name[0]);
                                                  // old timelevel
         conc = pcs_vector[l]->GetNodeValue(i, idxC);
         mass += conc * molar_weights_vector[j];
      }
      if(mass > 0)
         rho_N_old = mass;                        // [kg/m?3REV]
      else
         rho_N_old = 0;
      // determine the new NAPL density rho_N_neu at current node
      mass = 0;
      for(j = 0; j < nNAPLcomps; j++)
      {
         l=pcs_napl_comps_vector[j];
         idxC = pcs_vector[l]->GetNodeValueIndex(pcs_vector[l]->pcs_primary_function_name[0]);
                                                  // +1 new timelevel
         conc = pcs_vector[l]->GetNodeValue(i, idxC+1);
         if(fabs(conc)<1e-19)
            conc = 0.0;
         mass += conc * molar_weights_vector[j];
                                                  // this is required calculating the napl fluid density
         volume += conc / molar_densities_vector[j];
      }
      if (mass > 0)
         rho_N_new = mass;                        // [kg/m?3REV]
      else
         rho_N_new = 0;

      // get the old SATURATION2 of NAPL after flow / transport step
      satu_N_old = m_pcs->GetNodeValue(i, idx0);
      // calculate new NAPL Saturation: dSatu = Satu(t+dt)-Satu(t) = Satu(t)*(1-rho(t)/rho(t+dt))
      if(satu_N_old * rho_N_new * rho_N_old > 0)
                                                  //fast
         satu_N_new = satu_N_old + satu_N_old * (1 - rho_N_old / rho_N_new);
      //satu_N_new = satu_N_old + satu_N_old * (rho_N_new / rho_N_old - 1);   //slow
      else
         satu_N_new = satu_N_old ;
      if (satu_N_new <0)
      {
         cout << " Warning in fct CalcNewNAPLSat: NAPL-Sat: " << satu_N_new << endl;
         satu_N_new = MRange(0.0, satu_N_new, 1.0);
      }
      // set new SATURATION2
      //if(satu_N_new > 0.00001)  satu_N_new = 1-0.95; // 0.985 CB 14.01.09 Hansen+Kueper BM

      //m_pcs->SetNodeValue(i, idx0, satu_N_new);  // idx0 for timelevel 0 ??
      m_pcs->SetNodeValue(i, idx0+1, satu_N_new); // idx0+1 for timelevel 1 ??
      // set new SATURATION1
      //m_pcs->SetNodeValue(i, idx2, 1-satu_N_new);     // idx2 for timelevel 0 ??
                                                  // idx2+1 for timelevel 1 ??
      m_pcs->SetNodeValue(i, idx2+1, 1-satu_N_new);

      // finally determine the new napl fluid density
      if(mass * volume > 0)
         rho_N_fluid = mass / volume;             // [kg/m?3N] = [kg/m?3REV] / [m?3N/m?3REV]
      else
         rho_N_fluid = m_mfp->Density();          // use NAPL phase fluid density as defined in .mfp
      // set new DENSITY2
      m_pcs->SetNodeValue(i, idx1, rho_N_fluid);

   }

   pcs_napl_comps_vector.clear();
   molar_weights_vector.clear();
   molar_densities_vector.clear();

   //m_pcs->WriteAllVariables();
}


/**************************************************************************
Task: Calculates the NAPL and Water phase saturations
      in case of NAPL dissolution
Programing:
   01/2008   CB   Implementation                                          */
/**************************************************************************/
double CalcNAPLDens(CRFProcess*m_pcs, int node)
{
   long i, j, k, l;
   int idx1, idxC;
   int nNAPLcomps;
   double conc, rho_N_new, mass, volume;
   string var_name;

   vector<int>pcs_napl_comps_vector;
   vector<double>molar_weights_vector;
   vector<double>molar_densities_vector;

   int no_processes;
   CRFProcess *n_pcs = NULL;
   CFluidProperties* m_mfp = NULL;

   // Get indices of node value variable DENSITY2 for NAPL phase
   if(m_pcs->pcs_type_number==0)
      m_pcs = pcs_vector[m_pcs->pcs_number+1];
   idx1 = m_pcs->GetNodeValueIndex("DENSITY2");

   m_mfp = mfp_vector[m_pcs->pcs_type_number];
   //m_mfp->mode = 1; // CB ??

   i = j = k = l = 0;
   no_processes =(int)pcs_vector.size();

   //collect parameters and concentrations for NAPL-components
   for (i = 0; i < no_processes; i++)
   {
      n_pcs = pcs_vector[i];
      //	  if(n_pcs->_pcs_type_name.compare("MASS_TRANSPORT")==0){
      if (n_pcs->getProcessType () == MASS_TRANSPORT)
      {
         j = n_pcs->GetProcessComponentNumber();
         if (cp_vec[j]->transport_phase == 3)     // is in NAPL component
         {
            pcs_napl_comps_vector.push_back(i);   // store processes that are NAPL dissolution components
                                                  // get the corresponding molar weight
            molar_weights_vector.push_back(cp_vec[j]->molar_weight);
            molar_densities_vector.push_back(
               cp_vec[j]->molar_density);         // get the corresponding densities
            l++;
         }
      }
   }
   nNAPLcomps = l;

   // determine the NAPL density rho_N_neu at current node
   conc = rho_N_new = mass = volume = 0;
   for(j = 0; j < nNAPLcomps; j++)
   {
      l=pcs_napl_comps_vector[j];
      idxC = pcs_vector[l]->GetNodeValueIndex(pcs_vector[l]->pcs_primary_function_name[0]);
      conc = pcs_vector[l]->GetNodeValue(node, idxC);
      if(fabs(conc)<1e-19)
         conc = 0.0;
      mass += conc * molar_weights_vector[j];     // [kg/m3REV] = [molN/m3REV] * [kg/molN]
      volume += conc / molar_densities_vector[j]; // [m3N/m3REV] = [molN/m3REV] / [molN/m3REV]
   }
   if(mass * volume > 0)
      rho_N_new = mass / volume;                  // [kg/m3N] = [kg/m3REV] / [m3N/m3REV]
   else
      rho_N_new = m_mfp->Density();               // use NAPL phase fluid density as defined in .mfp

   pcs_napl_comps_vector.clear();
   molar_weights_vector.clear();
   molar_densities_vector.clear();

   return rho_N_new;
}


/**************************************************************************
Task: Preprocessing function set flag pcs->flow_pcs_type for
      calculation of velocity and Saturation in mass transport element
      matrices

Programing:
   01/2008   CB   Implementation                                          */
/**************************************************************************/
void SetFlowProcessType()
{
   int i;
   int no_processes, flowtype;
   CRFProcess *m_pcs = NULL;

   m_pcs = PCSGetFlow();
   flowtype = m_pcs->type;
   no_processes =(int)pcs_vector.size();

   for(i=0;i<no_processes;i++)
   {
      m_pcs = pcs_vector[i];
      //if(m_pcs->_pcs_type_name.compare("MASS_TRANSPORT")==0)
      //m_pcs->twophaseflow=true;
      m_pcs->flow_pcs_type = flowtype;
   }
}


/**************************************************************************
Task: Postprocessing function copies the new time step node values of
      secondary variables PRESSURE1 and SATURATION1 to the old time level

Programing:
   13/2008   CB   Implementation                                          */
/**************************************************************************/
void CopyTimestepNODValuesSVTPhF()
{
   long i, j;
   int idx0, idx1;
   long nnodes;
   double val = 0;
   CRFProcess *m_pcs = NULL;

   nnodes = (long) fem_msh_vector[0]->nod_vector.size();

   for (j = 0; j < 2; j++)                        //pcs 1 and 2
   {
      m_pcs = pcs_vector[j];
      //		if (m_pcs->_pcs_type_name.compare("TWO_PHASE_FLOW") != 0)
      if (m_pcs->getProcessType () == TWO_PHASE_FLOW)
         break;
      if (j == 0)
                                                  // old timelevel
         idx0 = m_pcs->GetNodeValueIndex("PRESSURE2");
      else
                                                  // old timelevel
         idx0 = m_pcs->GetNodeValueIndex("SATURATION1");
      idx1 = idx0 + 1;
      for (i = 0; i < nnodes; i++)
      {
         val = m_pcs->GetNodeValue(i, idx1);
         m_pcs->SetNodeValue(i, idx0, val);
      }
      //m_pcs->WriteAllVariables();
   }
}


/**************************************************************************
PCSLib-Method:
based on WriteSolution by WW
01/2007 OK Implementation
**************************************************************************/
void CRFProcess::WriteAllVariables()
{
   std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));
   string m_file_name = FileName + "_" + pcs_type_name + "_"
      + pcs_primary_function_name[0] + ".asc";
   ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
   if (!os.good())
   {
      cout << "Failure to open file: " << m_file_name << endl;
      abort();
   }
   //
   int j;
   int idx[20];
   for (j = 0; j < pcs_number_of_primary_nvals; j++)
   {
      os << pcs_primary_function_name[j] << " ";
      idx[j] = GetNodeValueIndex(pcs_primary_function_name[j]);
      os << pcs_primary_function_name[j] << " ";
      idx[j + pcs_number_of_primary_nvals] = idx[j] + 1;
   }
   if (this->type == 12)                          //output of old & new time step for secondary variables PRESSURE2 and SATURATION1
   {
      os << pcs_secondary_function_name[0] << " ";
      idx[2* pcs_number_of_primary_nvals + 0] = GetNodeValueIndex(
         pcs_secondary_function_name[0]);
      os << pcs_secondary_function_name[0] << " ";
      idx[2* pcs_number_of_primary_nvals + 1] = GetNodeValueIndex(
         pcs_secondary_function_name[0]) + 1;
                                                  //other secondary variables
      for (j = 2; j < pcs_number_of_secondary_nvals; j++)
      {
         os << pcs_secondary_function_name[j] << " ";
         idx[2* pcs_number_of_primary_nvals + j] = GetNodeValueIndex(
            pcs_secondary_function_name[j]);
      }
   }
   else
   {
      for (j = 0; j < pcs_number_of_secondary_nvals; j++)
      {
         os << pcs_secondary_function_name[j] << " ";
         idx[2* pcs_number_of_primary_nvals + j] = GetNodeValueIndex(
            pcs_secondary_function_name[j]);
      }
   }
   os << endl;
   long i;
   for (i = 0; i < m_msh->GetNodesNumber(false); i++)
   {
      for (j = 0; j < 2* pcs_number_of_primary_nvals ; j++)
         os << GetNodeValue(i, idx[j]) << "  ";
      for (j = 0; j < pcs_number_of_secondary_nvals; j++)
         os << GetNodeValue(i, idx[2* pcs_number_of_primary_nvals + j])
            << "  ";
      os << endl;
   }
   os.close();
}


/**************************************************************************
PCSLib-Method:
based on MMPCalcSecondaryVariables
01/2007 OK Implementation
08/2008 CB NAPLdissolution
**************************************************************************/
void MMPCalcSecondaryVariablesNew(CRFProcess*m_pcs, bool NAPLdiss)
{
   long i;

   //----------------------------------------------------------------------
   int ndx_density_phase;
   int ndx_viscosity_phase;
   ndx_density_phase = -1;                        //WW
   ndx_viscosity_phase = -1;                      //WW
   CFEMesh* m_msh = m_pcs->m_msh;                 // PCH
   //----------------------------------------------------------------------
   m_pcs->SetBC();

   // For accessing the other process
   CRFProcess *cpl_pcs = NULL;
   if(m_pcs->pcs_type_number==0)
      cpl_pcs = pcs_vector[m_pcs->pcs_number+1];
   else if(m_pcs->pcs_type_number==1)
      cpl_pcs = pcs_vector[m_pcs->pcs_number-1];

   int ndx_pressure1, ndx_p_cap, ndx_pressure2, ndx_s_wetting, ndx_s_nonwetting;
   //======================================================================
   //----------------------------------------------------------------------
   // Capillary pressure - p_c (S) <- This is always the secondary variable
   // in both phase1 and phase2	// PCH
   CMediumProperties *mmp = NULL;

   //======================================================================
   switch(m_pcs->pcs_type_number)
   {
      case 0:
         //..................................................................
         //..................................................................
         // PCH
         // The primary variable is PRESSURE1
         // From PRESSURE1, we are assigning PRESSURE2 which is
         // the secondary variables of PRESSURE1.
         ndx_pressure1 = m_pcs->GetNodeValueIndex("PRESSURE1");
         ndx_pressure2 = m_pcs->GetNodeValueIndex("PRESSURE2");
         ndx_p_cap = m_pcs->GetNodeValueIndex("PRESSURE_CAP");
         double pressure1, pressure2, p_cap;
         for(i=0;i<(long)m_pcs->m_msh->nod_vector.size();i++)
         {
                                                  // New
            pressure1 = m_pcs->GetNodeValue(i,ndx_pressure1+1);
                                                  // Old
            pressure2 = m_pcs->GetNodeValue(i,ndx_pressure1);

            // Let's get capillary pressure before updating pressure2
            // by accessing the primary variable of the saturation equation
            // not the secondary variable of it.
            int cpl_ndx_sat2 = cpl_pcs->GetNodeValueIndex("SATURATION2");
            double cpl_sat2 = cpl_pcs->GetNodeValue(i,cpl_ndx_sat2+1);

            if(mmp_vector.size() > 1)
            {
               double sum = 0.0;
               CNode* thisNode = m_msh->nod_vector[i];
               int NumOfNeighborElements = (int)thisNode->connected_elements.size();
               // Harmonic mean
               for(int p=0; p< NumOfNeighborElements; ++p)
               {
                  // Mount neighboring elemenets and get the corresponding material group one by one.
                  int eleIdx = thisNode->connected_elements[p];
                  CElem* thisEle = m_msh->ele_vector[eleIdx];
                  int matgrp = thisEle->GetPatchIndex();
                  mmp = mmp_vector[matgrp];
                  mmp->mode = 2;
                  sum += 1.0/MRange(mmp->saturation_res[1],cpl_sat2,1.0-mmp->saturation_res[0]);
               }
               cpl_sat2 = (double)NumOfNeighborElements/sum;
            }
            else
            {
               mmp = mmp_vector[0];
               cpl_sat2 = MRange(mmp->saturation_res[1],cpl_sat2,1.0-mmp->saturation_res[0]);
            }
            // Assigning the secondary variable, Pc
            if(mmp_vector.size() > 1)
               p_cap = m_pcs->GetCapillaryPressureOnNodeByNeighobringElementPatches(i, 2, 1.0-cpl_sat2);
            else
               p_cap = mmp->CapillaryPressureFunction(i,NULL,1.0,0,1.0-cpl_sat2);

            m_pcs->SetNodeValue(i,ndx_p_cap,p_cap);
            m_pcs->SetNodeValue(i,ndx_p_cap+1,p_cap);

            pressure2 = pressure1 + p_cap;
            // Assigning the secondary variables
                                                  // Previous
            m_pcs->SetNodeValue(i,ndx_pressure2,pressure2);
                                                  // Now
            m_pcs->SetNodeValue(i,ndx_pressure2+1,pressure2);
         }
         //......................................................................
         ndx_density_phase = m_pcs->GetNodeValueIndex("DENSITY1");
         ndx_viscosity_phase = m_pcs->GetNodeValueIndex("VISCOSITY1");
         printf("Pressure2 from the known Pressure1 is updated for Process 0\n");
         break;

      case 1:
         //..................................................................
         // PCH
         // Calc secondary variable saturation Snonwetting = 1-Swetting
         // Don't forget here the primary variable is SATURATION2
         // From SATURATION2, we are assigning SATURATION1 which is
         // the secondary variables of SATURATION2.
         ndx_s_wetting = m_pcs->GetNodeValueIndex("SATURATION1");
         ndx_s_nonwetting = m_pcs->GetNodeValueIndex("SATURATION2");
         ndx_p_cap = cpl_pcs->GetNodeValueIndex("PRESSURE_CAP");

         double s_wetting,s_nonwetting;
         for(i=0;i<(long)m_pcs->m_msh->nod_vector.size();i++)
         {
            s_nonwetting = m_pcs->GetNodeValue(i,ndx_s_nonwetting+1);
            // Due to the iterative solution scheme in solving Snw with no
            // explicit boundary condition for non-zero flux condition,
            // Snw may become negative particularly the density difference
            // between two fluids is big. To prevent negative Snw, the
            // saturation restriction added.
            if(mmp_vector.size() > 1)
            {
               double sum = 0.0;
               CNode* thisNode = m_msh->nod_vector[i];
               int NumOfNeighborElements = (int)thisNode->connected_elements.size();
               // Harmonic mean
               for(int p=0; p< NumOfNeighborElements; ++p)
               {
                  // Mount neighboring elemenets and get the corresponding material group one by one.
                  int eleIdx = thisNode->connected_elements[p];
                  CElem* thisEle = m_msh->ele_vector[eleIdx];
                  int matgrp = thisEle->GetPatchIndex();
                  mmp = mmp_vector[matgrp];
                  mmp->mode = 2;
                  sum += 1.0/MRange(mmp->saturation_res[1],s_nonwetting,1.0-mmp->saturation_res[0]);
               }
               s_nonwetting = (double)NumOfNeighborElements/sum;
            }
            else
            {
               mmp = mmp_vector[0];
               s_nonwetting = MRange(mmp->saturation_res[1],s_nonwetting,1.0-mmp->saturation_res[0]);
            }
            // Assigning the secondary variable, Pc
            if(mmp_vector.size() > 1)
               p_cap = m_pcs->GetCapillaryPressureOnNodeByNeighobringElementPatches(i, 2, 1.0-s_nonwetting);
            else
               p_cap = mmp->CapillaryPressureFunction(i,NULL,1.0,0,1.0-s_nonwetting);

            m_pcs->SetNodeValue(i,ndx_s_nonwetting,s_nonwetting);
            m_pcs->SetNodeValue(i,ndx_s_nonwetting+1,s_nonwetting);
            s_wetting = 1.0-s_nonwetting;

            // Assigning the secondary variables
                                                  // Previous
            m_pcs->SetNodeValue(i,ndx_s_wetting,s_wetting);
                                                  // Now
            m_pcs->SetNodeValue(i,ndx_s_wetting+1,s_wetting);

            cpl_pcs->SetNodeValue(i,ndx_p_cap,p_cap);
            cpl_pcs->SetNodeValue(i,ndx_p_cap+1,p_cap);
         }
         //......................................................................
         ndx_density_phase = m_pcs->GetNodeValueIndex("DENSITY2");
         ndx_viscosity_phase = m_pcs->GetNodeValueIndex("VISCOSITY2");
         printf("Saturation1 from the known Saturation2 is updated for Process 1\n");
         break;
   }

   //----------------------------------------------------------------------
   // Fluid properties
   double density;
   double viscosity;
   CFluidProperties* m_mfp = NULL;
   m_mfp = mfp_vector[m_pcs->pcs_type_number];
   m_mfp->mode = 1;
   for(i=0;i<(long)m_pcs->m_msh->nod_vector.size();i++)
   {
      // CB NAPL dissolution reqiuires update of Density based on new composition of NAPL phase
                                                  //CB phase 2
      if(NAPLdiss==true && m_pcs->pcs_type_number==1)
         density = CalcNAPLDens(m_pcs, i);
      else
         density = m_mfp->Density();
      m_pcs->SetNodeValue(i,ndx_density_phase,density);
      viscosity = m_mfp->Viscosity();
      m_pcs->SetNodeValue(i,ndx_viscosity_phase,viscosity);
   }
   m_mfp->mode = 0;
   //----------------------------------------------------------------------
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
bool CRFProcess::OBJRelations()
{
   bool succeed = true;
   std::cout << "OBJ->PCS relations" << '\n';

   std::string pcs_type_name (convertProcessTypeToString (this->getProcessType()));

   // NUM
   std::cout << " - NUM->PCS" << '\n';
   m_num = NUMGet(pcs_type_name);
   if (!m_num)
   {
      std::cout << "Warning in CRFProcess::Create() - no NUM data - default"
         << std::endl;
      succeed = false;
   }
   else
   {
      pcs_nonlinear_iterations = m_num->nls_max_iterations;
      pcs_nonlinear_iteration_tolerance = m_num->nls_error_tolerance;
      num_type_name = "NEW";                      //OK
   }

   // TIM
   cout << " - TIM->PCS" << '\n';
   Tim = TIMGet(pcs_type_name);
   if (Tim)
   {
      // Time unit factor //WW OK: -> TIM
      if (Tim->time_unit.find("MINUTE") != string::npos)
         time_unit_factor = 60.0;
      else if (Tim->time_unit.find("HOUR") != string::npos)
         time_unit_factor = 3600.0;
      else if (Tim->time_unit.find("DAY") != string::npos)
         time_unit_factor = 86400.0;
      else if (Tim->time_unit.find("MONTH") != string::npos)
         time_unit_factor = 2592000.0;
      else if (Tim->time_unit.find("YEAR") != string::npos)
         time_unit_factor = 31536000;
   }
   else
   {
      cout << "Warning in CRFProcess::Create() - no TIM data - default"
         << endl;
   }

   // OUT

   // MSH
   if (fem_msh_vector.size() == 1)
      m_msh = fem_msh_vector[0];
   else
      m_msh = MSHGet(pcs_type_name);
   if (!m_msh)
   {
      cout << "Warning in CRFProcess::Create() - no MSH data" << endl;
      succeed = false;
   }

   return succeed;
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
bool CRFProcess::NODRelations()
{
   bool succeed = true;
   size_t DOF = GetPrimaryVNumber();              //OK should be PCS member variable

   cout << "NOD->PCS relations" << '\n';

   // BC
   cout << " - BC->PCS" << '\n';
   if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC")
      != string::npos)                            //WW
   {
      setBC_danymic_problems();
   }
   else
   {
      // create BC groups for each process
      CBoundaryConditionsGroup *m_bc_group = NULL;
      std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));
      for (size_t i = 0; i < DOF; i++)
      {
         BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
         m_bc_group = new CBoundaryConditionsGroup();
                                                  //OK
         m_bc_group->setProcessTypeName(pcs_type_name);
         m_bc_group->setProcessPrimaryVariableName(
            pcs_primary_function_name[i]);        //OK
         m_bc_group->Set(this, Shift[i]);
      }
   }

   // ST
   if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC")
      != string::npos)                            //WW
   {
      setST_danymic_problems();
   }

   CSourceTermGroup *m_st_group = NULL;
   if (WriteSourceNBC_RHS == 2)                   // Read from file
      ReadRHS_of_ST_NeumannBC();
   else                                           // WW
   {                                              // Calculate directly
      std::string pcs_type_name (convertProcessTypeToString(this->getProcessType()));
      for (size_t i = 0; i < DOF; i++)
      {
         m_st_group = STGetGroup(pcs_type_name, pcs_primary_function_name[i]);
         if (!m_st_group)
         {
            m_st_group = new CSourceTermGroup();
                                                  //OK
            m_st_group->pcs_type_name = pcs_type_name;
                                                  //OK
            m_st_group->pcs_pv_name = pcs_primary_function_name[i];
            m_st_group->Set(this, Shift[i]);
         }
      }
      if (WriteSourceNBC_RHS == 1)                // WW
         WriteRHS_of_ST_NeumannBC();
   }

   // NOD values
   cout << "->Config NOD values" << '\n';

   // Names
   nod_val_name_vector.clear();
   //
   for (int i = 0; i < pcs_number_of_primary_nvals; i++)
   {
                                                  // new time
      nod_val_name_vector.push_back(pcs_primary_function_name[i]);
                                                  // old time //need this MB!
      nod_val_name_vector.push_back(pcs_primary_function_name[i]);
   }
   for (int i = 0; i < pcs_number_of_secondary_nvals; i++)
                                                  // new time
      nod_val_name_vector.push_back(pcs_secondary_function_name[i]);
   if ((int) nod_val_name_vector.size() != (2* pcs_number_of_primary_nvals
      + pcs_number_of_secondary_nvals))
      succeed = false;

   // Values
   double* nod_values = NULL;
   const size_t nod_val_vector_size (nod_val_vector.size());
   for (size_t i = 0; i < nod_val_vector_size; i++)
   {
      delete nod_val_vector[i];
      nod_val_vector[i] = NULL;
   }
   nod_val_vector.clear();
   //
   //OK m_msh->NodesNumber_Quadratic;
   number_of_nvals = 2* DOF + pcs_number_of_secondary_nvals;
   for (long j = 0; j < (long) m_msh->nod_vector.size(); j++)
   {
      nod_values = new double[number_of_nvals];
      for (int i = 0; i < number_of_nvals; i++)
         nod_values[i] = 0.0;
      nod_val_vector.push_back(nod_values);
   }
   if ((long) nod_val_vector.size() != (long) m_msh->nod_vector.size())
      succeed = false;

   // IC
   cout << "->Assign IC" << '\n';
   if (reload == 2 && type != 4 && type != 41)
      ReadSolution();                             //WW
   SetIC();

   if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC")
      != string::npos)                            //WW
      setIC_danymic_problems();
   return succeed;
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
bool CRFProcess::ELERelations()
{
   bool succeed = true;
   //OK->MB please shift to Config()
   //	if (_pcs_type_name.compare("GROUNDWATER_FLOW") == 0)
   if (this->getProcessType() == GROUNDWATER_FLOW)
      MSHDefineMobile(this);
   //
   if (type == 4 || type == 41)
      m_msh->SwitchOnQuadraticNodes(true);
   else
      m_msh->SwitchOnQuadraticNodes(false);
   CheckMarkedElement();

   // ELE - GP values
   AllocateMemGPoint();

   // ELE values
   double* ele_values = NULL;                     // PCH
   int number_of_evals = 2* pcs_number_of_evals ; //PCH, increase memory
   if (number_of_evals > 0)                       // WW added this "if" condition
   {
      // Names
      for (int i = 0; i < pcs_number_of_evals; i++)
      {
                                                  // new time
         ele_val_name_vector.push_back(pcs_eval_name[i]);
                                                  // old time
         ele_val_name_vector.push_back(pcs_eval_name[i]);
      }
      if (ele_val_name_vector.size() == 0)
         succeed = false;

      // Values
      size_t m_msh_ele_vector_size (m_msh->ele_vector.size());
      if (ele_val_vector.size() == 0)
      {
         for (size_t j = 0; j < m_msh_ele_vector_size; j++)
         {
            ele_values = new double[number_of_evals];
            size_eval += number_of_evals;         //WW
            for (int i = 0; i < number_of_evals; i++)
               ele_values[i] = 0.0;
            ele_val_vector.push_back(ele_values);
         }
      }
      else
      {
         for (size_t j = 0; j < m_msh_ele_vector_size; j++)
         {
            ele_values = ele_val_vector[j];
            /* //Comment by WW
             #ifndef SX
             #ifdef GCC
             size = malloc_usable_size( ele_values )/sizeof(double);
             #elif HORIZON
             //KG44: malloc_usable_size and _msize are not available
             #else
             size= _msize( ele_values )/sizeof(double);
             #endif
             #endif
             */
            ele_values = resize(ele_values, size_eval, size_eval
               + number_of_evals);
            size_eval += number_of_evals;
            ele_val_vector[j] = ele_values;
         }
      }
      if (ele_val_vector.size() != m_msh->ele_vector.size())
         succeed = false;
   }

   // ELE matrices
   if (Memory_Type != 0)
   {
      AllocateLocalMatrixMemory();
      if ((long) Ele_Matrices.size() != (long) m_msh->ele_vector.size())
         succeed = false;
   }

   // Element matrix output. WW
   if (Write_Matrix)
   {
      cout << "->Write Matrix" << '\n';
      string m_file_name = FileName + "_" + convertProcessTypeToString (this->getProcessType())
         + "_element_matrix.txt";
      matrix_file = new fstream(m_file_name.c_str(), ios::trunc | ios::out);
      if (!matrix_file->good())
         cout << "Warning in GlobalAssembly: Matrix files are not found"
            << endl;
   }

   // FEM
   if (type == 4 || type == 41)
   {
      // Set initialization function
      CRFProcessDeformation *dm_pcs = (CRFProcessDeformation *) this;
      dm_pcs->Initialization();
      if (!dm_pcs->fem_dm)
         succeed = false;
   } else                                         // Initialize FEM calculator
   {
      int Axisymm = 1;                            // ani-axisymmetry
      if (m_msh->isAxisymmetry())
         Axisymm = -1;                            // Axisymmetry is true
                                                  //OK4801 needs NUM
      fem = new CFiniteElementStd(this, Axisymm * m_msh->GetCoordinateFlag());
      fem->SetGaussPointNumber(m_num->ele_gauss_points);
      if (!fem)
         succeed = false;
   }

   return succeed;
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
#ifndef NEW_EQS                                   //WW 07.11.2008
bool CRFProcess::CreateEQS()
{
   if(!m_num) return false;                       //OK46
   bool succeed = true;
   //----------------------------------------------------------------------------
   if(eqs)
      return false;
   //----------------------------------------------------------------------------
   int DOF = GetPrimaryVNumber();                 //OK should be PCS member variable
   //----------------------------------------------------------------------------
   // EQS - create equation system
   cout << "->Create EQS" << '\n';
   //----------------------------------------------------------------------------
   if(type==4)
   {
      eqs = CreateLinearSolverDim(m_num->ls_storage_method,DOF,DOF*m_msh->GetNodesNumber(true));
      InitializeLinearSolver(eqs,m_num);
      PCS_Solver.push_back(eqs);                  //WW
   }
   //----------------------------------------------------------------------------
   else if(type==41)
   {
      if(num_type_name.find("EXCAVATION")!=string::npos)
         eqs = CreateLinearSolverDim(m_num->ls_storage_method,DOF-1,DOF*m_msh->GetNodesNumber(true));
      else
         eqs = CreateLinearSolverDim(m_num->ls_storage_method,DOF,
            (DOF-1)*m_msh->GetNodesNumber(true)+m_msh->GetNodesNumber(false));
      InitializeLinearSolver(eqs,m_num);
      PCS_Solver.push_back(eqs);                  //WW
   }
   //----------------------------------------------------------------------------
   else
   {
      /*
          // If there is a solver exsiting. WW
          CRFProcess* m_pcs = NULL;
          for(int i=0; i<(int)pcs_vector.size(); i++)
         {
            m_pcs = pcs_vector[i];
            if(m_pcs&&m_pcs->eqs)
           {
              if(m_pcs->_pcs_type_name.find("DEFORMATION")==string::npos)
                break;
           }
      }
      // If unique mesh
      if(m_pcs&&m_pcs->eqs&&(fem_msh_vector.size()==1))
      eqs = m_pcs->eqs;
      //
      else
      {
      */
      eqs = CreateLinearSolver(m_num->ls_storage_method,m_msh->GetNodesNumber(false));
      InitializeLinearSolver(eqs,m_num);
      PCS_Solver.push_back(eqs);
      //}
   }
   //----------------------------------------------------------------------------
   std::string pcs_type_name (convertProcessTypeToString (this->getProcessType()));
   strcpy(eqs->pcs_type_name,pcs_type_name.data());
   //----------------------------------------------------------------------------
   if((int)PCS_Solver.size()==0)
      succeed = false;
   //----------------------------------------------------------------------------
   return succeed;
}
#endif
/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
#ifndef NEW_EQS                                   //WW. 07.11.2008
void PCSCreateNew()
{
   int i;
   CRFProcess* m_pcs = NULL;
   //----------------------------------------------------------------------
   for(i=0;i<(int)pcs_vector.size();i++)
   {
      m_pcs = pcs_vector[i];
      m_pcs->CreateNew();
      //----------------------------------------------------------------------
   }
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
void CRFProcess::CreateNew()
{
   pcs_type_number = (int)pcs_vector.size();
   Config();
   m_bCheckOBJ = OBJRelations();
   m_bCheckEQS = CreateEQS();
   m_bCheckNOD = NODRelations();
   m_bCheckELE = ELERelations();
   MMP2PCSRelation(this);
   ConfigureCouplingForLocalAssemblier();
}
#endif
/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
bool CRFProcess::Check()
{
   // MMP
   MSHTestMATGroups();
   return true;
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
bool PCSCheck()
{
   if((int)pcs_vector.size()==0)
      return false;
   CRFProcess* m_pcs = NULL;
   for(int i=0;i<(int)pcs_vector.size();i++)
   {
      m_pcs = pcs_vector[i];
      //if(m_pcs->m_bCheck)
      if(!m_pcs->Check())
         return false;
#ifdef MFC
      FMRead();                                   //WW
#endif
   }

   return true;
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
#ifndef NEW_EQS                                   //WW 07.11.2008
void EQSDelete()
{
   LINEAR_SOLVER *eqs = NULL;
   CRFProcess* m_pcs = NULL;
   //----------------------------------------------------------------------
   for(size_t i=0;i<PCS_Solver.size();i++)
   {
      eqs = PCS_Solver[i];
      ProcessType pcs_type (convertProcessType (eqs->pcs_type_name));
      m_pcs = PCSGet(pcs_type);
      if(eqs->unknown_vector_indeces)
         eqs->unknown_vector_indeces = \
            (int*) Free(eqs->unknown_vector_indeces);
      if(eqs->unknown_node_numbers)
         eqs->unknown_node_numbers = \
            (long*) Free(eqs->unknown_node_numbers);
      if(eqs->unknown_update_methods)
         eqs->unknown_update_methods = \
            (int*) Free(eqs->unknown_update_methods);
      eqs = DestroyLinearSolver(eqs);
      if(m_pcs)
         m_pcs->eqs = NULL;
      PCS_Solver.erase((PCS_Solver.begin()+i));
   }
   //PCS_Solver.clear();
}
#endif
/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
void CRFProcess::NODRelationsDelete()
{
   std::string pcs_type_name (convertProcessTypeToString (getProcessType()));
   // BC
   for (size_t i = 0; i < GetPrimaryVNumber(); i++)
   {
      BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
   }
   const size_t bc_node_value_size(bc_node_value.size());
   for (size_t i = 0; i < bc_node_value_size; i++)
   {
      delete bc_node_value[i];
      bc_node_value[i] = NULL;
   }
   bc_node_value.clear();

   // ST
   for (size_t i = 0; i < GetPrimaryVNumber(); i++)
   {
      STGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
   }

   CNodeValue* nod_val = NULL;
   const size_t st_node_value_size (st_node_value.size());
   for (size_t i = 0; i < st_node_value_size; i++)
   {
      nod_val = st_node_value[i];
      //OK delete st_node_value[i];
      //OK st_node_value[i] = NULL;
      if (nod_val->check_me)                      //OK
      {
         nod_val->check_me = false;
         delete nod_val;
         nod_val = NULL;
      }
   }
   st_node_value.clear();

   // NOD values
   nod_val_name_vector.clear();
   const size_t nod_val_vector_size(nod_val_vector.size());
   for (size_t i = 0; i < nod_val_vector_size; i++)
   {
      delete nod_val_vector[i];
      nod_val_vector[i] = NULL;
   }
   nod_val_vector.clear();
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
void CRFProcess::ELERelationsDelete()
{
   long i;
   //----------------------------------------------------------------------
   // FEM element
   if(fem) delete fem;                            //WW
   fem = NULL;
   //----------------------------------------------------------------------
   // ELE matrices
   ElementMatrix *eleMatrix = NULL;
   ElementValue* gp_ele = NULL;
   if(Ele_Matrices.size()>0)
   {
      for (i=0;i<(long)Ele_Matrices.size();i++)
      {
         eleMatrix = Ele_Matrices[i];
         delete eleMatrix;
         eleMatrix = NULL;
      }
      Ele_Matrices.clear();
   }
   //----------------------------------------------------------------------
   // ELE - GP values
   if(ele_gp_value.size()>0)
   {
      for(i=0;i<(long)ele_gp_value.size();i++)
      {
         gp_ele = ele_gp_value[i];
         delete gp_ele;
         gp_ele = NULL;
      }
      ele_gp_value.clear();
   }
   //----------------------------------------------------------------------
   // ELE values
   ele_val_name_vector.clear();
   for(i=0;i<(long)ele_val_vector.size();i++)
   {
      delete ele_val_vector[i];
      //delete[] ele_val_vector[i];
      ele_val_vector[i] = NULL;
   }
   ele_val_vector.clear();
   //----------------------------------------------------------------------
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
void CRFProcess::OBJRelationsDelete()
{
   //----------------------------------------------------------------------------
   cout << "OBJ->PCS relations delete" << '\n';
   //----------------------------------------------------------------------------
   m_num = NULL;
   Tim =  NULL;
   m_msh = NULL;
   //----------------------------------------------------------------------------
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
void CRFProcess::Delete()
{
   //----------------------------------------------------------------------------
   cout << "PCS  delete" << '\n';
   //----------------------------------------------------------------------------
   ELERelationsDelete();
   NODRelationsDelete();
#ifndef NEW_EQS                                //WW. 07.11.2008
   EQSDelete();
#endif
   OBJRelationsDelete();
   //MMP2PCSRelation(this);
   //ConfigureCouplingForLocalAssemblier();
   //----------------------------------------------------------------------------
}


/**************************************************************************
PCSLib-Method:
07/2007 OK Implementation
**************************************************************************/
#ifndef NEW_EQS                                   //WW 07.11.2008
void CRFProcess::EQSDelete()
{
   std::string pcs_type_name (convertProcessTypeToString (this->getProcessType()));
   LINEAR_SOLVER *eqs = NULL;
   for (size_t i = 0; i < PCS_Solver.size(); i++)
   {
      eqs = PCS_Solver[i];
      if (pcs_type_name.compare(eqs->pcs_type_name) == 0)
      {
         if (eqs->unknown_vector_indeces)
            eqs->unknown_vector_indeces = (int*) Free(
               eqs->unknown_vector_indeces);
         if (eqs->unknown_node_numbers)
            eqs->unknown_node_numbers = (long*) Free(
               eqs->unknown_node_numbers);
         if (eqs->unknown_update_methods)
            eqs->unknown_update_methods = (int*) Free(
               eqs->unknown_update_methods);
         eqs = DestroyLinearSolver(eqs);
         eqs = NULL;
      }
   }

   for (size_t i = 0; i < PCS_Solver.size(); i++)
   {
      eqs = PCS_Solver[i];
      if (pcs_type_name.compare(eqs->pcs_type_name) == 0)
         PCS_Solver.erase((PCS_Solver.begin() + i));
   }
}
#endif

/*************************************************************************
ROCKFLOW - Function: CRFProcess::
Task:
Programming:
09/2007 WW Implementation
**************************************************************************/

// New solvers WW
#ifdef NEW_EQS                                    //1.09.2007 WW

void CreateEQS_LinearSolver()
{
   int dof_DM=1;
   CRFProcess *m_pcs = NULL;
   CFEMesh *a_msh = NULL;
   SparseTable *sp = NULL;
   SparseTable *spH = NULL;
   Linear_EQS *eqs = NULL;
   Linear_EQS *eqsH = NULL;
   int dof = 1;
   //
   size_t dof_nonDM (1);

   for(size_t i=0;i<pcs_vector.size();i++)
   {
      m_pcs = pcs_vector[i];
      if(m_pcs->type==22)                         // Monolithic TH2
      {
         dof_nonDM = m_pcs->GetPrimaryVNumber();
         dof = dof_nonDM;
      }
      if(m_pcs->type==4||m_pcs->type==41)         // Deformation
      {
         dof_DM = m_pcs->GetPrimaryVNumber();
         dof = dof_DM;
      }
      else                                        // Monolithic scheme for the process with linear elements
      {
         if(dof_nonDM < m_pcs->GetPrimaryVNumber())
         {
            dof_nonDM = m_pcs->GetPrimaryVNumber();

            // PCH: DOF Handling for FLUID_MOMENTUM in case that the LIS and PARDISO solvers
            // are chosen.
            //				if(m_pcs->_pcs_type_name.compare("FLUID_MOMENTUM")==0)
            if(m_pcs->getProcessType() == FLUID_MOMENTUM)
               dof_nonDM = 1;
         }
      }
   }
   //
   for(size_t i=0; i<fem_msh_vector.size(); i++)
   {
      a_msh = fem_msh_vector[i];
#if defined(USE_MPI)
      eqs = new Linear_EQS(a_msh->GetNodesNumber(true)*dof);
      EQS_Vector.push_back(eqs);
      EQS_Vector.push_back(eqs);
#else
      a_msh->CreateSparseTable();
      // sparse pattern with linear elements exists
      sp = a_msh->GetSparseTable();
      // sparse pattern with quadratic elements exists
      spH = a_msh->GetSparseTable(true);
      //
      eqs=NULL;
      eqsH=NULL;
      if(sp)
         eqs = new Linear_EQS(*sp, dof_nonDM);
      if(spH)
         eqsH = new Linear_EQS(*spH, dof_DM);
      EQS_Vector.push_back(eqs);
      EQS_Vector.push_back(eqsH);
#endif
   }
}


#else                                             // NEW_EQS
/*************************************************************************
ROCKFLOW - Function: CRFProcess::
Task:
Programming:
09/2007 WW Implementation
**************************************************************************/
#include <iomanip>
void CRFProcess::DumpEqs(string file_name)
{
   int ii;
   long i, j, k, m;
   fstream eqs_out;
   CNode *node;
   eqs_out.open(file_name.c_str(), ios::out );
   eqs_out.setf(ios::scientific, ios::floatfield);
   setw(14);
   eqs_out.precision(14);
   //
   long nnode = eqs->dim/eqs->unknown_vector_dimension;
   for(i=0; i<eqs->dim; i++)
   {
      m = i%nnode;
      node = m_msh->nod_vector[m];
      for(ii=0; ii<eqs->unknown_vector_dimension; ii++)
      {
         for(j=0; j<(long)node->connected_nodes.size(); j++)
         {
            k = node->connected_nodes[j];
            k = ii*nnode+k;
            if(k>=eqs->dim)  continue;
            eqs_out<<i<<"  "<<k<<"  "<<MXGet(i,k)<<endl;
         }
      }
      eqs_out<<   i<<"  "<<   eqs->b[i]<<"  "<<eqs->x[i]<<endl;
   }
   eqs_out.close();
}
#endif                                            //ifdef NEW_QES

/*************************************************************************
ROCKFLOW - Function: CRFProcess::
Task:
Programming:
01/2008 WW Implementation
**************************************************************************/
void CRFProcess::WriteBC()
{
   const size_t size_bc (bc_node_value.size());
   const size_t size_st (st_node_value.size());

   if (size_bc == 0 && size_st == 0)
      return;

   std::string m_file_name = FileName + "_" + convertProcessTypeToString(this->getProcessType()) + "_BC_ST.asc";
   std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
   if (!os.good())
   {
      cout << "Failure to open file: " << m_file_name << endl;
      abort();
   }
   os.setf(ios::scientific, ios::floatfield);
   os.precision(12);
   CNode *anode = NULL;
   long nindex = 0;
   if (size_bc > 0)
   {
      os << "#Dirchilet BC  (from " << m_file_name << ".bc file) " << endl;
      os << "#Total BC nodes  " << size_bc << endl;
      os << "#Node index, name, x, y, z,   value: " << endl;
      for (size_t i = 0; i < size_bc; i++)
      {
         nindex = bc_node_value[i]->geo_node_number;
         anode = m_msh->nod_vector[nindex];
         os << nindex << "  " << bc_node_value[i]->pcs_pv_name << " "
            << std::setw(14) << anode->X() << " " << std::setw(14) << anode->Y()
            << " " << std::setw(14) << anode->Z() << " " << std::setw(14)
            << bc_node_value[i]->node_value << endl;
      }
   }
   if (size_st > 0)
   {
      os << "#Source term or Neumann BC  (from " << m_file_name
         << ".st file) " << endl;
      os << "#Total ST nodes  " << size_st << endl;
      os << "#Node index, x, y, z, name    value: " << endl;
      for (size_t i = 0; i < size_st; i++)
      {
         nindex = st_node_value[i]->geo_node_number;
         anode = m_msh->nod_vector[nindex];
         os << nindex << "  " << convertPrimaryVariableToString(
            st_node[i]->getProcessPrimaryVariable()) << " " << std::setw(14)
            << anode->X() << " " << std::setw(14) << anode->Y() << " "
            << std::setw(14) << anode->Z() << " " << std::setw(14)
            << st_node_value[i]->node_value << endl;
      }
   }
   os << "#STOP" << endl;
   os.close();
}


/*************************************************************************
GeoSys-Function:
Task: PI time step contorl
Programming:
08/2008 WW Implementation
10/2008 WW Node value criteria (test)
03/2009 WW Euclidean norm
**************************************************************************/
void CRFProcess::PI_TimeStepSize(double *u_n)
{
   //----------------------------------------------------------------------
   //----------------------------------------------------------------------
   // Time step control
   double hmin;
   double hmax;
   double factor1;                                // 1/hmin
   double factor2;                                // 1/hmax
   double sfactor = 0.9;
   //
   //
   hmax = Tim->GetMaximumTSizeRestric();
   hmin = Tim->GetMinimumTSizeRestric();
   //
   if (hmin<DBL_MIN) factor1 = 5.0;
   else factor1 = 1.0/hmin;
   if (hmax<DBL_MIN) factor2 = 0.166666666666666666667e+00;
   else  factor2 = 1.0/hmax;
   if (factor1<1.0e0) factor1 = 5.0;
   if (factor2>1.0e0) factor2 = 0.166666666666666666667e+00;
   //
   hmax = Tim->max_time_step;
   if(hmax<DBL_MIN)
      hmax = fabs(Tim->time_end-aktuelle_zeit);
   //
   //esitmate the error
   double hnew;
   double err, fac;
   double factorGus;
   double hacc = Tim->GetHacc();
   double erracc = Tim->GetErracc();
   //
#define aE_NORM
#ifdef E_NORM
   //
   long i;
   CElem* elem = NULL;
   bool Check2D3D;
   double norm_e, norm_en;
   double norm_e_rank, norm_en_rank;
   norm_e = norm_en = norm_e_rank = norm_en_rank = 0.;

   Check2D3D = false;
   if(type == 66)                                 //Overland flow
      Check2D3D = true;
   //----------------------------------------------------------------------
   // DDC
   if(dom_vector.size()>0)
   {
      cout << "      Domain Decomposition" << '\n';
      CPARDomain* m_dom = NULL;
      int j = 0;
      //
#if defined(USE_MPI)
      j = myrank;
#else
      for(j=0;j<(int)dom_vector.size();j++)
      {
#endif
         m_dom = dom_vector[j];
         for(int ii=0;ii<(int)continuum_vector.size();ii++)
         {
            continuum = ii;
            //
            for(i=0;i<(long)m_dom->elements.size();i++)
            {
               elem = m_msh->ele_vector[m_dom->elements[i]];
               if(elem->GetMark())
               {
                  elem->SetOrder(false);
                  fem->SetElementNodesDomain(m_dom->element_nodes_dom[i]);
                  fem->ConfigElement(elem,Check2D3D);
                  fem->m_dom = m_dom;
                  fem->CalcEnergyNorm(u_n, norm_e_rank, norm_en_rank);
                  // _new
                  if(ii==1)
                     fem->CalcEnergyNorm_Dual(u_n, norm_e_rank, norm_en_rank);
               }
            }
         }
#if defined(USE_MPI)
         MPI_Allreduce(&norm_e_rank, &norm_e, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
         MPI_Allreduce(&norm_en_rank, &norm_en, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else                                    //USE_MPI
         norm_e += norm_e_rank;
         norm_en += norm_en_rank;
      }
      //....................................................................
#endif
   }
   //----------------------------------------------------------------------
   // STD
   else
   {
      for(int ii=0;ii<(int)continuum_vector.size();ii++)
      {
         continuum = ii;
         for(i=0;i<(long)m_msh->ele_vector.size();i++)
         {
            elem = m_msh->ele_vector[i];
            if (elem->GetMark())                  // Marked for use
            {
               elem->SetOrder(false);
               fem->ConfigElement(elem,Check2D3D);
               fem->CalcEnergyNorm(u_n, norm_e, norm_en);
               // _new
               if(ii==1)
                  fem->CalcEnergyNorm_Dual(u_n, norm_e, norm_en);
            }
         }
      }
   }
   //compute energy norm as the error
   err = sqrt(fabs(norm_e/norm_en));
#else
   err = 0.0;
   //
   int ii, nidx1;
   long g_nnodes, j, k, l, size_x;
   double x0, x1;
   double Rtol = Tim->GetRTol();
   double Atol = Tim->GetATol();
   double *u_n0 = NULL;
#if defined(NEW_EQS)
   u_n0 = eqs_new->x;
#else
   u_n0 = eqs->x;
#endif

   size_x = 0;
   for(ii=0; ii<pcs_number_of_primary_nvals; ii++)
   {
      nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii])+1;
      g_nnodes =m_msh->GetNodesNumber(false);     //DOF>1, WW
      size_x += g_nnodes;
      for(j=0;j<g_nnodes;j++)
      {
         k = m_msh->Eqs2Global_NodeIndex[j];
         l = j+ii*g_nnodes;
         x0 = u_n[l];
         x1 = GetNodeValue(k,nidx1);
         err += pow( (x1- u_n0[l])/(Atol+Rtol*max(fabs(x0),fabs(x1))),2);
      }
   }
   err = sqrt(err/(double)size_x);
#endif
   //----------------------------------------------------------------------
   //
   //compute hnew with the resriction: 0.2<=hnew/h<=6.0;
   fac = max(factor2,min(factor1,pow(err,0.25)/sfactor));
   hnew = dt/fac;

   //determine if the error is small enough
   if(err<=1.0e0)                                 //step is accept
   {
      accept_steps++;
      if(Tim->GetTimeStepCrtlType()==2)           //Mod. predictive step size controller (Gustafsson)
      {
         if(accept_steps>1)
         {
            factorGus = (hacc/dt)*pow(err*err/erracc,0.25)/sfactor;
            factorGus = max(factor2, min(factor1,factorGus));
            fac = max(fac, factorGus);
            hnew=dt/fac;
         }
         hacc = dt;
         erracc = max(1.0e-2,err);
         Tim->SetHacc(hacc);
         Tim->setErracc(erracc);
      }
      if(fabs(hnew)>hmax) hnew = hmax;
      if(!accepted)
         hnew = min(fabs(hnew), dt);
      Tim->SetTimeStep(fabs(hnew));
      accepted = true;
   }                                              //end if(err<=1.0e0)
   else
   {
      reject_steps++;
      accepted = false;
      Tim->SetTimeStep(hnew);
      if(reject_steps>100&&accept_steps==0)
      {
         cout<<"!!! More than 100 steps rejected and none of steps accepted. Quit the simulation now"<<endl;
         exit(1);
      }
      // Recover solutions
   }

}


#ifdef NEW_EQS                                    //WW
/*************************************************************************
ROCKFLOW - Function: CRFProcess::
Task:  //For fluid momentum, WW
Programming:
01/2008 WW Implementation
**************************************************************************/
void CRFProcess::EQSInitialize()
{eqs_new->Initialize();}

#if defined(LIS)                                  //11.03.2008 WW
/*************************************************************************
ROCKFLOW - Function: CRFProcess::
Task:  //For fluid momentum,
Programming:
02/2008 PCH Implementation
03/2009 PCH option to tell if this is FLUID_MOMENTUM
**************************************************************************/
void CRFProcess::EQSSolver(double* x)
{
   eqs_new->Solver(this->m_num);                  //NW

                                                  //OK411
   for(int i=0; i < (int)m_msh->nod_vector.size(); ++i)
      x[i] = eqs_new->X(i);
}
#endif
#endif

#ifdef GEM_REACT
void CRFProcess::IncorporateSourceTerms_GEMS(void)
{
   // Initialization
   long it;                                       // iterator
   long N_Nodes;                                  // Number of Nodes
   int nDC;                                       // Number of mass transport components.
   int i =0 ;                                     // index of the component

   // Get a vector pointing to the REACT_GEM class
   if (m_vec_GEM)
   {
      // Get needed informations.----------------------------
      // Number of Nodes
      N_Nodes = (long) m_msh->GetNodesNumber(false);
      // Number of DC
      nDC = m_vec_GEM->nDC;
      // Identify which PCS it is and its sequence in mcp.
      i = this->pcs_component_number;
      // ----------------------------------------------------

      // Loop over all the nodes-----------------------------
      for ( it=0 ; it <  N_Nodes/*Number of Nodes*/; it++ )
      {
         // Adding the rate of concentration change to the right hand side of the equation.
#ifdef NEW_EQS                           //15.12.2008. WW
         eqs_new->b[it] -= m_vec_GEM->m_xDC_Chem_delta[it*nDC+i] / dt ;
#else
         eqs->b[it] -= m_vec_GEM->m_xDC_Chem_delta[it*nDC+i] / dt ;
#endif
      }
      // ----------------------------------------------------
   }
}
#endif
/*************************************************************************
ROCKFLOW - Function: CRFProcess::
Task:
Programming:
10/2008 //WW/CB Implementation
07/2010 TF substituted GEOGetPLYByName
03/2010  WW Read binary file of precipitation
06/2010  WW Output top surface flux and head to a raster file
12/2016  MJ Add the case of hexahedral mesh
**************************************************************************/
void CRFProcess::UpdateTransientBC()
{
   //--------------- 24.03.2010. WW
   long i;
   CSourceTerm *precip;
   CSourceTerm *a_st;
   precip  = NULL;

   for(i=0;i<(long)st_vector.size(); i++)
   {
      a_st = st_vector[i];
      if(  a_st->getProcessDistributionType() == FiniteElement::PRECIPITATION
         &&a_st->getProcessType() == getProcessType() )
      {
         precip = a_st;
         if(!m_msh->top_surface_checked)          // 07.06--19.08.2010. WW
         {
            if(m_msh->GetCoordinateFlag()/10==3)
               m_msh->MarkInterface_mHM_Hydro_3D();
            m_msh->top_surface_checked = true;
         }
         break;
      }
   }
   if(precip)                                     // 08-07.06.2010.  WW
   {
      string ofile_name;
      ofile_name = precip->DirectAssign_Precipitation(aktuelle_zeit);

      if( m_msh->GetCoordinateFlag()/10==3 )      //19.08.2010. WW
      {
         /// Remove .bin from the file name
         i = (int)ofile_name.find_last_of(".");
         if(i>0)
            ofile_name.erase(ofile_name.begin()+i,ofile_name.end());
         i = (int)ofile_name.find_last_of(".");
         if(i>0)
            ofile_name.erase(ofile_name.begin()+i,ofile_name.end());
         string of_name = ofile_name + ".flx.asc";
         ofstream of_flux(of_name.c_str(), ios::trunc|ios::out);

         of_name = ofile_name + ".pri.asc";
         ofstream of_primary(of_name.c_str(), ios::trunc|ios::out);

         /// GIS_shape_head[0]:  ncols
         /// GIS_shape_head[1]:  nrows
         /// GIS_shape_head[2]:  xllcorner
         /// GIS_shape_head[3]:  yllcorner
         /// GIS_shape_head[4]:  cellsize
         /// GIS_shape_head[5]:  NONDATA_value
         double *g_para = precip->GIS_shape_head;
         long size = (long)(g_para[0]*g_para[1]);
         vector<double> cell_data_p(size);
         vector<double> cell_data_v(size);
         for(i=0; i<size; i++)
         {
            cell_data_p[i] = g_para[5];
            cell_data_v[i] = g_para[5];
         }

         int j, k, nnodes;
         int node_xmin, node_xmax, node_ymin, node_ymax;
         long m, n, mm, nn, l;
         CElem *elem;
         CNode *node;
         double *cent;
         double vel_av[3], x1[3], x2[3], x3[3], x4[3], sub_area[3], area, inarea, tol_a;

         double x_min, y_min, x_max, y_max;
         long  row_min, col_min, row_max, col_max;

         int *vel_idx;
         int idx = fem->idx0+1;
         vel_idx = fem->idx_vel;

         long n_idx, irow, icol, nrow, ncol;
         nrow = (long)g_para[1];
         ncol = (long)g_para[0];

         node_xmin = node_xmax = node_ymin = node_ymax = 0;

         tol_a = 1.e-8;

         of_flux.setf(ios::fixed, ios::floatfield);
         of_primary.setf(ios::fixed, ios::floatfield);
         of_flux.precision(1);
         of_primary.precision(1);

         of_flux<<"ncols"<<setw(19)<<ncol<<endl;
         of_flux<<"nrows"<<setw(19)<<nrow<<endl;
         of_flux<<"xllcorner"<<setw(15)<<g_para[2]<<endl;
         of_flux<<"yllcorner"<<setw(15)<<g_para[3]<<endl;
         of_flux<<"cellsize"<<setw(16)<<(long)g_para[4]<<endl;
         of_flux<<"NODATA_value"<<setw(11)<<(long)g_para[5]<<endl;

         of_primary<<"ncols"<<setw(19)<<ncol<<endl;
         of_primary<<"nrows"<<setw(19)<<nrow<<endl;
         of_primary<<"xllcorner"<<setw(15)<<g_para[2]<<endl;
         of_primary<<"yllcorner"<<setw(15)<<g_para[3]<<endl;
         of_primary<<"cellsize"<<setw(16)<<(long)g_para[4]<<endl;
         of_primary<<"NODATA_value"<<setw(11)<<(long)g_para[5]<<endl;

         for(i=0; i<(long)m_msh->face_vector.size(); i++)
         {
            elem = m_msh->face_vector[i];
            if(!elem->GetMark())
               continue;

            if(elem->GetElementType() !=MshElemType::QUAD && elem->GetElementType() !=MshElemType::TRIANGLE )        /// If not triangle
                continue;

            //// In element
            nnodes = elem->GetNodesNumber(false);
            cent = elem->gravity_center;

            /// Find the range of this element （找包含每个element的cell的最小范围（不一定只有一个cell））
            x_min = y_min = 1.e+20;
            x_max = y_max = -1.e+20;
            for(k=0; k<nnodes; k++)
            {
               node = elem->nodes[k];

               if(node->X()<x_min)
               {
                  x_min = node->X();
                  node_xmin = k;
               }
               if(node->X()>x_max)
               {
                  x_max = node->X();
                  node_xmax = k;
               }
               if(node->Y()<y_min)
               {
                  y_min = node->Y();
                  node_ymin = k;
               }
               if(node->Y()>y_max)
               {
                  y_max = node->Y();
                  node_ymax = k;
               }
            }

            /// Determine the cells that this element covers（找到asc中对应的cell的行数和列数）. 05.10. 2010
            col_min = (long)((x_min-g_para[2])/g_para[4]); //这是asc里面的行数最小值
            row_min = nrow -(long)((y_max-g_para[3])/g_para[4]);
            col_max = (long)((x_max-g_para[2])/g_para[4]);
            row_max = nrow -(long)((y_min-g_para[3])/g_para[4]);

            if (elem->GetElementType() ==MshElemType::TRIANGLE)
            {
                node = elem->nodes[0];
                x1[0] = node->X();
                x1[1] = node->Y();
                node = elem->nodes[1];
                x2[0] = node->X();
                x2[1] = node->Y();
                node = elem->nodes[2];
                x3[0] = node->X();
                x3[1] = node->Y();

		x3[2] = x2[2] = x1[2] = 0.;
	        cent[2] = 0.;
            }

            else if (elem->GetElementType() ==MshElemType::QUAD)
            {
                node = elem->nodes[0];
                x1[0] = node->X();
                x1[1] = node->Y();
                node = elem->nodes[1];
                x2[0] = node->X();
                x2[1] = node->Y();
                node = elem->nodes[2];
                x3[0] = node->X();
                x3[1] = node->Y();
		        node = elem->nodes[3];
                x4[0] = node->X();
                x4[1] = node->Y();

		    x4[2] = x3[2] = x2[2] = x1[2] = 0.;
	        cent[2] = 0.;
            }

            for(m=col_min; m<=col_max; m++)
            {
               mm = m; //m 代表asc文件中的列数
               if(m>ncol-1) mm = ncol-1;
               if(m<0) mm = 0;
               cent[0] = g_para[2] + g_para[4]*(mm + 0.5); //cent[0]代表cell中心的x

               for(n=row_min; n<=row_max; n++)
               {
                  nn = n; //n 代表asc文件中的行数
                  if(n>nrow-1) nn = nrow-1;
                  if(nn<0) nn = 0;
                  cent[1] = g_para[3] + g_para[4]*(nrow-nn + 0.5); //cent[1]代表cell中心的y

                  if(cent[0]<x_min)
                  {
                     node = elem->nodes[node_xmin];
                     cent[0] = node->X();
                     cent[1] = node->Y();
                  }
                  if(cent[0]>x_max)
                  {
                     node = elem->nodes[node_xmax];
                     cent[0] = node->X();
                     cent[1] = node->Y();
                  }
                  if(cent[1]<y_min)
                  {
                     node = elem->nodes[node_ymin];
                     cent[0] = node->X();
                     cent[1] = node->Y();
                  }
                  if(cent[1]<y_max)
                  {
                     node = elem->nodes[node_ymax];
                     cent[0] = node->X();
                     cent[1] = node->Y();
                  }
            if (elem->GetElementType() == MshElemType::TRIANGLE)
		    {
                 	 /// Check whether this point is in this element.
		       	  sub_area[0] = ComputeDetTri(cent, x2, x3);
		          sub_area[1] = ComputeDetTri(cent, x3, x1);
		          sub_area[2] = ComputeDetTri(cent, x1, x2);
		          area =  ComputeDetTri(x1, x2, x3);
                  inarea = fabs(area-sub_area[0]-sub_area[1]-sub_area[2]);

                  /// This point locates within the element
                  if(inarea<tol_a)
                  {
                     /// Use sub_area[k] as shape function
                     for(k=0; k<3; k++)
                        sub_area[k] /= area;

                     l = nn*ncol + mm;
                     cell_data_p[l] = 0.0;

                     for(k=0; k<3; k++)
                        vel_av[k] = 0.;

                     for(j=0; j<nnodes; j++)
                     {
                        n_idx = elem->GetNodeIndex(j);

                        cell_data_p[l] +=  sub_area[j]*GetNodeValue(n_idx,idx);

                        for(k=0; k<3; k++) // k is dimension
                           vel_av[k] += sub_area[j]*GetNodeValue(n_idx,vel_idx[k]);
                     }
                     cell_data_v[l] =  (vel_av[0]*(*elem->tranform_tensor)(0,2)
                        +vel_av[1]*(*elem->tranform_tensor)(1,2)
                                                  // 1000*:  m-->mm
                        +vel_av[2]*(*elem->tranform_tensor)(2,2));

                  }
		    }
		    else if (elem->GetElementType() == MshElemType::QUAD)
		    {
		    	/// Check whether this point is in this element.
    			sub_area[0] = ComputeDetTri(cent, x2, x3);
		        sub_area[1] = ComputeDetTri(cent, x3, x4);
		        sub_area[2] = ComputeDetTri(cent, x4, x1);
	    		sub_area[3] = ComputeDetTri(cent, x1, x2);
		    	area = ComputeDetTri(x1, x2, x3) + ComputeDetTri(x1, x3, x4);
  			inarea = fabs(area-sub_area[0]-sub_area[1]-sub_area[2]-sub_area[3]);
                /// This point locates within the element
                  if(inarea<tol_a)
                  {
                     /// Use sub_area[k] as shape function
                     for(k=0; k<4; k++)
                        sub_area[k] /= area;

                     l = nn*ncol + mm;
                     cell_data_p[l] = 0.0;

                     for(k=0; k<3; k++)
                        vel_av[k] = 0.;

                     for(j=0; j<4; j++) //nnodes: number of nodes in an element
                     {
                        n_idx = elem->GetNodeIndex(j);

                        cell_data_p[l] += sub_area[j]*GetNodeValue(n_idx,idx); //GetNodeValue(n_idx,idx)指第n_idx个节点的idx(此处是Head)

                        cell_data_p[l] = 2 * cell_data_p[l];

                        for(k=0; k<3; k++)
                           vel_av[k] += sub_area[j]*GetNodeValue(n_idx,vel_idx[k]);
                     }
                     cell_data_v[l] =  -(vel_av[0]*(*elem->tranform_tensor)(0,2)
                        +vel_av[1]*(*elem->tranform_tensor)(1,2)
                                                  // 1000*:  m-->mm
                        +vel_av[2]*(*elem->tranform_tensor)(2,2));//通过张量变换求z方向的速度

                       cell_data_v[l] = 2 * cell_data_v[l];
                  }
            }
        }

    }

}

         of_flux.precision(4);
         //of_flux.setf(ios::scientific, ios::floatfield);
         of_primary.precision(2);
         for(irow=0; irow<nrow; irow++)
         {
            for(icol=0; icol<ncol; icol++)
            {
               m = irow*ncol+icol;
               of_flux<<" "<<setw(9)<< cell_data_v[m];
               of_primary<<setw(11)<< cell_data_p[m];
            }
            of_flux<<endl;
            of_primary<<endl;
         }

         cell_data_p.clear();
         cell_data_v.clear();
         of_flux.close();
         of_primary.close();
      }
   }
   //-------------------------------------------------------------------------------------
   if (bc_transient_index.size() == 0)
      return;
   bool valid = false;
   long end_i = 0;
   double t_fac = 0.;
   std::vector<double> node_value;

   for (size_t i = 0; i < bc_transient_index.size(); i++)
   {
      CBoundaryCondition *bc = bc_node[bc_transient_index[i]];
                                                  // TF
      CGLPolyline *polyline = GEOGetPLYByName (bc->getGeoName());
      long start_i = bc_transient_index[i];
      if (i == bc_transient_index.size() - 1)
         end_i = (long) bc_node.size();
      else
         end_i = bc_transient_index[i + 1];
      node_value.resize(end_i - start_i);
      // Piecewise linear distributed.
      for (size_t k = 0; k < (bc->getDistribedBC()).size(); k++)
      {
         for (size_t j = 0; j < polyline->point_vector.size(); j++)
         {
            if ((bc->getPointsWithDistribedBC())[k] == polyline->point_vector[j]->id)
            {
               if (fabs((bc->getDistribedBC())[k]) < MKleinsteZahl)
                  (bc->getDistribedBC())[k] = 1.0e-20;
               polyline->point_vector[j]->setPropert ((bc->getDistribedBC())[k]);
               CFunction* m_fct = FCTGet(bc->getPointsFCTNames()[k]);
               if (m_fct)
                  t_fac = m_fct->GetValue(aktuelle_zeit, &valid);
               else
                  cout
                     << "Warning in CBoundaryConditionsGroup - no FCT data"
                     << endl;
               if (valid)
                  polyline->point_vector[j]->setPropert (polyline->point_vector[j]->getPropert()*t_fac);
               break;
            }
         }
      }
      std::cout << "CRFProcess::UpdateTransientBC() for polyline: " << bc->getGeoName() << std::endl;
      std::cout << "polyline.ibuffer.size(): " << polyline->getIBuffer().size() << std::endl;
      std::cout << "polyline.getOrderedPoints().size(): " << polyline->getOrderedPoints().size() << std::endl;
      InterpolationAlongPolyline(polyline, node_value);
      for (long k = start_i; k < end_i; k++)
      {
         bc_node_value[k]->node_value = node_value[k - start_i];
      }
   }
}


/*************************************************************************
ROCKFLOW - Function: CRFProcess::Add_GEMS_Water_ST
Task: add one more record into the vector Water_ST_vec
Programming:
11/2008 //HS Implementation
**************************************************************************/

void CRFProcess::Add_GEMS_Water_ST(long idx, double val)
{
   Water_ST_GEMS tmp_st;

   tmp_st.index_node = idx;
   tmp_st.water_st_value = val;

   Water_ST_vec.push_back(tmp_st);
}


void CRFProcess::Clean_Water_ST_vec()
{
   Water_ST_vec.clear();
}


void CRFProcess::configMaterialParameters ()
{
   // Output material parameters
   //WW
   const size_t out_vector_size (out_vector.size());

   for (size_t i = 0; i < out_vector_size; i++)
   {
      COutput *out = out_vector[i];
      const size_t size (out->_nod_value_vector.size());
      for (size_t k = 0; k < size; k++)
      {
         if (out->_nod_value_vector[k].find("PERMEABILITY_X1") != string::npos)
         {
            additioanl2ndvar_print = 1;
            break;
         }
      }
      if (additioanl2ndvar_print == 1)
         break;
   }

   for (size_t i = 0; i < out_vector_size; i++)
   {
      COutput *out = out_vector[i];
      const size_t size (out->_nod_value_vector.size());
      for (size_t k = 0; k < size; k++)
      {
         if (out->_nod_value_vector[k].find("POROSITY") != string::npos)
         {
            if (additioanl2ndvar_print > 0)
               additioanl2ndvar_print = 2;
            else
               additioanl2ndvar_print = 3;
         }
         if (additioanl2ndvar_print > 1)
            break;
      }
      if (additioanl2ndvar_print > 1)
         break;
   }

   if (additioanl2ndvar_print > 0)                //WW
   {
      if (additioanl2ndvar_print < 3)
      {
         pcs_secondary_function_name[pcs_number_of_secondary_nvals]
            = "PERMEABILITY_X1";
         pcs_secondary_function_unit[pcs_number_of_secondary_nvals]
            = "1/m^2";
         pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
         pcs_number_of_secondary_nvals++;
         pcs_secondary_function_name[pcs_number_of_secondary_nvals]
            = "PERMEABILITY_Y1";
         pcs_secondary_function_unit[pcs_number_of_secondary_nvals]
            = "1/m^2";
         pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
         pcs_number_of_secondary_nvals++;
         if (max_dim == 2)                        // 3D
         {
            pcs_secondary_function_name[pcs_number_of_secondary_nvals]
               = "PERMEABILITY_Z1";
            pcs_secondary_function_unit[pcs_number_of_secondary_nvals]
               = "1/m^2";
            pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals]
               = 1;
            pcs_number_of_secondary_nvals++;
         }
      }
      if (additioanl2ndvar_print > 1)             //WW
      {
         pcs_secondary_function_name[pcs_number_of_secondary_nvals]
            = "POROSITY";
         pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
         pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
         pcs_number_of_secondary_nvals++;
      }
   }
}


/*************************************************************************
GeoSys - Function:
06/2009 OK Implementation
**************************************************************************/
bool PCSConfig()
{
   bool some_thing_done = false;
   for(int i=0;i<(int)pcs_vector.size();i++)      //OK
   {
      pcs_vector[i]->Config();
      some_thing_done = true;
   }
   return some_thing_done;
}
