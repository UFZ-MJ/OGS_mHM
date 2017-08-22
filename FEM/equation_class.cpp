/**************************************************************************
Task: Linear equation
Programing:
11/2007 WW/
**************************************************************************/
// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#if defined(USE_MPI)
#include <mpi.h>
#include "par_ddc.h"
#endif

#include "Configure.h"
#include "makros.h"

#include <cfloat>
// NEW_EQS To be removed
#ifdef NEW_EQS                                    //1.11.2007 WW
#include <iomanip>

#ifdef LIS                                        // 07.02.2008 PCH
#include "lis.h"
#include <omp.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif

#ifdef MKL
#include "mkl.h"

/* PARDISO prototype. */
//#define PARDISO pardiso_
extern int omp_get_max_threads();
extern int PARDISO
(int *, int *, int *, int *, int *, int *,
double *, int *, int*, int *, int *, int *,
int *, double *, double*, int *);
#endif

#include"rf_num_new.h"
#include"matrix_class.h"
#include"equation_class.h"

std::vector<Math_Group::Linear_EQS *> EQS_Vector;
using namespace std;

//
namespace Math_Group
{
   /**************************************************************************
   Task: Linear equation::Constructor
   Programing:
   10/2007 WW/
   **************************************************************************/
   Linear_EQS::Linear_EQS(const SparseTable &sparse_table,
      const int dof, bool messg):message(messg)
   {
      long i, size;
      A = new CSparseMatrix(sparse_table, dof);
      size = A->Dim();
#if defined(USE_MPI)
      prec_M =NULL;
      x = NULL;
      cpu_time = 0.0;
      // WW
      border_buffer0 = NULL;
      border_buffer1 = NULL;
#else
      x = new double[size];
#endif
      b = new double[size];
      //
      for(i=0; i<size; i++)
      {
#ifndef USE_MPI
         x[i] = 0.;
#endif
         b[i] = 0.;
      }
      iter = 0;
      bNorm = 1.0;
      error = 1.0e10;
      size_global = 0;

      //NW moved below to rf.cpp. Initialization and finalization of LIS solver
      //must be called only once but they were called several times with multiple meshes
#ifdef LIS
      ////	A->Write();
      //
      int argc=0;
      char** argv = NULL;
      // Initialization of the lis solver.
      lis_initialize(&argc, &argv);
#endif
   }
#if defined(USE_MPI)
   /**************************************************************************
   Task: Linear equation::Constructor
   Programing:
   12/2007 WW/
   **************************************************************************/
   Linear_EQS::Linear_EQS(const long size)
   {
      long i;
      A = NULL;
      b = NULL;
      size_global = size;
      x = new double[size];
      //
      for(i=0; i<size; i++)
         x[i] = 0.;
      iter = 0;
      bNorm = 1.0;
      error = 1.0e10;
   }
#endif
   /**************************************************************************
   Task: Linear equation::Destructor
   Programing:
   10/2007 WW/
   **************************************************************************/
   Linear_EQS::~Linear_EQS()
   {
      if(A) delete A;
      if(x) delete [] x;
      if(b) delete [] b;
      //
      A = NULL;
      x = NULL;
      b = NULL;

      //NW moved to rf.cpp
#ifdef LIS
      lis_finalize();
#endif
      //
   }
   /**************************************************************************
   Task: Linear equation::
   Programing:
   10/2007 WW/
   **************************************************************************/
   void Linear_EQS::ConfigNumerics(CNumerics *m_num, const long n)
   {
      (void)n;
      int i, nbuffer = 0;                         // Number of temperary float arrays
      precond_type = m_num->ls_precond;
      solver_type = m_num->ls_method;
      switch(solver_type)
      {
         case 1:
            solver_name = "Gauss";
            break;
         case 2:
            solver_name = "BiCGSTab";
            nbuffer = 8;
            break;
         case 3:
            solver_name = "BiCG";
            nbuffer = 8;                          //20.10.2010. WW
            break;
         case 4:
            solver_name = "QMRCGStab";
            break;
         case 5:
            solver_name = "CG";
            nbuffer = 3;
            break;
         case 6:
            solver_name = "CGNR";
            break;
         case 7:
            solver_name = "CGS";
            nbuffer = 9;
            break;
         case 8:
            solver_name = "Richardson";
            break;
         case 9:
            solver_name = "JOR";
            break;
         case 10:
            solver_name = "SOR";
            break;
         case 11:
            solver_name = "AMG1R5";
            break;
         case 12:
            solver_name = "UMF";
            break;
      }
      // Buffer
      /*
      #if defined(USE_MPI)
      long size = A->Dim()+A->Dof()*dom->BSize();
      #else
      */
      long size = A->Dim();
      //#endif
      for(i=0; i<nbuffer; i++)
      {
         double *new_array = new double[size];
         f_buffer.push_back(new_array);
      }
#if defined(USE_MPI)
      if(!x)
         x = new double[size];
#if defined(NEW_BREDUCE)
      // For concatenate border entries
      double *catb = NULL;
      f_buffer.push_back(catb);
#endif
      // Buffer for a global array
      double *x_g_buffer = new double[n];
      f_buffer.push_back(x_g_buffer);
      // Buffer for border array
      // Will be removed if topo is ready
      border_buffer0 = new double[A->Dof()*dom->BSize()];
      border_buffer1 = new double[A->Dof()*dom->BSize()];
      //
#endif
      //---------------------------------------------
      switch(precond_type)
      {
         case 1:
            precond_name = "Jacobi";
#if defined(USE_MPI)
            prec_M = new double[size];
#endif
            break;
         case 100:
            // precond_name = "ILU"; break;
            // If ILU is ready, remove follows
            // ----------------------------------------------
            precond_name = "ILU not available. Use Jacobi";
            precond_type = 1;
#if defined(USE_MPI)
            prec_M = new double[size];
#endif
            // ----------------------------------------------
            break;
         default:
            precond_name = "No preconditioner"; break;
      }
      //
      //
      max_iter = m_num->ls_max_iterations;
      tol = m_num->ls_error_tolerance;
      //
   }
   /**************************************************************************
   Task: Linear equation::Alocate memory for solver
   Programing:
   11/2007 WW/
   **************************************************************************/
   void Linear_EQS::Initialize()
   {
      (*A) = 0.;
      for(long i=0; i<A->Dim(); i++)
         b[i] = 0.;
      error = 1.0e10;
   }
   /**************************************************************************
   Task: Linear equation::Alocate memory for solver
   Programing:
   11/2007 WW/
   **************************************************************************/
   void Linear_EQS::Clean()
   {
#if defined(USE_MPI)
      double cpu_time_local = -MPI_Wtime();
#endif
      for(int i=0; i<(int)f_buffer.size(); i++)
      {
         if(f_buffer[i]) delete [] f_buffer[i];
         f_buffer[i] = NULL;
      }
      f_buffer.clear();
#if defined(USE_MPI)
      if(prec_M) delete [] prec_M;
      //
      if(border_buffer0) delete [] border_buffer0;
      border_buffer0 = NULL;
      if(border_buffer1) delete [] border_buffer1;
      border_buffer1 = NULL;
      //
      prec_M = NULL;
      cpu_time_local += MPI_Wtime();
      cpu_time += cpu_time_local;
#endif
   }
   /**************************************************************************
   Task: Linear equation::Write
   Programing:
   11/2007 WW/
   **************************************************************************/
   void Linear_EQS::Write(std::ostream &os)
   {
      A->Write(os);
      //
      os<<" b ( RHS): " <<endl;
      os.width(10);
      os.precision(6);
      //
      for(long i=0; i<A->Dim(); i++)
         os<<setw(10)<<i<<" "<<setw(15)<<b[i]<<endl;
   }
   /**************************************************************************
   Task: Linear equation::Write
   Programing:
   07/2010 NW
   **************************************************************************/
   void Linear_EQS::WriteRHS(ostream &os)
   {
      os.width(10);
      os.precision(6);
      //
      for(long i=0; i<A->Dim(); i++)
         os<<setw(15)<<b[i]<<endl;
   }
   /**************************************************************************
   Task: Linear equation::Write
   Programing:
   07/2010 NW
   **************************************************************************/
   void Linear_EQS::WriteX(ostream &os)
   {
      os.width(10);
      os.precision(6);
      //
      for(long i=0; i<A->Dim(); i++)
         os<<setw(15)<<x[i]<<endl;
   }
   /**************************************************************************
   Task: Linear equation::Solver
   Programing:

   PARDISO openmp-paralle direct solver: 805

   LIS matrix solver options
   CG -i {cg|1}
   BiCG -i {bicg|2}
   CGS -i {cgs|3}
   BiCGSTAB -i {bicgstab|4}
   BiCGSTAB(l) -i {bicgstabl|5} -ell [2] Value for l
   GPBiCG -i {gpbicg|6}
   TFQMR -i {tfqmr|7}
   Orthomin(m) -i {orthomin|8} -restart [40] Value for Restart m
   GMRES(m) -i {gmres|9} -restart [40] Value for Restart m
   Jacobi -i {jacobi|10}
   Gauss-Seidel -i {gs|11}
   SOR -i {sor|12} -omega [1.9] Value for Relaxation Coefficient  (0 <  < 2)
   BiCGSafe -i {bicgsafe|13}
   CR -i {cr|14}
   BiCR -i {bicr|15}
   CRS -i {crs|16}
   BiCRSTAB -i {bicrstab|17}
   GPBiCR -i {gpbicr|18}
   BiCRSafe -i {bicrsafe|19}
   FGMRES(m) -i {fgmres|20} -restart [40] Value for Restart m
   IDR(s) -i {idrs|21} -restart [40] Value for Restart s

   Preconditioner Option Auxiliary Option
   None -p {none|0}
   Jacobi -p {jacobi|1}
   ILU(k) -p {ilu|2} -ilu_fill [0] Fill level k
   SSOR -p {ssor|3} -ssor_w [1.0] Relaxation Coefficient  (0 <  < 2)
   Hybrid -p {hybrid|4} -hybrid_i [sor] Iterative method
   -hybrid_maxiter [25] Maximum number of iterations
   -hybrid_tol [1.0e-3] Convergence criteria
   -hybrid_w [1.5] Relaxation Coefficient  for
   the SOR method (0 <  < 2)
   -hybrid_ell [2] Value for l of the BiCGSTAB(l) method
   -hybrid_restart [40] Restart values for GMRES and Orthomin
   I+S -p {is|5} -is_alpha [1.0] Parameter ?for preconditioner
   of a I + ?(m) type
   -is_m [3] Parameter m for preconditioner
   of a I + ?(m) type
   SAINV -p {sainv|6} -sainv_drop [0.05] Drop criteria
   SA-AMG -p {saamg|7} -saamg_unsym [false] Selection of asymmetric version
   Crout ILU -p {iluc|8} -iluc_drop [0.05] Drop criteria
   -iluc_rate [5.0] Ratio of Maximum fill-in
   ILUT -p {ilut|9} -ilut_drop [0.05] Drop criteria
   -ilut_rate [5.0] Ratio of Maximum fill-in
   additive Schwarz -adds true -adds_iter [1] Number of iterations

   02/2008 PCH OpenMP parallelization by LIS
   03/2009 PCH Solver type and precondition options added for .num file
   **************************************************************************/
#if defined(USE_MPI)
   int Linear_EQS::Solver(double *xg, const long n)
   {
      //
      double cpu_time_local = -MPI_Wtime();
      iter = 0;
      ComputePreconditioner();
      size_global = n;
      switch(solver_type)
      {
         case 2:
            iter = BiCGStab(xg, n);
            break;
         case 3:
            iter = BiCG(xg, n);
            break;
         case 5:
            iter = CG(xg, n);
            break;
         case 7:
            iter = CGS(xg, n);
            break;
      }
      cpu_time_local += MPI_Wtime();
      cpu_time += cpu_time_local;
      return iter;
   }

#else
#ifdef LIS                                     // PCH 02.2008
   //NW 01.2010 Use configuration in NUM file
   int Linear_EQS::Solver(CNumerics *num)
   {

      // Check the openmp solver type iterative and directive
      CNumerics* m_num;
      if (num!=NULL)
         m_num = num;                             //NW
      else
         m_num = num_vector[0];

      if(m_num->ls_method==805)                   // Then, PARDISO parallel direct solver
      {
#ifdef MKL                               // PCH 10.03.2009: Requires the system platform where Math Kernel Library is properly configured.
         int i, iter, ierr;
         // Assembling the matrix
         // Establishing CRS type matrix from GeoSys Matrix data storage type
         int nonzero=A->nnz();
         int numOfNode=A->Size()*A->Dof();
         double *value;
         value = new double [nonzero];
         ierr = A->GetCRSValue(value);
         int* ptr=NULL; int *index=NULL;
         ptr= (int *)malloc((numOfNode+1)*sizeof( int));
         index = (int *)malloc((nonzero)*sizeof( int));

         // Reindexing ptr according to Fortran-based PARDISO
         for(i=0; i < numOfNode; ++i)
            ptr[i] = A->ptr[i]+1;
         //ptr needs one more storage
         ptr[i] = A->ptr[i]+1;
         // Reindexing index according to Fortran-based PARDISO
         // and zonzero of Matrix A
         for(i=0; i < nonzero; ++i)
            index[i] = A->col_idx[i]+1;

         int mtype = 11;                          /* Real unsymmetric matrix */
         int nrhs = 1;                            /* Number of right hand sides. */
         /* Internal solver memory pointer pt, */
         /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
         /* or void *pt[64] should be OK on both architectures */
         void *pt[64];
         /* Pardiso control parameters.*/
         int iparm[64];
         int maxfct, mnum, phase, error, msglvl;

         /* Auxiliary variables.*/
         double ddum;                             /* Double dummy */
         int idum;                                /* Integer dummy. */

         /* --------------------------------------------------------------------*/
         /* .. Setup Pardiso control parameters.*/
         /* --------------------------------------------------------------------*/
         for (i = 0; i < 64; i++)
         {
            iparm[i] = 0;
         }
         iparm[0] = 1;                            /* No solver default */
         iparm[1] = 2;                            /* Fill-in reordering from METIS */
         /* Numbers of processors, value of MKL_NUM_THREADS */
         iparm[2] = mkl_get_max_threads();
         iparm[3] = 0;                            /* No iterative-direct algorithm */
         iparm[4] = 0;                            /* No user fill-in reducing permutation */
         iparm[5] = 0;                            /* Write solution into x */
         iparm[6] = 0;                            /* Not in use */
         iparm[7] = 2;                            /* Max numbers of iterative refinement steps */
         iparm[8] = 0;                            /* Not in use */
         iparm[9] = 13;                           /* Perturb the pivot elements with 1E-13 */
         iparm[10] = 1;                           /* Use nonsymmetric permutation and scaling MPS */
         iparm[11] = 0;                           /* Not in use */
         iparm[12] = 0;                           /* Not in use */
         iparm[13] = 0;                           /* Output: Number of perturbed pivots */
         iparm[14] = 0;                           /* Not in use */
         iparm[15] = 0;                           /* Not in use */
         iparm[16] = 0;                           /* Not in use */
         iparm[17] = -1;                          /* Output: Number of nonzeros in the factor LU */
         iparm[18] = -1;                          /* Output: Mflops for LU factorization */
         iparm[19] = 0;                           /* Output: Numbers of CG Iterations */
         maxfct = 1;                              /* Maximum number of numerical factorizations. */
         mnum = 1;                                /* Which factorization to use. */
         msglvl = 0;                              /* Print statistical information in file */
         error = 0;                               /* Initialize error flag */

         /* --------------------------------------------------------------------*/
         /* .. Initialize the internal solver memory pointer. This is only */
         /* necessary for the FIRST call of the PARDISO solver. */
         /* --------------------------------------------------------------------*/
         for (i = 0; i < 64; i++)
         {
            pt[i] = 0;
         }

         /* --------------------------------------------------------------------*/
         /* .. Reordering and Symbolic Factorization. This step also allocates */
         /* all memory that is necessary for the factorization. */
         /* --------------------------------------------------------------------*/
         phase = 11;
         PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &numOfNode, value, ptr, index , &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

         if (error != 0)
         {
            printf("\nERROR during symbolic factorization: %d", error);
            exit(1);
         }

         /* --------------------------------------------------------------------*/
         /* .. Numerical factorization.*/
         /* --------------------------------------------------------------------*/
         phase = 22;
         PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &numOfNode, value, ptr, index, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
         if (error != 0)
         {
            printf("\nERROR during numerical factorization: %d", error);
            exit(2);
         }

         /* --------------------------------------------------------------------*/
         /* .. Back substitution and iterative refinement. */
         /* --------------------------------------------------------------------*/
         phase = 33;
         iparm[7] = 2;                            /* Max numbers of iterative refinement steps. */

         /* Set right hand side to one. */

         PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &numOfNode, value, ptr, index , &idum, &nrhs,
            iparm, &msglvl, b, x, &error);
         if (error != 0)
         {
            printf("\nERROR during solution: %d", error);
            exit(3);
         }

         phase = -1;                              /* Release internal memory. */
         PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &numOfNode, value, ptr, index, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

         // Releasing the local memory
         delete [] value; free(ptr); free(index);
         //		MKL_FreeBuffers();
#endif
      }
      else                                        // LIS parallel solver
      {
         int i, iter, ierr, size;
         // Fix for the fluid_momentum Dof
         size = A->Size()*A->Dof();

         // Assembling the matrix
         // Establishing CRS type matrix from GeoSys Matrix data storage type
         int nonzero=A->nnz();
         double *value;
         value = new double [nonzero];
         ierr = A->GetCRSValue(value);

         // Creating a matrix.
         ierr = lis_matrix_create(0,&AA);
         ierr = lis_matrix_set_type(AA,LIS_MATRIX_CRS);
         ierr = lis_matrix_set_size(AA,0,size);

         // Matrix solver and Precondition can be handled better way.
         char solver_options[MAX_ZEILE], tol_option[MAX_ZEILE];
         sprintf(solver_options, "-i %d -p %d",m_num->ls_method, m_num->ls_precond);
         // tolerance and other setting parameters are same
                                                  //NW add max iteration counts
         sprintf(tol_option, "-tol %e -maxiter %d",m_num->ls_error_tolerance, m_num->ls_max_iterations);

         ierr = lis_matrix_set_crs(nonzero,A->ptr,A->col_idx, value,AA);
         ierr = lis_matrix_assemble(AA);

         // Assemble the vector, b, x
         //OK411 int iflag = 0;
         ierr = lis_vector_duplicate(AA,&bb);
         ierr = lis_vector_duplicate(AA,&xx);
#pragma omp parallel for
         for(i=0; i < size; ++i)
         {
            ierr = lis_vector_set_value(LIS_INS_VALUE,i,x[i],xx);
            ierr = lis_vector_set_value(LIS_INS_VALUE,i,b[i],bb);
         }

         // Create solver
         ierr = lis_solver_create(&solver);

         ierr = lis_solver_set_option(solver_options,solver);
         ierr = lis_solver_set_option(tol_option,solver);
         ierr = lis_solve(AA,bb,xx,solver);
         ierr = lis_solver_get_iters(solver,&iter);
         //NW
         printf("\t iteration: %d/%d\n", iter, m_num->ls_max_iterations);
         double resid = 0.0;
         ierr = lis_solver_get_residualnorm(solver,&resid);
         printf("\t residuals: %e\n", resid);
         //	lis_vector_print(xx);
         //	lis_vector_print(bb);

         // Update the solution (answer) into the x vector
#pragma omp parallel for
         for(i=0; i<size; ++i)
         {
            lis_vector_get_value(xx,i,&(x[i]));
         }

         // Clear memory
         delete [] value;
         //	lis_matrix_destroy(AA);
         lis_vector_destroy(bb);
         lis_vector_destroy(xx);
         lis_solver_destroy(solver);
      }

      return -1;                                  // This right now is meaningless.
   }
#else
   int Linear_EQS::Solver()
   {
      //
      iter = 0;
      ComputePreconditioner();
      switch(solver_type)
      {
         case 1:
            return Gauss();
         case 2:
            iter= BiCGStab();
            return iter;                          //kg44 only to make sure here is iter returned
         case 3:
            return BiCG();
         case 4:
            return QMRCGStab();
         case 5:
            return CG();
         case 6:
            return CGNR();
         case 7:
            return CGS();
         case 8:
            return Richardson();
         case 9:
            return JOR();
         case 10:
            return SOR();
         case 11:
            return AMG1R5();
         case 12:
            return UMF();
      }
      return -1;
   }
#endif
#endif
   // Preconditioners
   /**************************************************************************
   Task: Preconditioners
   Programing:
   10/2007 WW
   **************************************************************************/
   void Linear_EQS::ComputePreconditioner()
   {
      switch(precond_type)
      {
         case 1:
#if defined(USE_MPI)
            ComputePreconditioner_Jacobi();
#endif
            return;
         case 100:
            ComputePreconditioner_ILU();
            return;
         default:
            return;
      }
   }
   /**************************************************************************
   Task: Linear equation::SetKnownXi
         Configure equation system when one entry of the vector of
         unknown is given
   Programing:
   10/2007 WW/
   **************************************************************************/
   void Linear_EQS::SetKnownX_i(const long i, const double x_i)
   {
      A->Diagonize(i, x_i, b);
   }
   /**************************************************************************
   Task: Linear equation::Preconditioner
   Programing:
   08/2007 WW/
   **************************************************************************/
   void Linear_EQS::Precond(double *vec_s, double *vec_r)
   {
      bool pre = true;
      switch(precond_type)
      {
         case 1:
#if defined(USE_MPI)
            Precond_Jacobi(vec_s, vec_r);
#else
            A->Precond_Jacobi(vec_s, vec_r);
#endif
            break;
         case 100:
            pre = false;                          //A->Precond_ILU(vec_s, vec_r);
            break;
         default:
            pre = false;                          //A->Precond_ILU(vec_s, vec_r);
            break;
      }
      if(!pre)
      {
         for(long i=0; i<A->Dim(); i++)
            vec_r[i] = vec_s[i];
      }
   }
   /**************************************************************************
   Task: Linear equation:: M^T x
     Transpose of preconditioner times a vector
   Programing:
   02/2010 WW/
   **************************************************************************/
   void Linear_EQS::TransPrecond(double *vec_s, double *vec_r)
   {
      Precond(vec_s, vec_r);
   }
   /*\!
    ********************************************************************
      Dot production of two vectors
      Programm:
      10/2007 WW
      12/2007 WW  Parallel
   ********************************************************************/
   double Linear_EQS::dot (const double *xx,  const double *yy)
   {
      long i;
      double val = 0.;
#if defined(USE_MPI)
      long size = A->Dim();
      double val_i = dom->Dot_Interior(xx,  yy);
      val_i += dom->Dot_Border_Vec(xx,  yy);
      //
      MPI_Allreduce(&val_i, &val, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
      for(i=0; i<A->Dim(); i++)
         val += xx[i]*yy[i];
#endif
      return val;
   }
   /*\!
    ********************************************************************
      Dot production of two vectors
      Programm:
      01/2008 WW
   ********************************************************************/
   double  Linear_EQS::NormX()
   {
#if defined(USE_MPI)
      return sqrt(dot(x, x, size_global ));
#else
      return sqrt(dot(x, x ));
#endif
   }
   //
#if defined(USE_MPI)
   /*\!
    ********************************************************************
      Dot production of two vectors
      Programm:
      12/2007 WW
   ********************************************************************/
   double Linear_EQS::dot (const double *xx,  const double *yy, const long n)
   {
      long i;
      double val = 0.;
      for(i=0; i<n; i++)
         val += xx[i]*yy[i];
      return val;
   }
   /*\!
    ********************************************************************
      Dot production of two vectors
      Programm:
      12/2007 WW
   ********************************************************************/
   inline void Linear_EQS::MatrixMulitVec(double *xx,  double *yy)
   {
      long i, size = A->Dim();
      //
      A->multiVec(xx, yy);
#if defined(NEW_BREDUCE)
      dom->ReduceBorderV(yy);
#else
      dom->Local2Border(yy, border_buffer0);
      MPI_Allreduce(border_buffer0, border_buffer1, A->Dof()*dom->BSize(),
         MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Border2Local(border_buffer1, yy);
#endif
   }
   /*!
    ********************************************************************
      Dot production of two vectors
      Programm:
      12/2007 WW
      02/2010 WW Revise
   ********************************************************************/
   inline void Linear_EQS::TransMatrixMulitVec(double *xx,  double *yy)
   {
      //
      A->Trans_MultiVec(xx, yy);
#if defined(NEW_BREDUCE)
      dom->ReduceBorderV(yy);
#else
      dom->Local2Border(yy, border_buffer0);
      MPI_Allreduce(border_buffer0, border_buffer1, A->Dof()*dom->BSize(),
         MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Border2Local(border_buffer1, yy);
#endif
   }
#endif
   /*!
    ********************************************************************
      ConvergeTest
      Programm:
      09/2007 WW
   ********************************************************************/
   void Linear_EQS::Message()
   {
#ifdef USE_MPI
      if(myrank>0)
         return;
#endif
      if (!message) return;
      cout.width(10);
      cout.precision(3);
      cout.setf(ios::scientific);
      //
      //system("color 0B");
      cout<<"\n================================================\n";
      cout<<"Linear solver "<<solver_name<<" with "<<precond_name<<":\n";
      //cout<<"\n------------------------------------------------ \n";
      cout<<"Iterations |";
      cout<<"Max Iters  |";
      cout<<"Norm of b  |";
      cout<<"Error      |\n";
      //cout<<"\n------------------------------------------------ \n";
      cout<<setw(11)<<iter<<"|"<<setw(11)<<max_iter<< "|"
         <<setw(11)<<bNorm<< "|"<<setw(11)<<error<<"|\n";
      cout<<"================================================\n";
      if (iter==max_iter) cout << " Maximum iteration reached !!! \n";
      cout.flush();
      //
   }
   /*\!
    ********************************************************************
      Check if the norm of b is samll enough for convengence.
      normb_new is given to bNorm;
      Programm:
      09/2007 WW
   ********************************************************************/
   inline bool Linear_EQS::CheckNormRHS(const double normb_new)
   {
      if(bNorm>0.0)
      {
         if((normb_new/bNorm)<tol)
         {
            error = normb_new/bNorm;
            bNorm = normb_new;
            Message();
            return true;
         }
      }
      bNorm = normb_new;
      if(bNorm<DBL_MIN)
      {
         error = 0.;
         Message();
         return true;
      }
      return false;
   }
#ifndef USE_MPI
   /**************************************************************************
   Task: Linear equation::CG
   Programing:
   11/2007 WW/
   **************************************************************************/
   int Linear_EQS::CG()
   {
      long i, size;
      double rrM1;
      double *p, *r, *s;
      //
      size = A->Dim();
      p = f_buffer[0];
      r = f_buffer[1];
      s = f_buffer[2];
      //
      double bNorm_new = Norm(b);
      // Check if the norm of b is samll enough for convengence
      if(CheckNormRHS(bNorm_new))
         return 0;
      //
      // r0 = b-Ax
      A->multiVec(x,s);
      for(i=0; i<size; i++)
         r[i] = b[i]-s[i];
      //
      // Preconditioning: M^{-1}r
      Precond(r, s);
      for(i=0; i<size; i++)
         p[i] = s[i];
      // Check the convergence
      if ((error = Norm(r)/bNorm) < tol)
      {
         Message();
         return 1;
      }
      //
      double rr = dot(r , s);
      //
      for (iter=1; iter<=max_iter; ++iter)
      {
         A->multiVec(p, s);
         const double alpha = rr / dot(p,s);
         // Update
         for(i=0; i<size; i++)
         {
            x[i] += alpha*p[i];
            r[i] -= alpha*s[i];
         }
         if ((error = Norm(r)/bNorm) < tol)
         {
            Message();
            return iter <= max_iter;
         }
         //
         Precond(r, s);
         //
         rrM1 = rr;
         rr   = dot(s, r);
         const double beta = rr / rrM1;
         for(i=0; i<size; i++)
            p[i] = s[i] + beta*p[i];
      }
      //
      Message();
      return iter <= max_iter;
   }
   /**************************************************************************
   Task: Linear equation::BiCG
   Programing:
   10/2010 WW/
   **************************************************************************/
   int Linear_EQS::BiCG()
   {
      long i, size;
      double rho1, rho2 = 0., alpha, beta;
      double *z, *zt, *p, *pt, *q, *qt, *r, *rt;
      //
      size = A->Dim();
      z = f_buffer[0];
      zt = f_buffer[1];
      p = f_buffer[2];
      pt = f_buffer[3];
      q = f_buffer[4];
      qt = f_buffer[5];
      r = f_buffer[6];
      rt = f_buffer[7];
      //
      double bNorm_new = Norm(b);
      // Check if the norm of b is samll enough for convengence
      if(CheckNormRHS(bNorm_new))
         return 0;
      //
      // r0 = b-Ax
      A->multiVec(x,rt);
      for(i=0; i<size; i++)
      {
         r[i] = b[i]-rt[i];
         rt[i] = r[i];
      }
      //
      // Check the convergence
      if ((error = Norm(r)/bNorm) < tol)
      {
         Message();
         return 1;
      }
      //
      //
      for (iter=1; iter<=max_iter; ++iter)
      {
         Precond(r, z);
         TransPrecond(rt, zt);
         rho1 = dot(z, rt);
         //
         if (fabs(rho1) < DBL_MIN)
         {
            Message();
            return iter <= max_iter;
         }
         //
         if(iter == 1)
         {
            for(i=0; i<size; i++)
            {
               p[i] = z[i];
               pt[i] = zt[i];
            }
         }
         else
         {
            beta = rho1/rho2;
            for(i=0; i<size; i++)
            {
               p[i] = z[i] + beta*p[i];
               pt[i] = zt[i] + beta*pt[i];
            }
         }
         //
         A->multiVec(p, q);
         A->Trans_MultiVec(pt, qt);
         alpha = rho1/dot(pt,q);
         //
         for(i=0; i<size; i++)
         {
            x[i] += alpha*p[i];
            r[i] -= alpha*q[i];
            rt[i] -= alpha*qt[i];
         }
         //
         rho2 = rho1;
         if ((error = Norm(r)/bNorm) < tol)
         {
            Message();
            return iter <= max_iter;
         }
         //
      }
      //
      Message();
      return iter <= max_iter;
   }

   /*************************************************************************
   GeoSys-Function:
   Task: BiCGStab solver
   Programming:
   10/2007 WW
   **************************************************************************/
   int Linear_EQS::BiCGStab()
   {
      long i, size;
      double rho_0, rho_1, alpha, beta, omega, tt=0., norm_r = 0.;
      double *r0, *r, *s, *s_h, *t, *v, *p, *p_h;
      //
      size = A->Dim();
      r0 = f_buffer[0];
      r = f_buffer[1];
      s = f_buffer[2];
      s_h = f_buffer[3];
      t = f_buffer[4];
      v = f_buffer[5];
      p = f_buffer[6];
      p_h = f_buffer[7];
      //
      rho_0 = alpha = omega = 1.0;
      //
      double bNorm_new = Norm(b);
      // Check if the norm of b is small enough for convengence
      if(CheckNormRHS(bNorm_new))
         return 0;
      //
      A->multiVec(x,s);                           // s as buffer
      for(i=0; i<size; i++)
      {
         r0[i] = b[i]-s[i];                       // r = b-Ax
         r[i] = r0[i];
         v[i] = 0.;
         p[i] = 0.;
      }
      if ((error = Norm(r)/bNorm) < tol)
      {
         Message();
         return 0;
      }
      //
      for (iter = 1; iter <= max_iter; iter++)
      {
         rho_1 = dot(r0, r);
         if (fabs(rho_1) < DBL_MIN)               // DBL_EPSILON
         {
            Message();
            return 0;
         }
         if (iter == 1)
         {
            for(i=0; i<size; i++)
               p[i] = r[i];
         }
         else
         {
            beta = (rho_1/rho_0) * (alpha/omega);
            for(i=0; i<size; i++)
               p[i] = r[i] + beta * (p[i] - omega * v[i]);
         }
         // Preconditioner
         Precond(p, p_h);
         // A M^{-1}p-->v
         A->multiVec(p_h, v);
         //
         alpha = rho_1 / dot(r0, v);
         //
         for(i=0; i<size; i++)
            s[i] = r[i] - alpha * v[i];
         if ((error = Norm(s)/bNorm) < tol)
         {
            for(i=0; i<size; i++)
               x[i] += alpha * p_h[i];
            Message();
            return iter;
         }
         //  M^{-1}s,
         Precond(s, s_h);
         // A* M^{-1}s
         A->multiVec(s_h, t);
         //
         tt = dot(t,t);
         if(tt>DBL_MIN)
            omega = dot(t,s) / tt;
         else
            omega = 1.0;
         // Update solution
         for(i=0; i<size; i++)
         {
            x[i] += alpha * p_h[i] + omega * s_h[i];
            r[i] = s[i] - omega * t[i];
         }
         rho_0 = rho_1;
         //
         norm_r = Norm(r);
         if ((error = norm_r / bNorm) < tol)
         {
            Message();
            return iter;
         }
         if (fabs(omega) < DBL_MIN)
         {
            error = norm_r / bNorm;
            Message();
            return iter;
         }
      }
      //
      Message();
      //
      return iter;
   }

   /*************************************************************************
   GeoSys-Function:
   Task: CGS solver
   Programming:
   11/2007 WW
   **************************************************************************/
   int Linear_EQS::CGS()
   {
      long i, size;
      double rho_1, rho_2, alpha, beta;
      double *r0, *r, *p, *p_h, *q, *q_h, *v, *u, *u_h;
      //
      size = A->Dim();
      r0 = f_buffer[0];
      r = f_buffer[1];
      p = f_buffer[2];
      p_h = f_buffer[3];
      q = f_buffer[4];
      q_h = f_buffer[5];
      v = f_buffer[6];
      u = f_buffer[7];
      u_h = f_buffer[8];
      //
      rho_1 = rho_2 = 1.0;
      //
      double bNorm_new = Norm(b);
      // Check if the norm of b is samll enough for convengence
      if(CheckNormRHS(bNorm_new))
         return 0;
      //
      A->multiVec(x,v);                           // v as buffer
      for(i=0; i<size; i++)
      {
         r0[i] = b[i]-v[i];                       // r = b-Ax
         r[i] = r0[i];
         v[i] = 0.;
      }
      if ((error = Norm(r)/bNorm) < tol)
      {
         Message();
         return 0;
      }
      //
      for (iter = 1; iter <= max_iter; iter++)
      {
         rho_1 = dot(r0, r);
         if (fabs(rho_1) < DBL_MIN)               //  DBL_EPSILON
         {
            Message();
            return 0;
         }
         if (iter == 1)
         {
            for(i=0; i<size; i++)
               p[i] = u[i] = r[i];
         }
         else
         {
            beta = rho_1/rho_2;
            for(i=0; i<size; i++)
            {
               u[i] = r[i] + beta * q[i];
               p[i] = u[i] + beta * (q[i] + beta * p[i]);
            }
         }
         // Preconditioner
         Precond(p, p_h);
         // A M^{-1}p-->v
         A->multiVec(p_h, v);
         //
         alpha = rho_1 / dot(r0, v);
         //
         for(i=0; i<size; i++)
         {
            q[i] = u[i] - alpha * v[i];
            q_h[i] = u[i] + q[i];
         }
         // Preconditioner
         Precond(q_h, u_h);
         for(i=0; i<size; i++)
            x[i] += alpha * u_h[i];
         //
         A->multiVec(u_h, q_h);
         //
         for(i=0; i<size; i++)
            r[i] -= alpha * q_h[i];
         rho_2 = rho_1;
         if ((error = Norm(r) / bNorm) < tol)
         {
            Message();
            return iter <= max_iter;
         }
      }
      //
      Message();
      //
      return iter <= max_iter;
   }
#endif                                         // If not defined USE_MPI
   //
   //------------------------------------------------------------------------
#if defined(USE_MPI)
   /*************************************************************************
   GeoSys-Function:
   Task: Parallel preconditioner
   Programming:
   12/2007 WW
   **************************************************************************/
   void Linear_EQS::ComputePreconditioner_Jacobi()
   {
      //
      A->DiagonalEntries(prec_M);
      //
#if defined(NEW_BREDUCE)
      dom->ReduceBorderV(prec_M);
#else
      dom->Local2Border(prec_M, border_buffer0);  // to buffer
      MPI_Allreduce(border_buffer0, border_buffer1, A->Dof()*dom->BSize(), MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
      dom->Border2Local(border_buffer1, prec_M);
   }
   /*************************************************************************
   GeoSys-Function:
   Task: Parallel preconditioner, inverse
   Programming:
   12/2007 WW
   **************************************************************************/
   void Linear_EQS::Precond_Jacobi(const double *vec_s, double *vec_r)
   {
      double val;
      for(long i=0; i<A->Dim(); i++)
      {
         val = prec_M[i];
         if(val>DBL_EPSILON)
            vec_r[i] = vec_s[i]/val;
         else
            vec_r[i] = vec_s[i];
      }
   }
#define TEST_MPII
   /**************************************************************************
   Task: Linear equation::CG
   Programing:
   01/2008 WW/
   **************************************************************************/
   int Linear_EQS::CG(double *xg, const long n)
   {
      long i, size, size_b, size_t;
      double rrM1;
      double *p, *r, *s;
      double *p_b, *r_b, *s_b;
      double *x_b;
      //
      size = A->Dim();
      size_b = dom->BSize()*A->Dof();
      //size_t = size + size_b; //[0, size_i): internal; [size_i, size_i+size_b): border.
      //
      p = f_buffer[0];
      r = f_buffer[1];
      s = f_buffer[2];
      //

      //*** Norm b
      double bNorm_new;
      double buff_fl = dom->Dot_Interior(b, b);
      MPI_Allreduce(&buff_fl, &bNorm_new, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Local2Border(b, border_buffer0);       // p_b s_b as buffer
      MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                                                  // (rhs on border)
      buff_fl = bNorm_new + dot(border_buffer1, border_buffer1, size_b);
      bNorm_new = sqrt(buff_fl);
      // Check if the norm of b is samll enough for convengence
      if(CheckNormRHS(bNorm_new))
         return 0;
      //*** r = b-Ax
      //    A*x
      dom->Global2Local(xg, x, n);
      A->multiVec(x,s);                           // s as buffer
      for(i=0; i<size; i++)
         r[i] = b[i]-s[i];                        // r = b-Ax
      //   Collect border r
      dom->Local2Border(r, border_buffer0);       //
      MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Border2Local(border_buffer1, r);
      //
      // Preconditioning: M^{-1}r
      Precond(r, s);
      for(i=0; i<size; i++)
         p[i] = s[i];
      // Check the convergence
      if ((error = Norm(r)/bNorm) < tol)
      {
         Message();
         return 1;
      }
      //
      double rr = dot(r , s);
      //
      for (iter=1; iter<=max_iter; ++iter)
      {
         MatrixMulitVec(p, s);
         const double alpha = rr / dot(p,s);
         // Update
         for(i=0; i<size; i++)
         {
            x[i] += alpha*p[i];
            r[i] -= alpha*s[i];
         }
         if ((error = Norm(r)/bNorm) < tol)
            break;
         // Preconditioner
         Precond(r, s);
         //
         rrM1 = rr;
         rr   = dot(s, r);
         //
         const double beta = rr / rrM1;
         for(i=0; i<size; i++)
            p[i] = s[i] + beta*p[i];

      }
      //
      // concancert internal x
      dom->CatInnerX(xg, x, n);
      //
      Message();
      return iter <= max_iter;
   }

   /*************************************************************************
   GeoSys-Function:
   Task: CGS solver
   Programming:
   02/2008 WW
   **************************************************************************/
   int Linear_EQS::CGS(double *xg, const long n)
   {
      long i, size, size_b, size_t;
      double rho_1, rho_2, alpha, beta;
      double *r0, *r, *p, *p_h, *q, *q_h, *v, *u, *u_h;
      double *r0_b, *r_b, *p_b, *p_h_b, *q_b, *q_h_b, *v_b, *u_b, *u_h_b;
      double *x_b;
      //
      size = A->Dim();
      size_b = dom->BSize()*A->Dof();
      //size_t = size + size_b; //[0, size_i): internal; [size_i, size_i+size_b): border.
      r0 = f_buffer[0];
      r = f_buffer[1];
      p = f_buffer[2];
      p_h = f_buffer[3];
      q = f_buffer[4];
      q_h = f_buffer[5];
      v = f_buffer[6];
      u = f_buffer[7];
      u_h = f_buffer[8];
      //
      //
      rho_1 = rho_2 = 1.0;
      //
      //*** Norm b
      double bNorm_new;
      double buff_fl = dom->Dot_Interior(b, b);
      MPI_Allreduce(&buff_fl, &bNorm_new, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Local2Border(b, border_buffer0);       //
      MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                                                  //  (rhs on border)
      buff_fl = bNorm_new + dot(border_buffer1, border_buffer1, size_b);
      bNorm_new = sqrt(buff_fl);
      // Check if the norm of b is samll enough for convengence
      if(CheckNormRHS(bNorm_new))
         return 0;
      //*** r = b-Ax
      //    A*x
      dom->Global2Local(xg, x, n);
      A->multiVec(x,v);                           // v as buffer
      for(i=0; i<size; i++)
         r0[i] = b[i]-v[i];                       // r = b-Ax
      //   Collect border r
      dom->Local2Border(r0, border_buffer0);      //  buffer
      MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Border2Local(border_buffer1, r0);
      //
      for(i=0; i<size; i++)
      {
         r[i] = r0[i];
         v[i] = 0.;
      }
      if ((error = Norm(r)/bNorm) < tol)
      {
         Message();
         return 0;
      }
      //
      for (iter = 1; iter <= max_iter; iter++)
      {
         rho_1 = dot(r0, r);
         if (fabs(rho_1) < DBL_MIN)               //  DBL_EPSILON
            break;
         //
         if (iter == 1)
         {
            for(i=0; i<size; i++)
               p[i] = u[i] = r[i];
         }
         else
         {
            beta = rho_1/rho_2;
            for(i=0; i<size; i++)
            {
               u[i] = r[i] + beta * q[i];
               p[i] = u[i] + beta * (q[i] + beta * p[i]);
            }
         }
         // Preconditioner
         Precond(p, p_h);
         // A M^{-1}p-->v
         MatrixMulitVec(p_h, v);
         //
         alpha = rho_1 / dot(r0, v);
         //
         for(i=0; i<size; i++)
         {
            q[i] = u[i] - alpha * v[i];
            q_h[i] = u[i] + q[i];
         }
         // Preconditioner
         Precond(q_h, u_h);
         for(i=0; i<size; i++)
            x[i] += alpha * u_h[i];
         //
         MatrixMulitVec(u_h, q_h);
         //
         for(i=0; i<size; i++)
            r[i] -= alpha * q_h[i];
         rho_2 = rho_1;
         if ((error = Norm(r) / bNorm) < tol)
            break;
      }
      // concancert internal x
      dom->CatInnerX(xg, x, n);
      //
      Message();
      //
      return iter <= max_iter;
   }
#define  TEST_MPII
   /*************************************************************************
   GeoSys-Function:
   Task: Parallel BiCGStab solver
       xg:  Global solution
       n:  Size of x
   Programming:
   12/2007 WW
   **************************************************************************/
   int Linear_EQS::BiCGStab(double *xg, const long n)
   {
      long i, size, size_b;                       //, size_t;
      double rho_0, rho_1, alpha, beta, omega, tt=0., norm_v=0.;
      double *r0, *r, *s, *s_h, *t, *v, *p, *p_h;
      double *r0_b, *r_b, *s_b, *s_h_b, *t_b, *v_b, *p_b, *p_h_b;
      double *x_b;
      //
      size = A->Dim();
      //
      size_b = dom->BSize()*A->Dof();
      // size_t = size + size_b; //[0, size_i): internal; [size_i, size_i+size_b): border.
      r0 = f_buffer[0];
      r = f_buffer[1];
      s = f_buffer[2];
      s_h = f_buffer[3];
      t = f_buffer[4];
      v = f_buffer[5];
      p = f_buffer[6];
      p_h = f_buffer[7];
      //
      //
      rho_0 = alpha = omega = 1.0;

#ifdef  TEST_MPI
      //TEST
      string test = "rank";
      char stro[64];
      sprintf(stro, "%d",myrank);
      string test1 = test+(string)stro+"Assemble.txt";
      ofstream Dum(test1.c_str(), ios::out);
      Dum.width(20);
      Dum.precision(15);
      Dum.setf(ios::scientific);

      Dum<<"Time step: "<<aktueller_zeitschritt<<endl;
      Dum<<"Norm b inner  " <<dom->Dot_Interior(b, b)<<endl;
      //  if(A->Dof()==1)
      //  Dum.close();
#endif

      //*** Norm b
      double bNorm_new;
      double buff_fl = dom->Dot_Interior(b, b);
      MPI_Allreduce(&buff_fl, &bNorm_new, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Local2Border(b, border_buffer0);       // buffer
      MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                                                  // (rhs on border)
      buff_fl = bNorm_new + dot(border_buffer1, border_buffer1, size_b);
      bNorm_new = sqrt(buff_fl);

      // Check if the norm of b is samll enough for convengence
      if(CheckNormRHS(bNorm_new))
         return 0;
      //*** r = b-Ax
      //    A*x
      dom->Global2Local(xg, x, n);
      A->multiVec(x,s);                           // s as buffer
      for(i=0; i<size; i++)
         r0[i] = b[i]-s[i];                       // r = b-Ax
      //   Collect border r
      dom->Local2Border(r0, border_buffer0);      // buffer
      MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Border2Local(border_buffer1, r0);

#ifdef  TEST_MPI
      //TEST
      //if(A->Dof()>1)
      {
         Dum<<" |r_0|= "<<Norm(r0)<<endl;
         Dum<<" |b_0|= "<< bNorm<<endl;
         Dum<<" x  "<<endl;
         for(i=0; i<n; i++)
            Dum<<xg[i]<<endl;
         /*
         Dum<<" b  "<<endl;
         for(i=0; i<size; i++)
           Dum<<b[i]<<endl;
         Dum<<"inter: r  Ax  "<<endl;
         for(i=0; i<size; i++)
          Dum<<r0[i]<<endl;
         Dum<<"border: r  Ax  "<<endl;
         for(i=0; i<size_b; i++)
          Dum<<r0_b[i]<<endl;

         Dum<<"inter: x  "<<endl;
         for(i=0; i<size; i++)
         Dum<<x[i]<<endl;

         Dum<<"border: x  "<<endl;
         for(i=0; i<size_b; i++)
         Dum<<x_b[i]<<endl;
         */
      }
#endif

      // Initial. [0, size_i): internal; [size_i, size_i+size_b): border.
      for(i=0; i<size; i++)
      {
         r[i] = r0[i];
         v[i] = 0.;
         p[i] = 0.;
      }
      if ((error = Norm(r)/bNorm) < tol)
      {
         if(myrank==0)                            // Make screen output only by processor 0
            Message();
         return 0;
      }

      /*
        //TEST
      //if(A->Dof()>1)
        {
          Dum<<" bNorm  "<<bNorm <<"  Norm(r) "<<  Norm(r) <<endl;

        }
      */

      //
      for (iter = 1; iter <= max_iter; iter++)
      {
         rho_1 = dot(r0, r);
         if (fabs(rho_1) < DBL_MIN)
            break;

#ifdef  TEST_MPI
         //TEST
         //  if(A->Dof()>1)
         Dum<<" rho_1  "<<rho_1 <<endl;
#endif

         if (iter == 1)
         {
            // p[0, size_i): internal; p[size_i, size_i+size_b)-->p_b: border.
            for(i=0; i<size; i++)
               p[i] = r[i];
         }
         else
         {
            beta = (rho_1/rho_0) * (alpha/omega);
            // [0, size_i): internal; [size_i, size_i+size_b): border.
            for(i=0; i<size; i++)
               p[i] = r[i] + beta * (p[i] - omega * v[i]);
         }
         // Preconditioner
         Precond(p, p_h);
         // A M^{-1}p-->v
         MatrixMulitVec(p_h, v);
         //
         alpha = rho_1 / dot(r0, v);

#ifdef  TEST_MPI
         //TEST
         //  if(A->Dof()>1)
         {
            //TEST
            Dum<<"  alpha  "<<alpha <<" dot(r0, v) "<<dot(r0, v)<<" dot(r0, r) "<<dot(r0, r)  <<endl;

            Dum<<"\n r0, r,  v   "<<endl;
            //  for(i=0; i<size_t; i++)
            //   Dum<<r0[i]<<"       "<<r[i]<<"       "<<v[i]<<endl;
         }
#endif

         //
         for(i=0; i<size; i++)
            s[i] = r[i] - alpha * v[i];
         norm_v = sqrt(dot(s,s));
         if ((error = norm_v/bNorm) < tol)
         {
            for(i=0; i<size; i++)
               x[i] += alpha * p_h[i];
            break;
         }

#ifdef  TEST_MPI
         /*
         //TEST
         if(A->Dof()>1)
           {
         //TEST
         Dum<<"\n  norm_v/bNorm  "<<error <<endl;
         for(i=0; i<size_t; i++)
          Dum<<s[i]<<endl;
         exit(0);
           }
         */
#endif

         //  M^{-1}s,
         Precond(s, s_h);
         // A* M^{-1}s
         MatrixMulitVec(s_h, t);
         //
         tt = dot(t,t);

#ifdef  TEST_MPI
         //TEST
         //TEST
         //  if(A->Dof()>1)
         Dum<<"  tt  "<<tt <<endl;
#endif

         if(tt>DBL_MIN)
            omega = dot(t,s) / tt;
         else
            omega = 1.0;
         // Update solution
         for(i=0; i<size; i++)
         {
            x[i] += alpha * p_h[i] + omega * s_h[i];
            r[i] = s[i] - omega * t[i];
         }
         rho_0 = rho_1;
         //
         norm_v = sqrt(dot(r,r));

#ifdef  TEST_MPI
         //TEST
         //TEST
         // if(A->Dof()>1)
         {
            Dum<<" sqrt(dot(r,r))  "<<norm_v  <<endl;
            // exit(0);
         }
#endif

         if ((error = norm_v / bNorm) < tol)
            break;
         if (fabs(omega) < DBL_MIN)
         {
            error = norm_v / bNorm;
            break;
         }

#ifdef  TEST_MPI
         //TEST
         // Dum.close();
         //MPI_Finalize();
         // exit(0);
#endif

      }
      //
      // concancert internal x
      dom->CatInnerX(xg, x, n);
      // Form the local to global
      // dom->Border2Global(x_b, xg, n);

#ifdef  TEST_MPI
      // if(A->Dof()>1)
      {
         Dum<<" x "<<endl;
         for(i=0; i<n; i++)
            Dum<<xg[i]<<endl;
         Dum.close();
         // exit(0);
      }
#endif

      //
      Message();
      //
      //  return iter <= max_iter;
      return iter ;
   }

   //#define  CG_test
   /*************************************************************************
   GeoSys-Function:
   Task: Parallel BiCG solver
       xg:  Global solution
       n:  Size of x
   Programming:
   08/2008 WW
   **************************************************************************/
   int Linear_EQS::BiCG(double *xg, const long n)
   {
      long i, size, size_b;
      double rho1, rho2, alpha, beta;
      double *z, *zt, *p, *pt, *q, *qt, *r, *rt;
      //
      size = A->Dim();
      //
      size_b = dom->BSize()*A->Dof();
      z = f_buffer[0];
      zt = f_buffer[1];
      p = f_buffer[2];
      pt = f_buffer[3];
      q = f_buffer[4];
      qt = f_buffer[5];
      r = f_buffer[6];
      rt = f_buffer[7];
      //
      //*** Norm b
      double bNorm_new;
      double buff_fl = dom->Dot_Interior(b, b);
      MPI_Allreduce(&buff_fl, &bNorm_new, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Local2Border(b, border_buffer0);       //
      MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                                                  //  (rhs on border)
      buff_fl = bNorm_new + dot(border_buffer1, border_buffer1, size_b);
      bNorm_new = sqrt(buff_fl);
      // Check if the norm of b is samll enough for convengence

      // Check if the norm of b is samll enough for convengence
      if(CheckNormRHS(bNorm_new))
         return 0;
      //*** r = b-Ax
      //    A*x
      dom->Global2Local(xg, x, n);
      A->multiVec(x,rt);                          // rt as buffer
      for(i=0; i<size; i++)
         r[i] = b[i]-rt[i];                       // r = b-Ax
      //   Collect border r
      dom->Local2Border(r, border_buffer0);       //  buffer
      MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      dom->Border2Local(border_buffer1, r);
      //
      // Initial.
      for(i=0; i<size; i++)
      {
         r[i] = b[i]-rt[i];
         rt[i] = r[i];
      }
      if ((error = Norm(r)/bNorm) < tol)
      {
         if(myrank==0)                            // Make screen output only by processor 0
            Message();
         return 0;
      }

#ifdef CG_test
      //TEST
      string test = "rank";
      char stro[64];
      sprintf(stro, "%d",myrank);
      string test1 = test+(string)stro+"_Assemble.txt";
      ofstream Dum(test1.c_str(), ios::out);
      Dum.width(20);
      Dum.precision(15);
      Dum.setf(ios::scientific);

      Dum<<" Norm(r) "<< Norm(r)<<" bNorm "<<bNorm<<endl;
#endif

      //
      for (iter = 1; iter <= max_iter; iter++)
      {
         Precond(r, z);
         TransPrecond(rt, zt);
         rho1 = dot(z, rt);

#ifdef CG_test
         Dum<<" rho1 "<< rho1<<endl;
#endif

         //
         if (fabs(rho1) < DBL_MIN)
         {
            Message();
            break;
         }
         //
         if (iter == 1)
         {
            for(i=0; i<size; i++)
            {
               p[i] = z[i];
               pt[i] = zt[i];
            }
         }
         else
         {
            beta = rho1/rho2;
            for(i=0; i<size; i++)
            {
               p[i] = z[i] + beta*p[i];
               pt[i] = zt[i] + beta*pt[i];
            }
         }
         MatrixMulitVec(p, q);
         TransMatrixMulitVec(pt, qt);
         alpha = rho1/dot(pt,q);

#ifdef CG_test
         Dum<<" alpha "<< alpha<<endl;
#endif

         for(i=0; i<size; i++)
         {
            x[i] += alpha*p[i];
            r[i] -= alpha*q[i];
            rt[i] -= alpha*qt[i];
         }
         //
         rho2 = rho1;
         if ((error = Norm(r)/bNorm) < tol)
            break;

#ifdef CG_test
         Dum<<" error = Norm(r)/bNorm "<< error<<endl;
#endif

         //
      }
      //
      // concancert internal x
      dom->CatInnerX(xg, x, n);

      //
      Message();
#ifdef CG_test
      exit(0);
#endif

      //
      return iter <= max_iter;
   }
#endif
   //------------------------------------------------------------------------
}                                                 // namespace
#endif                                            // if defined(NEW_EQS)
