/*========================================================================
 GeoSys - class Matrix, Sparse matrix (Declaration)
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation.
 Function:   See the declaration below
 Design and programm WW
 03/2010 some improvements TF
==========================================================================*/
#ifndef matrix_class_INC
#define matrix_class_INC

#include "Configure.h"

#include <iostream>
#include <fstream>

#ifdef NEW_EQS
namespace Mesh_Group {class CFEMesh;}
//08.2007 WW
class CPARDomain;
using Mesh_Group::CFEMesh;
#endif
//#define OverLoadNEW_DELETE

namespace Math_Group
{

   class Matrix
   {
      public:
         Matrix(size_t rows, size_t cols=1);
         Matrix();
         explicit Matrix(const Matrix& m);
         //
         void resize(size_t rows, size_t cols=1);
         //
         virtual ~Matrix();

         // Operators
         virtual void operator= (double a);
         virtual void operator*= (double a);
         virtual void operator/= (double a);
         virtual void operator+= (double a);
         void operator= (const Matrix& m);
         void operator+= (const Matrix& m);
         void operator-= (const Matrix& m);

         void GetTranspose(Matrix& m);

         // vec_result = This*vec. vec_result must be initialized
         void multi(const double *vec, double *vec_result, double fac=1.0);
         // m_result = this*m. m_result must be initialized
         void multi(const Matrix& m, Matrix& m_result, double fac=1.0);
         // m_result = this*m1*m2. m_result must be initialized
         void multi(const Matrix& m1, const Matrix& m2, Matrix& m_result);

         // Access to members
         virtual double& operator() (size_t i, size_t j=0) const;
         void LimitSize(size_t nRows, size_t nCols=1);

         size_t Rows() const {return nrows;}
         size_t Cols() const {return ncols;}
         size_t Size() const {return size;}

         // Print
         void Write(std::ostream& os = std::cout);
         void Write_BIN(std::fstream& os);
         void Read_BIN(std::fstream& is);
      protected:
         size_t nrows, nrows0;
         size_t ncols, ncols0;
         size_t size;
         double *data;
         bool Sym;
   };

   // Symmetrical matrix. 12-01-2005. WW
   class SymMatrix : public Matrix
   {
      public:
         SymMatrix(size_t dim);
         SymMatrix();
         explicit SymMatrix(const SymMatrix& m);

         void resize(size_t dim);
         ~SymMatrix() {};

         // Operators
         void operator= (double a);
         void operator*= (double a);
         void operator+= (double a);
         void operator= (const SymMatrix& m);
         void operator+= (const SymMatrix& m);
         void operator-= (const SymMatrix& m);
         void LimitSize(size_t dim);

         // Access to members
         double& operator() (size_t i, size_t j) const;
   };

   typedef Matrix Vec;

   /*========================================================================
    GeoSys - class my_vector (Declaration)
    Task:       Carry out vector operation
    Function:   See the declaration below
    programming:
    05/2005  WW
   ==========================================================================*/
   template<class T> class vec
   {
      public:
         vec(int argSize);
         vec() : size(0) {}
         explicit vec(const vec<T>& v);

         virtual ~vec();
         // Operator
         virtual void operator = (T v) { for (size_t i=0; i<size; i++) entry[i] = v; }
         virtual void operator = (const vec<T>&);
         virtual void resize(int newh);
         virtual T& operator[] (int i) { return (T&) entry[i]; }
         virtual const T& operator[] (int i) const {return (const T&) entry[i]; }
         virtual size_t Size() const { return size; }

         T* Entry()        { return entry; }
         T* Entry()  const { return entry; }

         virtual void Write(std::ostream& os = std::cout) const;

      protected:
         T* entry;
         size_t size;

   };

   template<> class vec<void*>
   {
      public:
         vec(int argSize);
         vec() : size (0) {}
         explicit vec(const vec<void*>& v);

         virtual ~vec();
         // Operator
         void operator = (void* v)
            { for (size_t i=0; i<size; i++) entry[i] = v; }
         void operator = (const vec<void*>& v);
         void*& operator[] (int i) { return entry[i]; }
         const void*& operator[] (int i) const { return (const void*&) entry[i]; }

         // Access to memebers
         void** Entry()          { return entry; }
         const void** Entry()  const { return (const void**)entry; }

         virtual void resize(int newh);
         virtual size_t Size() const { return size; }
         virtual void Write(std::ostream& os = std::cout) const;

      protected:
         void** entry;
         size_t size;
   };

   template<class T> class vec<T*> : public vec<void*>
   {
      public:
         vec( int Size) : vec<void*>(Size) { }
         vec()               : vec<void*>()     { }
         explicit vec(const vec<T*>& v) : vec<void*>(v) { }

         ~vec() { }

         // Operator
         void operator = (T* v) { for (size_t i=0; i<size;i++) entry[i] = v; }
         void operator = (const vec<T*>& v);
         T*& operator[] (int i) { return (T*&) entry[i]; }
         const T*& operator[] (int i) const {return (const T*&) entry[i]; }

         T** Entry()        { return entry; }
         T** Entry()  const { return (const T**)entry; }

   };

#ifdef NEW_EQS
   //
   class SparseTable
   {
      public:
         SparseTable(CFEMesh *a_mesh, bool quadratic, bool symm=false);
         SparseTable(CPARDomain &m_dom, bool quadratic, bool symm=false);
         ~SparseTable();
         void Write(std::ostream &os=std::cout);
      private:
         bool symmetry;
         // Topology mapping from data array to matrix
         long *entry_column;
         long *num_column_entries;                // number of entries of each columns in sparse table
         long *row_index_mapping_n2o;             // Row index of sparse table to row index of matrix
         long *row_index_mapping_o2n;             // Inverse of last
         long *diag_entry;                        // Global index to the index of  entry_column
         long size_entry_column;
         long max_columns;
         long rows;
         friend class CSparseMatrix;
   };
   //08.2007 WW
   //
   class CSparseMatrix
   {
      public:
         CSparseMatrix(const SparseTable &sparse_table, const int dof);
         ~CSparseMatrix();
         // Preconditioner
         void Precond_Jacobi(double *vec_s, double *vec_r);
                                                  //TEMP
         void Precond_ILU(double *vec_s, double *vec_r)
         {
            vec_s = vec_r = NULL;
         }
         // Operator
         void operator = (const double a);
         void operator *= (const double a);
         void operator += (const double a);
         void operator = (const CSparseMatrix& m);
         void operator += (const CSparseMatrix& m);
         void operator -= (const CSparseMatrix& m);
         // Vector pass through augment and bring results back.
         void multiVec(double *vec_s, double *vec_r);
         void Trans_MultiVec(double *vec_s, double *vec_r);
         void Diagonize(const long idiag, const double b_given, double *b);
         //
         // Access to members
         double& operator() (const long i, const long j=0) const;
         //
         long Dim() const {return DOF*rows;}
         int Dof() const {return DOF;}
         void SetDOF(const int dof_n)             //_new. 02/2010. WW
         {
            DOF = dof_n;
         }
         long Size() const {return rows;}
#ifdef LIS                                  // These two pointers are in need for Compressed Row Storage
         int nnz() const                          // PCH
         {
            return DOF*DOF*size_entry_column;
         }
         int *ptr;
         int *col_idx;
         int* entry_index;
         int GetCRSValue(double* value);
#endif
         // Print
         void Write(std::ostream &os=std::cout);

         // Domain decomposition
#if defined(USE_MPI)
         void DiagonalEntries(double *diag_e);
#endif
      private:
         // Data
         double *entry;
         mutable double zero_e;
         //
         bool symmetry;
         // Topology mapping from data array to matrix. All are only pointers to the
         // correpinding members in SparseTable, and no memory are allocated for them
         long *entry_column;
         long *num_column_entries;                // number of entries of each columns in sparse table
         long *row_index_mapping_n2o;             // Row index of sparse table to row index of matrix
         long *row_index_mapping_o2n;             // Inverse of last
         long *diag_entry;
         long size_entry_column;
         long max_columns;
         long rows;
         //
         int DOF;
   };
   // Since the pointer to member funtions gives lower performance
#endif

   //
   // Cross production x^y. WW 12.01.2005
   //const Vec& operator ^ (Vec& x,  Vec& y);

   // End of class Matrix
}


//==========================================================================
#endif
