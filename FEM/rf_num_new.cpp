/**************************************************************************
FEMLib - Object: NUM
Task:
Programing:
11/2004 OK Implementation
last modified:
**************************************************************************/
// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#ifdef USE_MPI                                    //WW
#include "mpi.h"
#include "par_ddc.h"
//#undef SEEK_SET
//#undef SEEK_END
//#undef SEEK_CUR
#endif

#include "makros.h"
// C++ STL
#include <cfloat>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
// FEM-Makros
#include "makros.h"
#include "files0.h"
extern ios::pos_type GetNextSubKeyword(ifstream* file,string* line, bool* keyword);
// GeoSys-GeoLib
// GeoSys-FEMLib
#include "rf_num_new.h"
#ifndef NEW_EQS                                   //WW. 06.11.2008
#include "matrix.h"
#endif
#include "mathlib.h"
#include "rf_pcs.h"
// GeoSys-MSHLib

extern size_t max_dim;                            //OK411 todo

//==========================================================================
vector<CNumerics*>num_vector;
/**************************************************************************
FEMLib-Method:
Task: constructor
Programing:
11/2004 OK Implementation
10/2005 OK pcs_type_name
07/2007 OK DEFORMATION
**************************************************************************/
CNumerics::CNumerics(string name)
{
   pcs_type_name = name;                          //OK
   // GLOBAL
   renumber_method = 0;
   // LS - Linear Solver
   ls_method = 2;                                 //OK41
   ls_max_iterations = 1000;
   ls_error_method = 1;
   ls_error_tolerance = 1e-12;
   ls_theta = 1.0;
   ls_precond = 1;
   ls_storage_method = 2;                         //OK41
   // NLS - Nonlinear Solver
   nls_method_name = "PICARD";
   nls_method = 0;                                // Default, Picard. 1: Newton
   nls_max_iterations = 1;                        //OK
   nls_error_tolerance = 1.0e-4;
   nls_error_tolerance_local = 1.0e-10;           //For element level
   nls_relaxation = 0.0;
   // cpl WW
   cpl_iterations = 1;                            //OK
   cpl_tolerance = 1.0e-3;
   cpl_variable = "FLUX";
   // ELE
   ele_gauss_points = 3;
   ele_mass_lumping = 0;
   ele_upwind_method = 0;                         //CB
   ele_upwinding = 0;
   ele_supg_method = 0;                           //NW
   ele_supg_method_length = 0;                    //NW
   ele_supg_method_diffusivity = 0;               //NW
   //----------------------------------------------------------------------
   // Deformation
   GravityProfile = 0;
   DynamicDamping = NULL;                         //WW
   if(pcs_type_name.compare("DEFORMATION")==0)
   {
      ls_method = 2;
      ls_error_method = 2;
      ls_error_tolerance = 1e-12;
      ls_max_iterations = 2000;
      ls_precond = 100;
      ls_storage_method = 4;
   }
   //----------------------------------------------------------------------
   if(pcs_type_name.compare("RICHARDS_FLOW")==0)
   {
      ele_mass_lumping = 1;
      ele_upwinding = 0.5;
      ls_max_iterations = 2000;
      ls_error_method = 2;
      ls_error_tolerance = 1e-10;
      ls_precond = 4;
      ls_storage_method = 4;
      nls_max_iterations = 25;
      nls_error_tolerance = 1.0e-3;
   }
   //
}


/**************************************************************************
FEMLib-Method:
Task: deconstructor
Programing:
11/2004 OK Implementation
**************************************************************************/
CNumerics::~CNumerics(void)
{
   if(DynamicDamping) delete DynamicDamping;
   DynamicDamping = NULL;
}


/**************************************************************************
FEMLib-Method:
Task: OBJ read function
Programing:
11/2004 OK Implementation
**************************************************************************/
bool NUMRead(string file_base_name)
{
   //----------------------------------------------------------------------
   //OK  NUMDelete();
   //----------------------------------------------------------------------
   CNumerics *m_num = NULL;
   char line[MAX_ZEILE];
   string sub_line;
   string line_string;
   ios::pos_type position;
   //========================================================================
   // File handling
   string num_file_name = file_base_name + NUM_FILE_EXTENSION;
   ifstream num_file (num_file_name.data(),ios::in);
   if (!num_file.good())
      return false;
   num_file.seekg(0L,ios::beg);
   //========================================================================
   // Keyword loop
   cout << "NUMRead" << endl;
   while (!num_file.eof())
   {
      num_file.getline(line,MAX_ZEILE);
      line_string = line;
      if(line_string.find("#STOP")!=string::npos)
         return true;
      //----------------------------------------------------------------------
                                                  // keyword found
      if(line_string.find("#NUMERICS")!=string::npos)
      {
         m_num = new CNumerics("default");
         position = m_num->Read(&num_file);
         num_vector.push_back(m_num);
         num_file.seekg(position,ios::beg);
      }                                           // keyword found
   }                                              // eof
   return true;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 WW Implementation
**************************************************************************/
bool CNumerics::CheckDynamic()
{
   if(DynamicDamping)
      return true;
   else
      return false;
}


/**************************************************************************
FEMLib-Method:
Task: OBJ read function
Programing:
11/2004 OK Implementation
**************************************************************************/
ios::pos_type CNumerics::Read(ifstream *num_file)
{
   string line_string;
   bool new_keyword = false;
   bool new_subkeyword = false;
   ios::pos_type position;
   ios::pos_type position_subkeyword;
   std::stringstream line;
   //========================================================================
   // Schleife ueber alle Phasen bzw. Komponenten
   while(!new_keyword)
   {
      if(new_subkeyword)
         num_file->seekg(position,ios::beg);
      new_subkeyword = false;
      position = GetNextSubKeyword(num_file,&line_string,&new_keyword);
      if(new_keyword)
         return position;
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$PCS_TYPE")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));
         line >> pcs_type_name;
         line.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$RENUMBER")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));
         line >> renumber_method;
         if(renumber_method==2)
            line >> renumber_parameter;
         line.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$NON_LINEAR_SOLVER")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));
         line >> nls_method_name;
         line >> nls_error_tolerance;
         if(nls_method_name.find("NEWTON")!=string::npos)
            line>>nls_error_tolerance_local;
         nls_method = 1;
         line >> nls_max_iterations;
         line >> nls_relaxation;
         line.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$LINEAR_SOLVER")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));
         line >> ls_method;
         line >> ls_error_method >> ls_error_tolerance;
         line >> ls_max_iterations;
         line >> ls_theta;
         line >> ls_precond;
         line >> ls_storage_method;
         line.clear();
         continue;
      }
      //....................................................................
                                                  // subkeyword found
      if(line_string.find("$ELE_GAUSS_POINTS")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));
         line >> ele_gauss_points;                // probably element-type-wise
         line.clear();
         continue;
      }
                                                  // subkeyword found
      if(line_string.find("$ELE_MASS_LUMPING")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));
         line >> ele_mass_lumping;
         line.clear();
         continue;
      }
                                                  // subkeyword found
      if(line_string.find("$ELE_UPWINDING")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));
                                                  //CB now read also upwinding method
         line >> ele_upwinding >> ele_upwind_method  ;
         line.clear();
         continue;
      }
                                                  // subkeyword found
      if(line_string.find("$ELE_SUPG")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));
                                                  //NW
         line >> ele_supg_method >> ele_supg_method_length >> ele_supg_method_diffusivity;
         line.clear();
         cout << "->SUPG method is selected." << endl;
         continue;
      }
                                                  // subkeyword found
      if(line_string.find("$GRAVITY_PROFILE")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));    //WW
         line >> GravityProfile;
         line.clear();
         continue;
      }
                                                  // subkeyword found
      if(line_string.find("$DYNAMIC_DAMPING")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));    //WW
         DynamicDamping = new double[3];
         // Default
         DynamicDamping[0] = 0.515;
         DynamicDamping[1] = 0.51;
         DynamicDamping[2] = 0.51;
         line >> DynamicDamping[0] >> DynamicDamping[1] >> DynamicDamping[2];
         line.clear();
         continue;
      }
                                                  // WW subkeyword found
      if(line_string.find("$COUPLING_ITERATIONS")!=string::npos)
      {
         line.str(GetLineFromFile1(num_file));
         line >> cpl_variable;                    //pcs_name. WW MB
         line >> cpl_iterations;
         line >> cpl_tolerance;
         line.clear();
         continue;
      }
      //....................................................................
      /*
          if(line_string.find("$TIME_STEPS")!=string::npos) { // subkeyword found
            while((!new_keyword)||(!new_subkeyword)||(!num_file->eof())){
              position = num_file->tellg();
              line_string = GetLineFromFile1(num_file);
              if(line_string.find("#")!=string::npos){
                return position;
              }
              if(line_string.find("$")!=string::npos){
                new_subkeyword = true;
                break;
      }
      line.str(line_string);
      line >> no_time_steps;
      line >> time_step_length;
      for(i=0;i<no_time_steps;i++)
      time_step_vector.push_back(time_step_length);
      line.clear();
      }
      }
      */
      //....................................................................
   }
   return position;
}


/**************************************************************************
FEMLib-Method:
Task: master write function
Programing:
11/2004 OK Implementation
last modification:
**************************************************************************/
void NUMWrite(string base_file_name)
{
   CNumerics *m_num = NULL;
   string sub_line;
   string line_string;
   //========================================================================
   // File handling
   string num_file_name = base_file_name + NUM_FILE_EXTENSION;
   fstream num_file (num_file_name.data(),ios::trunc|ios::out);
   num_file.setf(ios::scientific,ios::floatfield);
   num_file.precision(12);
   if (!num_file.good()) return;
   num_file.seekg(0L,ios::beg);
   //========================================================================
   num_file << "GeoSys-NUM: Numerics ------------------------------------------------\n";
   //========================================================================
   // OUT vector
   int num_vector_size =(int)num_vector.size();
   int i;
   for(i=0;i<num_vector_size;i++)
   {
      m_num = num_vector[i];
      m_num->Write(&num_file);
   }
   num_file << "#STOP";
   num_file.close();
}


/**************************************************************************
FEMLib-Method:
Task: write function
Programing:
11/2004 OK Implementation
last modification:
**************************************************************************/
void CNumerics::Write(fstream* num_file)
{
   //KEYWORD
   *num_file  << "#NUMERICS" << endl;
   //--------------------------------------------------------------------
   /*OK
    *num_file << " $METHOD" << endl;
    *num_file << method_name << endl;
     if(method_name.find("LAGRANGE")!=string::npos){
       *num_file << lag_quality << " " << lag_max_steps << " " << lag_local_eps << " ";
       *num_file << lag_time_weighting << " " << lag_min_weight << " ";
       *num_file << lag_use_matrix << " " << lag_vel_method;
       *num_file << endl;
     }
   */
   //--------------------------------------------------------------------
   *num_file << " $PCS_TYPE" << endl;
   *num_file << "  " << pcs_type_name << endl;
   //--------------------------------------------------------------------
   *num_file << " $NON_LINEAR_SOLVER" << endl;
   *num_file << "  " << nls_method_name;
   *num_file << " "  << nls_error_tolerance;
   *num_file << " "  << nls_max_iterations;
   *num_file << " "  << nls_relaxation;
   *num_file << endl;
   //--------------------------------------------------------------------
   *num_file << " $LINEAR_SOLVER" << endl;
   *num_file << "  " << ls_method;
   *num_file << " "  << ls_error_method;
   *num_file << " "  << ls_error_tolerance;
   *num_file << " "  << ls_max_iterations;
   *num_file << " "  << ls_theta;
   *num_file << " "  << ls_precond;
   *num_file << " "  << ls_storage_method;
   *num_file << endl;
   //--------------------------------------------------------------------
   *num_file << " $ELE_GAUSS_POINTS" << endl;
   *num_file << "  " << ele_gauss_points;
   *num_file << endl;
   //--------------------------------------------------------------------
   *num_file << " $ELE_MASS_LUMPING" << endl;
   *num_file << "  " << ele_mass_lumping;
   *num_file << endl;
   //--------------------------------------------------------------------
   *num_file << " $ELE_UPWINDING" << endl;
   *num_file << "  " << ele_upwinding;
   *num_file << endl;
   //--------------------------------------------------------------------
}


//////////////////////////////////////////////////////////////////////////
// LINEAR_SOLVER
//////////////////////////////////////////////////////////////////////////
#ifndef NEW_EQS                                   //WW. 06.11.2008
/**************************************************************************
 ROCKFLOW - Funktion: NormOfUnkonwn

 Aufgabe:
   Compute the norm of RHS of a linear equation

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

 Ergebnis:
   - double - Eucleadian Norm

Programmaenderungen:
12/2002   WW   Erste Version

**************************************************************************/
double CalcNormOfRHS(LINEAR_SOLVER*ls)
{
   int i, j;
   int unknown_vector_dimension;
   long number_of_nodes;
   double NormW = 0.0;

   if (!ls) {printf(" \n Warning: solver not defined, exit from loop_ww.cc"); exit(1);}
   /* Ergebnisse eintragen */
   unknown_vector_dimension = GetUnknownVectorDimensionLinearSolver(ls);
   for (i = 0; i < unknown_vector_dimension; i++)
   {
      number_of_nodes=ls->unknown_node_numbers[i];
      for(j=0; j<number_of_nodes; j++)
         NormW += ls->b[number_of_nodes*i+j]*ls->b[number_of_nodes*i+j];

   }
   return sqrt(NormW);
}


/*************************************************************************
 ROCKFLOW - Funktion: SetZeroLinearSolver

 Aufgabe:
   Setzt Matrix, rechte Seite und oesungsvektor auf NULL.

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

 Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

Programmaenderungen:
02/1999     AH         Erste Version

*************************************************************************/
LINEAR_SOLVER *SetZeroLinearSolver(LINEAR_SOLVER*ls)
{
   /* Initialisieren der Felder */
   if (!ls)
      return NULL;
   SetLinearSolver(ls);
   if (ls->matrix)
   {
      MXInitMatrix();
   }
   if (ls->b)
      MNulleVec(ls->b, ls->dim);
   if (ls->x)
      MNulleVec(ls->x, ls->dim);
   /*if (ls->xx) MNulleVec(ls->xx,ls->dim); */
   /*if (ls->memory) MNulleVec(ls->memory,ls->dim); */

   return ls;
}


/*************************************************************************
 ROCKFLOW - Funktion: SetLinearSolver

 Aufgabe:
   Konstruktor fuer LINEAR_SOLVER

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

 Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

Programmaenderungen:
02/1999     AH         Erste Version
7/2000     CT         Implizite Kenntnis ueber Matrixspeichermodell

*************************************************************************/
LINEAR_SOLVER *SetLinearSolver(LINEAR_SOLVER * ls)
{
   if (!ls)
      return NULL;
   if (ls->matrix)
   {
      MXSetMatrixPointer(ls->matrix);
   }

   return ls;
}


/*************************************************************************
 ROCKFLOW - Funktion: GetUnknownVectorDimensionLinearSolver

 Aufgabe:
   Liefert die Dimension der Unbekannten im Vektors des lin. Loesers.

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

 Ergebnis:
- s. o. -

Programmaenderungen:
08/2000     MK         Erste Version

*************************************************************************/
int GetUnknownVectorDimensionLinearSolver(LINEAR_SOLVER*ls)
{
   if (!ls)
      return -1;
   else
      return ls->unknown_vector_dimension;
}


/*************************************************************************
 ROCKFLOW - Funktion: DestroyLinearSolver

 Aufgabe:
   Destruktor fuer LINEAR_SOLVER

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

 Ergebnis:
   - void -

Programmaenderungen:
02/1999     AH         Erste Version
08/1999     AH         Erweiterung memory  (neq)
10/1999     AH         Systemzeit
07/2000     C.Thorenz  Bugfix fuer MXDestroy

*************************************************************************/
LINEAR_SOLVER *DestroyLinearSolver(LINEAR_SOLVER * ls)
{
   if (!ls)
      return NULL;
   SetLinearSolver(ls);

   if (ls->system_time_name)
      ls->system_time_name = (char *) Free(ls->system_time_name);
   if (ls->system_time_assemble_function_name)
      ls->system_time_assemble_function_name = (char *) Free(ls->system_time_assemble_function_name);

   //OK    if (ls->name)
   //OK        ls->name = (char *) Free(ls->name);
   if (ls->b)
      ls->b = (double *) Free(ls->b);
   if (ls->x)
      ls->x = (double *) Free(ls->x);
   /* //WW
   if (ls->bc_name)
       ls->bc_name = (char *) Free(ls->bc_name);
   if (ls->bc_name2)
       ls->bc_name2 = (char *) Free(ls->bc_name2);
   if (ls->bc_name3)
       ls->bc_name3 = (char *) Free(ls->bc_name3);
   if (ls->sousin_name1)
       ls->sousin_name1 = (char *) Free(ls->sousin_name1);
   if (ls->sousin_name2)
       ls->sousin_name2 = (char *) Free(ls->sousin_name2);
   */
   if (ls->lsp_name)
      ls->lsp_name = (char *) Free(ls->lsp_name);
   if (ls->r)
      ls->r = (double *) Free(ls->r);
   if (ls->xx)
      ls->xx = (double *) Free(ls->xx);
   if (ls->memory)
      ls->memory = (double *) Free(ls->memory);
   if (ls->new_memory)
      ls->new_memory = (double **) DestroyMemoryLinearSolver(ls);
   if (ls->matrix)
   {
      MXSetMatrixPointer(ls->matrix);
      ls->matrix = MXDestroyMatrix();
   }
   /* Multiple unknowns (Dim) / multiple solvers (MultSolvers) */
   if (ls->name_ls)
      ls->name_ls = (char *) Free(ls->name_ls);

   if (ls->name_group_ls)
      ls->name_group_ls = (char *) Free(ls->name_group_ls);
   if (ls)
      ls = (LINEAR_SOLVER *) Free(ls);
   return NULL;
}


/*************************************************************************
 ROCKFLOW - Funktion: DestroyMemoryLinearSolver

 Aufgabe:
   Extra-Speicher-Destruktor einer Instanz vom Typ LINEAR_SOLVER.

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

 Ergebnis:
   - Ungleich NULL im Erfolgsfall -

Programmaenderungen:
03/2000    AH      Erste Version

*************************************************************************/
LINEAR_SOLVER *DestroyMemoryLinearSolver(LINEAR_SOLVER * ls)
{
   int i;

   if (!ls)
      return NULL;
   SetLinearSolver(ls);

   for (i = 0; i < ls->memory_number; i++)
   {
      if (ls->new_memory[i])
         ls->new_memory[i] = (double *) Free(ls->new_memory[i]);
   }
   if (ls->new_memory)
      ls->new_memory = (double **) Free(ls->new_memory);

   return ls;
}


/*************************************************************************
 ROCKFLOW - Funktion: InitLinearSolver

 Aufgabe:
   Konstruktor fuer LINEAR_SOLVER

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

 Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

Programmaenderungen:
02/1999     AH         Erste Version
08/1999     AH         Erweiterung memory  (neq)
7/2000     CT         Implizite Kenntnis ueber Matrixspeichermodell
Fehlerbereinigung

*************************************************************************/
LINEAR_SOLVER *InitLinearSolver(LINEAR_SOLVER * ls)
{
   if (!ls)
      return NULL;

   //#ifdef USE_MPI //WW
   //    ls->matrix = NULL;
   //#else
   if (ls->matrix)
   {
      MXSetMatrixPointer(ls->matrix);
      ls->matrix = MXDestroyMatrix();
   }
   ls->matrix = MXSetupMatrix(ls->dim, ls->store, 0l);
   //#endif
   if (ls->dim == 0)
      return ls;

   if (ls->b)
      ls->b = (double *) Free(ls->b);
   ls->b = (double *) Malloc(ls->dim * sizeof(double));
   if (ls->x)
      ls->x = (double *) Free(ls->x);
   ls->x = (double *) Malloc(ls->dim * sizeof(double));
   if (ls->r)
      ls->r = (double *) Free(ls->r);
   ls->r = (double *) Malloc(ls->dim * sizeof(double));

   if (ls->xx)
      ls->xx = (double *) Free(ls->xx);
   ls->xx = (double *) Malloc(ls->dim * sizeof(double));
   if (ls->memory)
      ls->memory = (double *) Free(ls->memory);
   ls->memory = (double *) Malloc(ls->dim * sizeof(double));

   return ls;
}


#ifdef USE_MPI                                    //WW
/////////////////////////////////////////////////////////////////////
/**************************************************************************
FEMLib-Method:
Task:
      Based on LINEAR_SOLVER *InitLinearSolver(LINEAR_SOLVER * ls)
07/2006 WW
last modified:
**************************************************************************/

LINEAR_SOLVER *InitVectorLinearSolver(LINEAR_SOLVER * ls)
{
   if (!ls)
      return NULL;
   if (ls->dim == 0)
      return ls;

   if (ls->b)
      ls->b = (double *) Free(ls->b);
   ls->b = (double *) Malloc(ls->dim * sizeof(double));
   if (ls->x)
      ls->x = (double *) Free(ls->x);
   ls->x = (double *) Malloc(ls->dim * sizeof(double));
   if (ls->r)
      ls->r = (double *) Free(ls->r);
   ls->r = (double *) Malloc(ls->dim * sizeof(double));

   if (ls->xx)
      ls->xx = (double *) Free(ls->xx);
   ls->xx = (double *) Malloc(ls->dim * sizeof(double));
   if (ls->memory)
      ls->memory = (double *) Free(ls->memory);
   ls->memory = (double *) Malloc(ls->dim * sizeof(double));

   return ls;
}
#endif
/////////////////////////////////////////////////////////////////////

/*************************************************************************
 ROCKFLOW - Funktion: ConfigSolverProperties

 Aufgabe:
   Konfiguration von Loeser-Parametern

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)

 Ergebnis:
   - void -

Programmaenderungen:
12/1999   OK   Implementierung

*************************************************************************/
void ConfigSolverProperties(void)
{
   /* Speichertechnik automatisch optimieren */
   switch (max_dim)
   {
      case 0:
         sp2_start = 3;
         sp2_inc = 1;
         break;
      case 1:
         sp2_start = 11;
         sp2_inc = 2;
         break;
      case 2:
         sp2_start = 33;
         sp2_inc = 6;
         break;
   }
   /* Umnummerierer-Methode auswaehlen */
   //WW ConfigRenumberProperties();
   /* umnummerierer muss definiert sein ! */
}


/*************************************************************************
 ROCKFLOW - Funktion: InitializeLinearSolver

 Aufgabe:
   Loesungsvektor vorbelegen.

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

 Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

Programmaenderungen:
02/1999     AH         Erste Version
08/1999     AH         Erweiterung memory  (neq)
09/2004      Remove BC/Sink incoorperation
*************************************************************************/
                                                  //WW
void SetLinearSolverType(LINEAR_SOLVER* ls ,CNumerics *m_num)
{
   if (!ls)
      return;
   /* Gruppeneingenschaften setzen */
   /*
     SetPropertiesLinearSolver(ls, prop_name);
     // Benutzer-Eigenschaften
     lsp = GetTLinearSolverProperties(ls->lsp_name, lsp, 0);
     if (lsp) {
       ls->lsp = lsp;
       ls->store = get_lsp_store(ls->lsp);
     }
     // Entwickler-Eigenschaften
     if (ls->lsp_name && !ls->lsp) {
       sprintf(string, "DEFAULT_");
   strcat(string, ls->lsp_name);
   lsp = GetLinearSolverProperties(string);
   if (lsp) {
   ls->lsp = lsp;
   ls->store = get_lsp_store(ls->lsp);
   }
   }
   // Master Default-Eigenschaften
   if (!ls->lsp) {
   ls->lsp = GetLinearSolverProperties(MASTER_DEFAULT_LINEAR_SOLVER_PROPERTIES);
   if (ls->lsp) {
   ls->store = get_lsp_store(ls->lsp);
   }
   }
   */
   ls->store = m_num->ls_storage_method;
   switch(m_num->ls_method)                       //OK (ls->lsp->type){
   {
      case 1:
         ls->LinearSolver = SpGauss;
         break;
      case 2:
#ifndef USE_MPI
         ls->LinearSolver = SpBICGSTAB;
#endif
         break;
      case 3:
         ls->LinearSolver = SpBICG;
         break;
      case 4:
         ls->LinearSolver = SpQMRCGSTAB;
         break;
      case 5:
         ls->LinearSolver = SpCG;
         break;
      case 6:
         ls->LinearSolver = SpCGNR;
         break;
      case 7:
         ls->LinearSolver = SpCGS;
         break;
      case 8:
         ls->LinearSolver = SpRichardson;
         break;
      case 9:
         ls->LinearSolver = SpJOR;
         break;
      case 10:
         ls->LinearSolver = SpSOR;
         break;
      case 11:
         ls->LinearSolver = SpAMG1R5;
         break;
      case 12:
         ls->LinearSolver = SpUMF;
         break;
   }

}


//WW
LINEAR_SOLVER *InitializeLinearSolver(LINEAR_SOLVER* ls ,CNumerics *m_num)
{
   SetLinearSolverType(ls,m_num);
#ifdef USE_MPI                                 //WW
   InitVectorLinearSolver(ls);
#else
   InitLinearSolver(ls);
#endif
   /* Internen Speicher allokieren */
   InitMemoryLinearSolver(ls,ls->memory_number);
   /* Speicher initialisieren */
   SetMemoryZeroLinearSolver(ls);
   return ls;
}


/*************************************************************************
 ROCKFLOW - Funktion: InitMemoryLinearSolver

 Aufgabe:
   Extra-Speicher-Konstruktor einer Instanz vom Typ LINEAR_SOLVER.

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

 Ergebnis:
   - Ungleich NULL im Erfolgsfall -

Programmaenderungen:
03/2000    AH      Erste Version

*************************************************************************/
LINEAR_SOLVER *InitMemoryLinearSolver(LINEAR_SOLVER * ls, int memory_number)
{
   int i;

   if (!ls || memory_number <= 0 || ls->dim <= 0)
      return NULL;
   SetLinearSolver(ls);

   if (ls->new_memory)
      DestroyMemoryLinearSolver(ls);

   ls->new_memory = (double **) Malloc(memory_number * sizeof(double *));
   if (ls->new_memory == NULL)
      return NULL;
   for (i = 0; i < memory_number; i++)
      ls->new_memory[i] = NULL;

   for (i = 0; i < memory_number; i++)
   {
      ls->new_memory[i] = (double *) Malloc(ls->dim * sizeof(double));
      if (ls->new_memory[i] == NULL)
         return DestroyMemoryLinearSolver(ls);
   }

   ls->memory_number = memory_number;
   return ls;
}


/*************************************************************************
 ROCKFLOW - Funktion: SetMemoryZeroLinearSolver

 Aufgabe:
   Setzt Speichervektor auf NULL.

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

 Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

Programmaenderungen:
08/1999     AH         Erste Version (ah neq)

*************************************************************************/
LINEAR_SOLVER *SetMemoryZeroLinearSolver(LINEAR_SOLVER * ls)
{
   if (!ls)
      return NULL;
   SetLinearSolver(ls);
   if (ls->xx)
      MNulleVec(ls->xx, ls->dim);

   return ls;
}


/*************************************************************************
 ROCKFLOW - Funktion: CreateLinearSolver

 Aufgabe:
   Konstruktor fuer LINEAR_SOLVER

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *name: Zeiger auf den Namen des LS's.

 Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

Programmaenderungen:
02/1999     AH         Erste Version
08/1999     AH         Erweiterung memory  (neq)
10/1999     AH         Systemzeit
10/1999     AH         Systemzeit
05/2000     OK         Erweiterung fuer vektorielle Groessen
08/2000     MK         Erweiterung fuer variable Dimension der Unbekannten
04/2003     MK         Erweiterung fuer variable Anzahl der Loeserobjekte

*************************************************************************/
LINEAR_SOLVER *CreateLinearSolver(long store, long dim)
{
   LINEAR_SOLVER *ls;

   ls = (LINEAR_SOLVER *) Malloc(sizeof(LINEAR_SOLVER));
   if (ls == NULL)
      return NULL;

   //OK    ls->name = NULL;
   ls->store = store;
   ls->dim = dim;
   ls->matrix = NULL;
   ls->b = NULL;
   ls->x = NULL;
   ls->lsp = NULL;
   ls->unknown_vector_dimension = 1;
   /* BC */
   //   ls->bc_names_dirichlet = NULL;
   //   ls->bc_names_neumann = NULL;
   ls->unknown_vector_indeces = NULL;
   ls->unknown_node_numbers = NULL;
   ls->unknown_update_methods = NULL;
   // wW   ls->bc_name = NULL;
   //    ls->bc_name2 = NULL;
   //    ls->bc_name3 = NULL;
   //    ls->boundary_conditions_function = NULL;
   /* SS */
   //    ls->sousin_name1 = NULL;
   //    ls->sousin_name2 = NULL;
   //    ls->sousin_name3 = NULL;
   //    ls->source_multiplicator_function = NULL;
   //    ls->source_multiplicator_function1 = NULL;
   //    ls->source_multiplicator_function2 = NULL;

   ls->lsp_name = NULL;

   ls->r = NULL;
   ls->xx = NULL;
   ls->memory = NULL;

   ls->memory_number = 0;
   ls->new_memory = NULL;

   ls->assemble_function = NULL;

   ls->system_time_name = NULL;
   ls->system_time_assemble_function_name = NULL;

   /* Multiple unknowns (Dim) / multiple solvers (MultSolvers) */
   ls->name_group_ls = NULL;
   ls->name_ls = NULL;
   ls->number_ls = -1;
   ls->num_of_unknowns_ls = 0;
   return ls;
}


/*************************************************************************
 ROCKFLOW - Funktion: CreateLinearSolverDim

 Aufgabe:
   Konstruktor fuer LINEAR_SOLVER

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *name: Zeiger auf den Namen des LS's.

 Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

Programmaenderungen:
02/1999     AH         Erste Version
08/1999     AH         Erweiterung memory  (neq)
10/1999     AH         Systemzeit
10/1999     AH         Systemzeit
05/2000     OK         Erweiterung fuer vektorielle Groessen
08/2000     MK         CreateLinearSolver->CreateLinearSolverDim
04/2003     MK         wird durch CreateLinearSolver wieder abgeloest
*************************************************************************/
LINEAR_SOLVER *CreateLinearSolverDim(long store, int unknown_vector_dimension, long dim)
{
   LINEAR_SOLVER *ls;

   ls = (LINEAR_SOLVER *) Malloc(sizeof(LINEAR_SOLVER));
   if (ls == NULL)
      return NULL;

   //OK    ls->name = NULL;
   ls->store = store;
   ls->dim = dim;
   ls->matrix = NULL;
   ls->b = NULL;
   ls->x = NULL;
   ls->lsp = NULL;
   ls->unknown_vector_dimension = unknown_vector_dimension;
   /* BC */
   /* //WW
   ls->bc_names_dirichlet = NULL;
   ls->bc_names_neumann = NULL;
   */
   ls->unknown_vector_indeces = NULL;
   ls->unknown_node_numbers = NULL;
   ls->unknown_update_methods = NULL;
   /* //WW
   ls->bc_name = NULL;
   ls->bc_name2 = NULL;
   ls->bc_name3 = NULL;
   ls->boundary_conditions_function = NULL;
   // SS
   ls->sousin_name1 = NULL;
   ls->sousin_name2 = NULL;
   ls->sousin_name3 = NULL;
   ls->source_multiplicator_function = NULL;
   ls->source_multiplicator_function1 = NULL;
   ls->source_multiplicator_function2 = NULL;
   */

   ls->lsp_name = NULL;

   ls->r = NULL;
   ls->xx = NULL;
   ls->memory = NULL;

   ls->memory_number = 0;
   ls->new_memory = NULL;

   ls->assemble_function = NULL;

   ls->system_time_name = NULL;
   ls->system_time_assemble_function_name = NULL;

   /* Multiple unknowns (Dim) / multiple solvers (MultSolvers) */
   ls->name_group_ls = NULL;
   ls->name_ls = NULL;
   ls->number_ls = -1;
   ls->num_of_unknowns_ls = 0;
   return ls;
}


//////////////////////////////////////////////////////////////////////////
// NUM
//////////////////////////////////////////////////////////////////////////

double GetNumericalTimeCollocation(char *name)
{
   CRFProcess* m_pcs = NULL;
   m_pcs = PCSGet(name);
   if(m_pcs)
      return m_pcs->m_num->ls_theta;
   else
      cout << "Fatal error in GetNumericalTimeCollocation: No valid PCS" << endl;
   //OK return np_array[GetNumericsIndex(name)].time_collocation;
   return 1.0;
}


/**************************************************************************/
/* ROCKFLOW - Funktion: CalcIterationError
 */
/* Aufgabe:
   Ermittelt den Fehler bei Iterationen

                                                                          */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *new_iteration : Vektor des neuen Iterationsschritts
   E double *old_iteration : Vektor des alten Iterationsschritts
   E double *reference     : Vektor des alten Zeitschritts (als Referenz)
   E longlength            : Laenge der Vektoren
   E int method            : Methode der Fehlerermittlung
                                                                          */
/* Ergebnis:
   - void -
                                                                          */
/* Programmaenderungen:
   1/1999     C.Thorenz  Zweite Version                                                                          */
/**************************************************************************/
double NUMCalcIterationError(double *new_iteration, double *old_iteration, double *reference, long length, int method)
{
   static long i;
   static double error, change, max_c, min_c;

   error = 0.;
   change = 0.;

   max_c = 0.;
   min_c = 1.e99;

   switch (method)
   {
      default:
      case 0:
         return 0.;
         /* Max. Unterschied zwischen altem und neuem Iterationsschritt */
      case 1:
         for (i = 0l; i < length; i++)
            error = max(error, fabs(new_iteration[i] - old_iteration[i]));
         return error;

         /* Max. Unterschied zwischen altem und neuem Iterationsschritt,
            jeweils normiert mit dem Mittelwert der Groesse des Wertes  */
      case 2:
         for (i = 0l; i < length; i++)
            error = max(error, 2. * fabs(new_iteration[i] - old_iteration[i]) / (fabs(new_iteration[i]) + fabs(old_iteration[i]) + MKleinsteZahl));
         return error;

         /* Max. Unterschied zwischen altem und neuem Iterationsschritt,
            normiert mit dem groessten Wert */
      case 3:
         for (i = 0l; i < length; i++)
         {
            error = max(error, fabs(new_iteration[i] - old_iteration[i]));
            max_c = max(max(max_c, fabs(new_iteration[i])),fabs(old_iteration[i]));
         }
         return error / (max_c + MKleinsteZahl);

         /* Max. Unterschied zwischen altem und neuem Iterationsschritt,
            normiert mit der Spanne der Werte */
      case 4:
         for (i = 0l; i < length; i++)
         {
            error = max(error, fabs(new_iteration[i] - old_iteration[i]));
            min_c = min(min_c, fabs(new_iteration[i]));
            max_c = max(max_c, fabs(new_iteration[i]));
         }
         return error / (max_c - min_c + MKleinsteZahl);

         /* Max. Unterschied zwischen altem und neuem Iterationsschritt,
            normiert mit dem Unterschied zum alten Zeitschritt. Die
            genaueste Methode, da die Fehlerberechnung dann Zeitschritt-
            unabhaengig wird! */
      case 5:
         for (i = 0l; i < length; i++)
            error = max(error, fabs(new_iteration[i] - old_iteration[i]) / (fabs(new_iteration[i] - reference[i]) + MKleinsteZahl));
         return error;

         /* Max. Unterschied zwischen altem und neuem Iterationsschritt,
            normiert mit dem maximalen Unterschied zum alten Zeitschritt */
      case 6:
         for (i = 0l; i < length; i++)
         {
            error = max(error, fabs(new_iteration[i] - old_iteration[i]));
            change = max(change, fabs(new_iteration[i] - reference[i]));
         }
         return error / (change + MKleinsteZahl);

         /* Der Vektorabstand */
      case 7:
         return MVekDist(old_iteration, new_iteration, length);

         /* Der Vektorabstand der Iteration, normiert mit dem Vektorabstand
            zur alten Zeitebene */
      case 8:
         return MVekDist(old_iteration, new_iteration, length) / (MVekDist(reference, new_iteration, length) + MKleinsteZahl);
   }
}
#endif                                            //ifndef NEW_EQS //WW. 06.11.2008

int GetNumericsGaussPoints(int element_dimension)
{
   int m_gaussian_points = 3;
   int g_gaussian_points = 3;
   switch (element_dimension)
   {
      case 1:
         m_gaussian_points = 1;
         break;
      case 2:
         m_gaussian_points = g_gaussian_points;
         break;
      case 3:
         m_gaussian_points = g_gaussian_points;
         break;
      case 4:
         m_gaussian_points = g_gaussian_points;
         break;
   }
   return m_gaussian_points;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
05/2005 OK Implementation
last modification:
**************************************************************************/
CNumerics* NUMGet(string num_name)
{
   CNumerics *m_num = NULL;
   for(int i=0;i<(int)num_vector.size();i++)
   {
      m_num = num_vector[i];
      if(m_num->pcs_type_name.compare(num_name)==0)
         return m_num;
   }
   return NULL;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
01/2005 OK Implementation
last modified:
**************************************************************************/
void NUMDelete()
{
   long i;
   int no_num =(int)num_vector.size();
   for(i=0;i<no_num;i++)
   {
      delete num_vector[i];
   }
   num_vector.clear();
}
