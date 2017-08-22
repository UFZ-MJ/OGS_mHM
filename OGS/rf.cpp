/**************************************************************************/
/* ROCKFLOW - Modul: rf.c
                                                                          */
/* Aufgabe:
   ROCKFLOW-FEM - Hauptprogramm
                                                                          */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   06/1998     AH         Konfigurationsdatei
   08/1999     OK         RF-FEM Applikation
   10/1999     AH         Systemzeit

   last modified: OK 14.12.1999
                                                                          */
/**************************************************************************/

/**
 * the preprocessor directive RFW_FRACTURE is only useable until version 4.11 of OGS
 * */

#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS)
#include <mpi.h>
#include "par_ddc.h"
#endif
#ifdef LIS
#include <omp.h>
#include "lis.h"
#endif

#include "BuildInfo.h"

/* Preprozessor-Definitionen */
#include "makros.h"
#define TEST
/* Benutzte Module */
#include "break.h"
#include "timer.h"
//16.12.2008. WW #include "rf_apl.h"
#include "files0.h"
#include "FileTools.h"
#ifdef SUPERCOMPUTER
// kg44 test for buffered outputh
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#endif
#include "problem.h"

Problem *aproblem = NULL;
/* Deklarationen */
int main ( int argc, char *argv[] );
void ShowSwitches ( void );
// LB,string FileName; //WW
// LB,string FilePath; //23.02.2009. WW
// ------  12.09.2007 WW:
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS)
double elapsed_time_mpi;
// ------
#endif
/* Definitionen */

/**************************************************************************/
/* ROCKFLOW - Funktion: main
                                                                          */
/* Aufgabe:
   Hauptprogramm
                                                                          */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int argc: Anzahl der Kommandozeilenparameter (incl. Programmname)
   E char *argv[]: Zeiger auf Feld der argc Kommandozeilenparameter
                                                                          */
/* Ergebnis:
   Fehlerfreie Bearbeitung: Exit-Code 0
                                                                          */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   08/1999     OK         RF-FEM Applikation
                                                                          */
/**************************************************************************/
int main ( int argc, char *argv[] )
{
  /* parse command line arguments */
  std::string anArg;
  std::string modelRoot;
  for( int i = 1; i < argc; i++ ) {
    anArg = std::string( argv[i] );
    if( anArg == "--help" || anArg == "-h")
      {
	std::cout << "Usage: ogs [MODEL_ROOT] [OPTIONS]\n"
		  << "Where OPTIONS are:\n"
		  << "  -h [--help]       print this message and exit\n"
		  << "  -b [--build-info] print build info and exit\n"
		  << "  --version         print ogs version and exit" << std::endl;
	continue;
      }
    if( anArg == "--build-info" || anArg == "-b" )
      {
		  std::cout << "ogs version: " << OGS_VERSION << std::endl
		  << "ogs date: " << OGS_DATE << std::endl
		  << "cmake command line arguments: " << CMAKE_CMD_ARGS << std::endl
		  << "git commit info: " << GIT_COMMIT_INFO << std::endl
		  << "subversion info: " << SVN_REVISION << std::endl
		  << "build timestamp: " << BUILD_TIMESTAMP << std::endl;
		  continue;
      }
    if( anArg == "--version" )
      {
		  std::cout << OGS_VERSION << std::endl;
		  continue;
      }
    if( anArg == "--model-root" || anArg == "-m" )
      {
		modelRoot = std::string( argv[++i] );
		continue;
      }
    // anything left over must be the model root, unless already found
    if ( modelRoot == "" ){ modelRoot = std::string( argv[i] ); }
  } // end of parse argc loop

  if( argc > 1 && modelRoot == "" ) // non-interactive mode and no model given
    exit(0);                         // e.g. just wanted the build info

  char *dateiname;
#ifdef SUPERCOMPUTER
// *********************************************************************
// buffered output ... important for performance on cray
// (unbuffered output is limited to 10 bytes per second)
// georg.kosakowski@psi.ch 11.10.2007

  char buf[1024*1024];
  int bsize;

	 bsize = 1024*1024; // question: what happens if buffer is full?
                         // according to documentation the buffer is flushed when full.
                         // If we have a lot of output, increasing buffer is usefull.
     if(bsize > 0) {
//        bufstd = malloc(bsize);
        setvbuf(stdout, buf, _IOFBF, bsize);
     }
//**********************************************************************
#endif
/*---------- MPI Initialization ----------------------------------*/
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS)
      printf("Before MPI_Init\n");
      MPI_Init(&argc,&argv);
      MPI_Barrier (MPI_COMM_WORLD); // 12.09.2007 WW
      elapsed_time_mpi = -MPI_Wtime(); // 12.09.2007 WW
      MPI_Comm_size(MPI_COMM_WORLD,&mysize);
      MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	  std::cout << "After MPI_Init myrank = " << myrank << '\n';
      time_ele_paral = 0.0;
#endif
/*---------- MPI Initialization ----------------------------------*/
/*---------- LIS solver -----------------------------------------*/
#ifdef LIS
  //int argc=0;
  //char** argv = NULL;
  // Initialization of the lis solver.
//  lis_initialize(&argc, &argv);	PCH: Undoing NW modification for compilation OGS-5
#endif
/*========================================================================*/
/* Kommunikation mit Betriebssystem */
  /* Ctrl-C ausschalten */
  NoBreak();
  /* Timer fuer Gesamtzeit starten */
#ifdef TESTTIME
    TStartTimer(0);
#endif
  /* Intro ausgeben */
#if defined(USE_MPI) //WW
   if(myrank==0)
#endif
  DisplayStartMsg();
  /* Speicherverwaltung initialisieren */
  if (!InitMemoryTest()) {
    DisplayErrorMsg("Fehler: Speicherprotokoll kann nicht erstellt werden!");
    DisplayErrorMsg("        Programm vorzeitig beendet!");
    return 1;	// LB changed from 0 to 1 because 0 is indicating success
  }
  if( argc == 1 )                     // interactive mode
    {
      dateiname = ReadString();
    }
  else                               // non-interactive mode
    {
      if ( argc == 2 )               // a model root was supplied 
	{
	  dateiname = (char *) Malloc((int)strlen(argv[1])+1);
	  dateiname = strcpy(dateiname,argv[1]);
	}
      else                          // several args supplied
	if( modelRoot != "")
	  {
	    dateiname = (char *) Malloc( (int) modelRoot.size() + 1 );
	    dateiname = strcpy( dateiname, modelRoot.c_str() );
	  }
      DisplayMsgLn(dateiname);
  }
  //WW  DisplayMsgLn("");
  //WW  DisplayMsgLn("");
  // ----------23.02.2009. WW-----------------

  // LB Check if file exists
  std::string tmpFilename = dateiname;
  tmpFilename.append(".pcs");
  if(!IsFileExisting(tmpFilename))
  {
	  std::cout << " Error: Cannot find file " << dateiname << std::endl;
	  return 1;
  }

  FileName = dateiname;
  size_t indexChWin, indexChLinux;
  indexChWin = indexChLinux = 0;
  indexChWin = FileName.find_last_of('\\');
  indexChLinux = FileName.find_last_of('/');
  //
  if(indexChWin!=std::string::npos)
     FilePath = FileName.substr(0,indexChWin)+"\\";
  else if(indexChLinux!=std::string::npos)
     FilePath = FileName.substr(0,indexChLinux)+"/";
  // ---------------------------WW
  aproblem = new Problem(dateiname);
  aproblem->Euler_TimeDiscretize();
  delete aproblem;
  aproblem = NULL;
#ifdef TESTTIME
  std::cout << "Simulation time: " << TGetTimer(0) << "s" << std::endl;
#endif
  /* Abspann ausgeben */
  /* Ctrl-C wieder normal */
  StandardBreak();
/*--------- MPI Finalize ------------------*/
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
   elapsed_time_mpi += MPI_Wtime(); // 12.09.2007 WW
   std::cout<<"\n *** Total CPU time of parallel modeling: "<< elapsed_time_mpi<<std::endl; //WW
   // Count CPU time of post time loop WW
   MPI_Finalize();
#endif
/*--------- MPI Finalize ------------------*/
/*--------- LIS Finalize ------------------*/
#ifdef LIS
//  lis_finalize();	//PCH: Undoing NW modification for compilation OGS-5
#endif
/*--------- LIS Finalize ------------------*/

  free(dateiname);
  return 0;
}
