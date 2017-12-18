
#ifndef __DELAUNAY_2D_H_
#define __DELAUNAY_2D_H_

/* Modules */
#include <stdio.h>
#include "msh_lib.h"
/* Declarationen */
#define KTJ 1000000						//Max Number Of Node
#define KTE (2*KTJ+1)					//Max Number Of Element
#define KBD 20							//Max Number Of Boundary
#define KCM 200							//Max Number Of Element Around One Node
#define DELAUN2D_STACK_SIZE 0x10000000	//Stack Size
#define OPEN_BOUNDARY_POSTFIX "_OPENB"

/* functions /prototypes */
extern int ReadInputData(int *nex, int *nbk, int *nob, int *nib, int *ibex, int **ibno, int *nidm, int *ibreak, int **nbreak, double *px, double *py, double *pd);
extern int WriteOutputFile(char* filepath, int *node, double *x, double *y, int *nelm, int *mtj, int *idm);
extern int WriteOutputFileBIN(char* filepath, int *node, double *x, double *y, int *nelm, int *mtj, int *idm);
extern int ExecuteDelaunay2D(char* outfilepath);
extern int ExecuteDelaunay2DProcess(char* outfilepath);
//extern int Start_Delaunay2D (int argc, char** argv);
//extern int WriteOutputData(CFEMesh* m_msh, int *node, double *x, double *y, int *nelm, int *mtj, int *idm);
//extern int ReadInputFile(char* filepath, int *nex, int *nbk, int *nob, int *nib, int *ibex, int **ibno, int *nidm, int *ibreak, int **nbreak, double *px, double *py, double *pd);
//extern int WriteOutputFileRFI(char* filepath, int *node, double *x, double *y, int *nelm, int *mtj, int *idm);
//extern void ExecuteDelaunay2D(char* infilepath, char* outfilepath);

#endif


