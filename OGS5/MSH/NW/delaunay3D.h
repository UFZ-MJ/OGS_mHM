
/* Modules */
#include <stdio.h>


/* Declarationen */
#define NODE_MAX 10000				//MAX NUMBER OF NODE
#define ELEM_MAX 20000				//MAX NUMBER OF ELEMENT
#define DELAUN3D_STACK_SIZE 0		//STACK SIZE

/* functions /prototypes */
extern int ReadInputFile(char* filepath, int *node, double *x, double *y, double *z);
extern int WriteOutputFile(char* filepath, int *node, double *x, double *y, double *z, int *nelm, int *mtj, int *jac, int *maxelm);
extern int ExecuteDelaunay3D(char* infilepath, char* outfilepath);
extern int ExecuteDelaunay3DProcess(char* infilepath, char* outfilepath);
extern unsigned __stdcall ExecuteDelaunay3DThread(void* arg);





