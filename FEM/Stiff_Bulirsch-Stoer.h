
/* Deklarationen */

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#define MAXSTEP 20000
#define TINY 1.0e-30
#define NR_END 0
#define FREE_ARG char*

/* externe Funktionen */
extern void derivs(double x, double y[], double dydx[], int n, long node);
extern void jacobn(double x, double y[], double dfdx[], double **dfdy, int n, long node);
