#include "Stiff_Bulirsch-Stoer.h"
#include "stdlib.h"
#include "stdio.h"
#include <iostream>

/* interne Deklarationen */
void nrerror(char *error_text);
double *dvector(long nl, long nh);
int *ivector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
double xs, double htot, int nstep, double yout[],
void (*derivs)(double, double [], double [], int, long), long);
void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
double yscal[], double *hdid, double *hnext,
void (*derivs)(double, double [], double [], int, long));
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
double hmin, double *hnext, int *nok, int *nbad,
void (*derivs)(double, double [], double [], int, long),
void (*stifbs)(double [], double [], int, double *, double, double, double [],
double *, double *, void (*)(double, double [], double [], int, long)), long);

/*************************************************************************************/
/* Numerical Recipes Signum                                                          */
/*************************************************************************************/

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/*************************************************************************************/
/* Numerical Recipes Square                                                          */
/*************************************************************************************/

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

/*************************************************************************************/
/* Numerical Recipes double-max                                                      */
/*************************************************************************************/

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))

/*************************************************************************************/
/* Numerical Recipes double-min                                                      */
/*************************************************************************************/

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
(dminarg1) : (dminarg2))

/*************************************************************************************/
/* Numerical Recipes standard error handler                                           */
/*************************************************************************************/

void nrerror(const char *error_text)
//char error_text[];
{
   //void exit();

   fprintf(stderr,"\n \n Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   std::cout.flush();
   exit(1);
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* allocate a double vector with subscript range v[nl..nh]                           */
/*************************************************************************************/

double *dvector(long nl, long nh)
{
   double *v;

   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v)                                        //CB
   {
      std::cout <<  "allocation failure in dvector()" << std::endl;
      nrerror("allocation failure in dvector()");
   }
   return v-nl+NR_END;
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* allocate an int vector with subscript range v[nl..nh]                             */
/*************************************************************************************/

int *ivector(long nl, long nh)
{
   int *v;

   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v)                                        //CB
   {
      std::cout <<  "allocation failure in ivector()" << std::endl;
      nrerror("allocation failure in ivector()");
   }
   return v-nl+NR_END;
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]               */
/*************************************************************************************/

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   /* allocate pointers to rows */
   m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
   if (!m)                                        //CB
   {
      std::cout <<  "allocation failure 1 in matrix()" << std::endl;
      nrerror("allocation failure 1 in matrix()");
   }
   m += NR_END;
   m -= nrl;

   /* allocate rows and set pointers to them */
   m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
   if (!m[nrl])                                   //CB
   {
      std::cout <<  "allocation failure 2 in matrix()" << std::endl;
      nrerror("allocation failure 2 in matrix()");
   }
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   /* return pointer to array of pointers to rows */
   return m;
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* free a double vector allocated with dvector()                                     */
/*************************************************************************************/

void free_dvector(double *v, long nl, long nh)
{
   long i = nh; i++;                              // to avoid warning on nh
   free((FREE_ARG) (v+nl-NR_END));
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* free an int vector allocated with ivector()                                       */
/*************************************************************************************/

void free_ivector(int *v, long nl, long nh)
{
   long i = nh; i++;                              // to avoid warning on nh
   free((FREE_ARG) (v+nl-NR_END));
}


/*************************************************************************************/
/* Numerical Recipes                                                                 */
/* free a double matrix allocated by dmatrix()                                       */
/*************************************************************************************/

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
   long i = nch; i++;                             // to avoid warning on nh
   i=nrh;                                         // to avoid warning on nrh
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}


/*************************************************************************************/
/* Numerical Recipes Chapter 2.3                                                     */
/* LU Decomposition of matrix a[1..n][1..n]                                          */

/* a[1..n][1..n] input-Matrix und output-Matrix										 */
/* n input: Dimension des arrays a													 */
/* indx output: vector zeichnet row permutation auf (kann weg)						 */
/* d output:																		 */

/* Diese Funktion entspricht der Rockflow-Routine LU-Decomposition (solver.c),       */
/* indiziert aber ALLE arrays von 1-nn (Rockflow 0-(nn-1))                           */

/*************************************************************************************/

extern double **d,*x;

void ludcmp(double **a, int n, int *indx, double *d)
{
   int i,imax=0,j,k;                              //imax=0 SB warning avoided
   double big,dum,sum,temp;
   double *vv;

   vv=dvector(1,n);
   *d=1.0;

   for (i=1;i<=n;i++)
   {
      big=0.0;
      for (j=1;j<=n;j++)
         if ((temp=fabs(a[i][j])) > big) big=temp;
      if (big == 0.0)                             //CB
      {
         std::cout << "Singular matrix in routine ludcmp" << std::endl;
         nrerror("Singular matrix in routine ludcmp");
      }
      vv[i]=1.0/big;
   }

   for (j=1;j<=n;j++)
   {
      for (i=1;i<j;i++)
      {
         sum=a[i][j];
         for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
      }
      big=0.0;

      for (i=j;i<=n;i++)
      {
         sum=a[i][j];
         for (k=1;k<j;k++)
            sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
         if ( (dum=vv[i]*fabs(sum)) >= big)
         {
            big=dum;
            imax=i;
         }
      }

      if (j != imax)
      {
         for (k=1;k<=n;k++)
         {
            dum=a[imax][k];
            a[imax][k]=a[j][k];
            a[j][k]=dum;
         }
         *d = -(*d);
         vv[imax]=vv[j];
      }

      indx[j]=imax;
      if (a[j][j] == 0.0) a[j][j]=TINY;
      if (j != n)
      {
         dum=1.0/(a[j][j]);
         for (i=j+1;i<=n;i++) a[i][j] *= dum;
      }
   }

   free_dvector(vv,1,n);
}


/*************************************************************************************/
/* Numerical Recipes Chapter 2.3                                                     */
/* Solves Set of n linear equations                                                  */

/* Diese Funktion entspricht der Rockflow-Routine lubksb_3 (solver.c),               */
/* indiziert aber ALLE arrays von 1-nn (Rockflow 0-(nn-1))                           */

/*************************************************************************************/

void lubksb(double **a, int n, int *indx, double b[])
{
   int i,ii=0,ip,j;
   double sum;

   for (i=1;i<=n;i++)
   {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii)
         for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
      else if (sum) ii=i;
      b[i]=sum;
   }

   for (i=n;i>=1;i--)
   {
      sum=b[i];
      for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
      b[i]=sum/a[i][i];
   }
}


/*************************************************************************************/
/* Numerical Recipes Chapter 16.4                                                    */
/* Use polynominal extrapolation to evaluate nv functions at x=0 by fitting a        */
/* polynominal to a sequence of estimates with progressively smaller values x=xest   */
/*************************************************************************************/

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
   int k1,j;
   double q,f2,f1,delta,*c;

   c=dvector(1,nv);
   x[iest]=xest;
   for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];

   if (iest == 1)
   {
      for (j=1;j<=nv;j++) d[j][1]=yest[j];

   } else
   {
      for (j=1;j<=nv;j++) c[j]=yest[j];
      for (k1=1;k1<iest;k1++)
      {
         delta=1.0/(x[iest-k1]-xest);
         f1=xest*delta;
         f2=x[iest-k1]*delta;

         for (j=1;j<=nv;j++)
         {
            q=d[j][k1];
            d[j][k1]=dy[j];
            delta=c[j]-q;
            dy[j]=f1*delta;
            c[j]=f2*delta;
            yz[j] += dy[j];
         }
      }

      for (j=1;j<=nv;j++) d[j][iest]=dy[j];
   }

   free_dvector(c,1,nv);
}


/*************************************************************************************/
/* Numerical Recipes Chapter 16.6                                                    */
/* Performs one step of semi-implicit midpoint rule                                  */
/*************************************************************************************/

void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
double xs, double htot, int nstep, double yout[],
void (*derivs)(double, double [], double [], int, long), long node)
{
   int i,j,nn,*indx;
   double d,h,x,**a,*del,*ytemp;

   indx=ivector(1,n);
   a=dmatrix(1,n,1,n);
   del=dvector(1,n);
   ytemp=dvector(1,n);
   h=htot/nstep;

   for (i=1;i<=n;i++)
   {
      for (j=1;j<=n;j++) a[i][j] = -h*dfdy[i][j];
      ++a[i][i];
   }

   /* Rockflow; void ludcmp_3(double *a, long n, long *indx) */
   /* achtung, 1D/2D array*/
   ludcmp(a,n,indx,&d);

   for (i=1;i<=n;i++)
      yout[i]=h*(dydx[i]+h*dfdx[i]);

   /* Rockflow: void lubksb_3(double *a, long n, long *indx, double *b)*/
   /* achtung, 1D/2D array*/
   lubksb(a,n,indx,yout);

   for (i=1;i<=n;i++)
      ytemp[i]=y[i]+(del[i]=yout[i]);
   x=xs+h;

   /* calculate derivatives with external user-provided function */
   (*derivs)(x,ytemp,yout,n,node);

   for (nn=2;nn<=nstep;nn++)
   {
      for (i=1;i<=n;i++)
         yout[i]=h*yout[i]-del[i];

      /* Rockflow: void lubksb_3(double *a, long n, long *indx, double *b)*/
      /* achtung, 1D/2D array*/
      lubksb(a,n,indx,yout);
      for (i=1;i<=n;i++)
         ytemp[i] += (del[i] += 2.0*yout[i]);
      x += h;
      (*derivs)(x,ytemp,yout,n,node);
   }
   for (i=1;i<=n;i++)
      yout[i]=h*yout[i]-del[i];

   /* Rockflow: void lubksb_3(double *a, long n, long *indx, double *b)*/
   /* achtung, 1D/2D array*/
   lubksb(a,n,indx,yout);
   for (i=1;i<=n;i++)
      yout[i] += ytemp[i];

   free_dvector(ytemp,1,n);
   free_dvector(del,1,n);
   free_dmatrix(a,1,n,1,n);
   free_ivector(indx,1,n);
}


/*************************************************************************************/
/* Numerical Recipes Chapter 16.6                                                    */
/* Semi-implicit extrapolation step for integrating stiff ODEs                       */
/*************************************************************************************/

#define KMAXX 7
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define SCALMX 0.1

double **d,*x;

void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
double yscal[], double *hdid, double *hnext,
void (*derivs)(double, double [], double [], int, long), long node)
{
   int i,iq,k,kk,km;
   static int first=1,kmax,kopt,nvold = -1;
   static double epsold = -1.0,xnew;
   double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
   double *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;
   static double a[IMAXX+1];
   static double alf[KMAXX+1][KMAXX+1];
   static int nseq[IMAXX+1]={0,2,6,10,14,22,34,50,70};
   //	static int nseq[IMAXX+1]={0,2,6,10,14,22,34,50,70,98,138,194,274,386};
   // Num Recip S. 744
   // Differenz zwischen zwei Werten muss ein Vielfaches von 4 sein
   // So wählen, dass Verhältnis der Werte <= 5/7 ist
   // z.B. 10, 14 -> 10/14 <= 5/7
   // nächster wäre 18, aber 14/18 > 5/7, deshalb 22 mit 14/22 <= 5/7
   int reduct,exitflag=0;
   km=0;                                          //SB avoid warnings
   red = 0.0;                                     //SB avoid warning
   errmax = 0.0;                                  // SB avoid warning
   scale = 0.0;                                   // SB avoid warning

   d=dmatrix(1,nv,1,KMAXX);
   dfdx=dvector(1,nv);
   dfdy=dmatrix(1,nv,1,nv);
   err=dvector(1,KMAXX);
   x=dvector(1,KMAXX);
   yerr=dvector(1,nv);
   ysav=dvector(1,nv);
   yseq=dvector(1,nv);
   if(eps != epsold || nv != nvold)
   {
      *hnext = xnew = -1.0e29;
      eps1=SAFE1*eps;
      a[1]=nseq[1]+1;
      for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
      for (iq=2;iq<=KMAXX;iq++)
      {
         for (k=1;k<iq;k++)
            alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
               ((a[iq+1]-a[1]+1.0)*(2*k+1))));
      }
      epsold=eps;
      nvold=nv;
      a[1] += nv;
      for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
      for (kopt=2;kopt<KMAXX;kopt++)
         if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
      kmax=kopt;
   }
   h=htry;
   for (i=1;i<=nv;i++) ysav[i]=y[i];

   jacobn(*xx,y,dfdx,dfdy,nv,node);

   if (*xx != xnew || h != (*hnext))
   {
      first=1;
      kopt=kmax;
   }
   reduct=0;

   for (;;)
   {
      for (k=1;k<=kmax;k++)
      {
         xnew=(*xx)+h;
         if (xnew == (*xx))                       //CB
         {
            std::cout << "step size underflow in stifbs" << std::endl;
            nrerror("step size underflow in stifbs");
         }
         simpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs,node);
         xest=DSQR(h/nseq[k]);

         pzextr(k,xest,yseq,y,yerr,nv);

         if (k != 1)
         {
            errmax=TINY;
            for (i=1;i<=nv;i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
            errmax /= eps;
            km=k-1;
            err[km]=pow(errmax/SAFE1,1.0/(double)(2*km+1));
         }

         if (k != 1 && (k >= kopt-1 || first))
         {
            if (errmax < 1.0)
            {
               exitflag=1;
               break;
            }
            if (k == kmax || k == kopt+1)
            {
               red=SAFE2/err[km];
               break;
            }
            else if (k == kopt && alf[kopt-1][kopt] < err[km])
            {
               red=1.0/err[km];
               break;
            }
            else if (kopt == kmax && alf[km][kmax-1] < err[km])
            {
               red=alf[km][kmax-1]*SAFE2/err[km];
               break;
            }
            else if (alf[km][kopt] < err[km])
            {
               red=alf[km][kopt-1]/err[km];
               break;
            }
         }
      }
      //		if (exitflag) std::cout << " Exitflag > 0 in stifbs of biodegradation" << std::endl;
      if (exitflag) break;

      red=DMIN(red,REDMIN);
      red=DMAX(red,REDMAX);
      h *= red;
      reduct=1;
   }
   *xx=xnew;
   *hdid=h;
   first=0;
   wrkmin=1.0e35;
   for (kk=1;kk<=km;kk++)
   {
      fact=DMAX(err[kk],SCALMX);
      work=fact*a[kk+1];
      if (work < wrkmin)
      {
         scale=fact;
         wrkmin=work;
         kopt=kk+1;
      }
   }
   *hnext=h/scale;
   if (kopt >= k && kopt != kmax && !reduct)
   {
      fact=DMAX(scale/alf[kopt-1][kopt],SCALMX);
      if (a[kopt+1]*fact <= wrkmin)
      {
         *hnext=h/fact;
         kopt++;
      }
   }

   free_dvector(yseq,1,nv);
   free_dvector(ysav,1,nv);
   free_dvector(yerr,1,nv);
   free_dvector(x,1,KMAXX);
   free_dvector(err,1,KMAXX);
   free_dmatrix(dfdy,1,nv,1,nv);
   free_dvector(dfdx,1,nv);
   free_dmatrix(d,1,KMAXX,1,KMAXX);
}


#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef SCALMX
#undef NRANSI

/*************************************************************************************/
/* Numerical Recipes Chapter 16.2                                                    */
/* Driver for Runge-Kutta Algorithm fitted to stiff Burlish-Stoer Method (stifbs)    */
/* Integrates starting values ystart[1..nvar] from x1 to x2 using an adaptive step   */
/* size control                                                                      */
/*************************************************************************************/

/*extern int kmax,kount;
extern double *xp,**yp,dxsav;
*/
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
double hmin, double *nexth, int *nok, int *nbad,
void (*derivs)(double, double [], double [], int, long),
void (*stifbs)(double [], double [], int, double *, double, double, double [],
double *, double *, void (*)(double, double [], double [], int, long), long), long node)
{
   int nstp,i;
   /*	double xsav,x,hnext,hdid,h;
    */
   double x,hdid,hnext,h;
   double *yscal,*y,*dydx;

   yscal=dvector(1,nvar);
   y=dvector(1,nvar);
   dydx=dvector(1,nvar);
   x=x1;
   h=SIGN(h1,x2-x1);

   /*	*nok = (*nbad) = kount = 0;
    */
   *nok = (*nbad) = 0;
   for (i=1;i<=nvar;i++) y[i]=ystart[i];

   /*	if (kmax > 0) xsav=x-dxsav*2.0;
    */
   for (nstp=1;nstp<=MAXSTEP;nstp++)
   {

      /* calculate derivatives with external user-provided function */
      (*derivs)(x,y,dydx,nvar,node);

      /* Preconditioning */
      /*			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY; */
      for (i=1;i<=nvar;i++)
         yscal[i]=DMAX(1.0,fabs(y[i]));

      /*		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
               xp[++kount]=x;
               for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
               xsav=x;
            }
      */
      if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;

      (*stifbs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs, node);
      if (hdid == h) ++(*nok); else ++(*nbad);
      if ((x-x2)*(x2-x1) >= 0.0)
      {
         for (i=1;i<=nvar;i++) ystart[i]=y[i];
         /*			if (kmax) {
                     xp[++kount]=x;
                     for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
                  }
         */
         free_dvector(dydx,1,nvar);
         free_dvector(y,1,nvar);
         free_dvector(yscal,1,nvar);
         *nexth=hnext;
         return;
      }
      if (fabs(hnext) <= hmin)                    //CB
      {
         std::cout << "Step size too small in odeint" << std::endl;
         nrerror("Step size too small in odeint");
      }
      h=hnext;
   }
                                                  //CB
   std::cout << "Too many steps in routine odeint" << std::endl;
   nrerror("Too many steps in routine odeint");
}


#undef MAXSTEP
#undef TINY
