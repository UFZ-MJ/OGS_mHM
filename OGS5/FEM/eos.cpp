/**********************************************************************
Module: Equation of State

Task: This file includes coefficients and functions for calculating the
thermal properties of liquids and gases in relation to density, pressure
and temperature.

Programming: NB
          Aug 2008
**********************************************************************/

#include "rf_mmp_new.h"
#include "rf_mfp_new.h"
//#include "rf_num_new.h"
#include "eos.h"
#include "tools.h"

using namespace std;

/**********************************************************************
Some functions and derivations as shown in [Span&Wagner 1994]
***********************************************************************/
double theta_fn (double tau, double A, double delta, double beta)
{
   return (1-tau)+A*pow(pow(delta-1,2),1/(2*beta));
}


double phi_fn (double C, double delta, double D, double tau)
{
   return exp (-C*pow(delta-1,2)-D*pow(tau-1,2));
}


double delta_fn (double theta_fn, double B, double delta, double alpha)
{
   return pow(theta_fn,2)+B*pow(pow(delta-1,2),alpha);
}


double dDELTApowb_ddelta (double b, double delta_fn, double dDELTA_deriv)
{
   return b*pow(delta_fn,b-1)*dDELTA_deriv;
}


double dDELTA_ddelta (double delta, double A, double theta_fn, double beta, double B, double a)
{
   return   ((delta-1)*(A*theta_fn*2/beta*pow(pow((delta-1),2),(1/(2*beta)-1))+2*B*a*pow(pow((delta-1),2),(a-1))));
}


double d2DELTA_ddelta2 (double delta, double dDELTA_deriv, double A, double theta_fn, double beta, double B, double a)
{
   return 1/(delta-1)*dDELTA_deriv+pow((delta-1),2)*
      (4*B*a*(a-1)*pow(pow((delta-1),2),(a-2))+2*pow(A,2)*pow((1/beta),2)*pow(pow(pow((delta-1),2),(1/(2*beta)-1)),2)
      +A*theta_fn*4/beta*(1/(2*beta)-1)*pow(pow((delta-1),2),(1/(2*beta)-2)));
}


double d2DELTApowb_ddelta2 (double b,double delta_fn,double d2DELTA_deriv,double dDELTA_deriv)
{
   return b*(pow(delta_fn,(b-1))*d2DELTA_deriv+(b-1)*pow(delta_fn,(b-2))*pow(dDELTA_deriv,2));
}


double dphi_ddelta (double C,double delta, double phi)
{
   return -2*C*(delta-1)*phi;
}


double dphi_dtau (double D,double tau, double phi_fn)
{
   return -2*D*(tau-1)*phi_fn;
}


double d2phi_dtau2 (double D,double tau, double phi_fn)
{
   return (2*D*pow((tau-1),2)-1)*2*D*phi_fn;
}


double d2phi_ddelta2 (double C,double delta, double phi_fn)
{
   return (2*C*pow((delta-1),2)-1)*2*C*phi_fn;
}


double dDELTA_dtau (double theta_fn, double b, double delta_fn)
{
   return -2*theta_fn*b*pow(delta_fn,(b-1));
}


double d2DELTA_dtau2 (double b, double delta_fn, double theta_fn)
{
   return 2*b*pow(delta_fn,(b-1))+4*pow(theta_fn,2)*b*(b-1)*pow(delta_fn,(b-2));
}


double d2DELTApowb_ddeltadtau (double A,double b,double beta,double delta_fn,double delta,double theta_fn,double dDELTA_deriv)
{
   return -A*b*2/beta*pow(delta_fn,(b-1))*(delta-1)*pow(pow((delta-1),2),(1/(2*beta)-1))
      -2*theta_fn*b*(b-1)*pow(delta_fn,(b-2))*dDELTA_deriv;
}


double d2phi_ddeltadtau (double C,double D,double delta,double tau,double phi_fn)
{
   return 4*C*D*(delta-1)*(tau-1)*phi_fn;
}


/**********************************************************************
A derivation of the free energy function phi
last change: NB JUN 09
***********************************************************************/
double CFluidProperties::phi_r_d (double rho, double T, int c)
{
   c = c;                                         //OK411

   double phi_a=0,phi_b=0,phi_c=0,phi_d=0;
   double delta,tau,DELTA,THETA,PHI,DPHI,dDELTA_deriv,dDELTApowbddelta;
   int i;

   tau = Tc/T;
   delta = rho/rhoc;

   for (i=0; i<limit[3]; i++)
   {
      if (i<limit[0])

         phi_a=phi_a+(K[0][i]*K[1][i]*pow(delta,(K[1][i]-1))*pow(tau,K[2][i]));

      else if (i<limit[1])

         phi_b=phi_b+(K[0][i]*exp(-pow(delta,K[3][i]))*(pow(delta,K[1][i]-1)*
               pow(tau,K[2][i])*(K[1][i]-K[3][i]*pow(delta,K[3][i]))));

      else if (i<limit[2])

         phi_c=phi_c+(K[0][i]*pow(delta,K[1][i])*pow(tau,K[2][i])*
               exp(-K[10][i]*pow(delta-K[13][i],2)-K[11][i]*pow(tau-K[12][i],2))*(K[1][i]
               /delta-2*K[10][i]*(delta-K[13][i])));

      else if (i<limit[3])
      {

         THETA  = theta_fn (tau,K[6][i],delta,K[11][i]);
         DELTA  = delta_fn (THETA,K[7][i],delta,K[4][i]);
         PHI    = phi_fn (K[8][i],delta,K[9][i],tau);
         dDELTA_deriv = dDELTA_ddelta(delta,K[6][i],THETA,K[11][i],K[7][i],K[4][i]);
         dDELTApowbddelta = dDELTApowb_ddelta(K[5][i],DELTA,dDELTA_deriv);
         DPHI   = dphi_ddelta (K[8][i],delta,PHI);

         phi_d=phi_d+K[0][i]*(pow(DELTA,K[5][i])*(PHI+delta*DPHI)+dDELTApowbddelta*delta*PHI);

      }
   }
   return phi_a+phi_b+phi_c+phi_d;
}


/**********************************************************************
A derivation of the free energy function phi
last change: NB JUN 09
***********************************************************************/
double CFluidProperties::phi_r_tt (double rho, double T, int c)
{
   c = c;                                         //OK411
   //CFluidProperties *FP;
   double phi_a=0,phi_b=0,phi_c=0,phi_d=0;
   double delta,tau,THETA,PHI,DELTA,DDELTA,D2DELTA,DPHI,D2PHI;
   int i=0;

   //FP=MFPGet(c);

   tau = Tc/T;
   delta = rho/rhoc;

   for (i=0; i<limit[3]; i++)
   {
      if (i<limit[0])

         phi_a = phi_a+ (K[0][i]*K[2][i]*(K[2][i]-1)*pow(delta,K[1][i])*
            pow(tau,(K[2][i]-2)));

      else if (i<limit[1])

         phi_b = phi_b+ (K[0][i]*K[2][i]*(K[2][i]-1)*
               pow(delta,K[1][i])*pow(tau,(K[2][i]-2))*exp(-pow(delta,K[3][i])));

      else if (i<limit[2])

         phi_c = phi_c+ (K[0][i]*pow(delta,K[1][i])*pow(tau,K[2][i])*
               exp(-K[10][i]*pow((delta-K[13][i]),2)-K[11][i]*pow((tau-K[12][i]),2))
               *(pow((K[2][i]/tau-2*K[11][i]*(tau-K[12][i])),2)-K[2][i]
               /pow(tau,2)-2*K[11][i]));

      else if (i<limit[3])
      {

         THETA  = theta_fn (tau,K[6][i],delta,K[11][i]);
         DELTA  = delta_fn (THETA,K[7][i],delta,K[4][i]);
         PHI    = phi_fn (K[8][i],delta,K[9][i],tau);

         D2DELTA = d2DELTA_dtau2 (K[5][i],DELTA,THETA);
         DDELTA  = dDELTA_dtau (THETA, K[5][i], DELTA);
         DPHI    = dphi_dtau (K[9][i],tau,PHI);
         D2PHI   = d2phi_dtau2 (K[9][i],tau,PHI);

         phi_d = phi_d + (K[0][i]*delta*(D2DELTA*PHI+2*DDELTA*DPHI+pow(DELTA,K[5][i])*D2PHI));
      }
   }
   return phi_a+phi_b+phi_c+phi_d;
}


/**********************************************************************
A derivation of the free energy function phi
last change: NB JUN 09
***********************************************************************/
double CFluidProperties::phi_0_t (double T,int c)
{
   c = c;                                         //OK411
   double phi_c=0,phi_d=0,phi_e=0;
   double tau;
   int i;

   tau = Tc/T;

   phi_c = k[0][1];
   phi_d = k[0][2]/tau;

   for (i=3; i<8; i++) phi_e= phi_e + (k[0][i]*k[1][i]*(1/(1-exp(-k[1][i]*tau))-1));

   return phi_c+phi_d+phi_e;
}


/**********************************************************************
A derivation of the free energy function phi
last change: NB JUN 09
***********************************************************************/
double CFluidProperties::phi_0_tt (double T, int c)
{
   c = c;                                         //OK411
   double phi_d=0,phi_e=0;
   double tau;
   int i;

   tau = Tc/T;
   phi_d = k[0][2]/pow(tau,2);
   for (i=3; i<8; i++) phi_e=phi_e+(k[0][i]*pow(k[1][i],2)*exp(-k[1][i]*tau)*pow(1-exp(-k[1][i]*tau),-2));

   return 0-phi_d-phi_e;
}


/**********************************************************************
A derivation of the free energy function phi
last change: NB JUN 09
***********************************************************************/
double CFluidProperties::phi_r_t (double rho, double T,int c)
{
   c = c;                                         //OK411

   double phi_a=0,phi_b=0,phi_c=0,phi_d=0,h;
   int i;
   double delta,tau;
   double thetafn,deltafn,ddeltatau,phifn,dphitau;

   tau = Tc/T;
   delta = rho/rhoc;

   for (i=0;i<limit[0];i++)
   {
      phi_a = phi_a + (K[0][i]*K[2][i]*pow(delta,K[1][i])*pow(tau,(K[2][i]-1)));
   }

   for  (i=limit[0];i<limit[1];i++)
   {
      phi_b = phi_b + (K[0][i]*K[2][i]*pow(delta,K[1][i])*pow(tau,(K[2][i]-1))*exp(-pow(delta,K[3][i])));
   }

   for (i=limit[1];i<limit[2];i++)
   {
      phi_c = phi_c + (K[0][i]*pow(delta,K[1][i])*pow(tau,(K[2][i]))*
         exp(-K[10][i]*pow((delta-K[13][i]),2)-K[11][i]*pow((tau-K[12][i]),2))*
         (K[2][i]/tau-2*K[11][i]*(tau-K[12][i])));

   }

   for (i=limit[2];i<limit[3];i++)
   {
      thetafn = theta_fn(tau,K[6][i],delta,K[11][i]);
      deltafn = delta_fn(thetafn,K[7][i],delta,K[4][i]);
      ddeltatau = dDELTA_dtau(thetafn,K[5][i],deltafn);
      phifn = phi_fn(K[8][i],delta,K[9][i],tau);
      dphitau = dphi_dtau(K[9][i],tau,phifn);

      phi_d = phi_d + (K[0][i]*delta*(ddeltatau*phifn+pow(deltafn,K[5][i])*dphitau));

   }
   h = phi_a+phi_b+phi_c+phi_d;
   return h;

}


/**********************************************************************
A derivation of the free energy function phi
last change: NB 4.9.05
***********************************************************************/
double CFluidProperties::phi_r_dt (double rho, double T, int c)
{
   c = c;                                         //OK411
   double phi_a=0,phi_b=0,phi_c=0,phi_d=0;
   int i;
   double delta,tau;
   double phifn,thetafn,deltafn,d2phideriv,dDELTAderiv,dDELTApowbddelta,dphidtau,
      dDELTApowbdtau,dphiddelta,d2DELTApowbddeltadtau;

   tau = Tc/T;
   delta = rho/rhoc;

   for (i=0;i<limit[0];i++)
   {
      phi_a = phi_a + (K[0][i]*K[1][i]*K[2][i]*pow(delta,(K[1][i]-1))*pow(tau,(K[2][i]-1)));

   }
   for (i=limit[0];i<limit[1];i++)
   {
      phi_b = phi_b + (K[0][i]*K[2][i]*pow(delta,(K[1][i]-1))*pow(tau,(K[2][i]-1))*(K[1][i]-K[3][i]*
         pow(delta,K[3][i]))*exp(-pow(delta,K[3][i])));

   }
   for (i=limit[1];i<limit[2];i++)
   {

      phi_c = phi_c + ((K[0][i]*pow(delta,K[1][i])*pow(tau,K[2][i])*exp(
         -K[10][i]*pow((delta-K[13][i]),2)-K[11][i]*pow((  tau-K[12][i]),2)))*

         (K[1][i]/delta-2*K[10][i]*(delta-K[13][i]))*
         (K[2][i]/tau-2*K[11][i]*(tau-K[12][i])));

   }
   for (i=limit[2];i<limit[3];i++)
   {
      phifn = phi_fn(K[8][i],delta,K[9][i],tau);
      thetafn = theta_fn(tau,K[6][i],delta,K[11][i]);
      deltafn = delta_fn(thetafn,K[7][i],delta,K[4][i]);
      d2phideriv = d2phi_ddeltadtau(K[8][i],K[9][i],delta,tau,phifn);

      dDELTAderiv = dDELTA_ddelta(delta,K[6][i],thetafn,K[11][i],K[7][i],K[4][i]);

      dDELTApowbddelta = dDELTApowb_ddelta(K[5][i],deltafn,dDELTAderiv);
      dphidtau = dphi_dtau(K[9][i],tau,phifn);
      dDELTApowbdtau = dDELTA_dtau(thetafn,K[5][i],deltafn);
      dphiddelta = dphi_ddelta(K[8][i],delta,phifn);
      d2DELTApowbddeltadtau = d2DELTApowb_ddeltadtau(K[6][i],K[5][i],K[11][i],deltafn,delta,thetafn,dDELTAderiv);

      phi_d = phi_d + (K[0][i]*(pow(deltafn,K[5][i])*(dphidtau+delta*d2phideriv)+
         delta*dDELTApowbddelta*dphidtau+dDELTApowbdtau*(phifn+delta*dphiddelta)+d2DELTApowbddeltadtau*
         delta*phifn));
   }

   return phi_a+phi_b+phi_c+phi_d;
}


/**********************************************************************
A derivation of the free energy function phi
last change: NB 4.9.05
***********************************************************************/
double CFluidProperties::phi_r_dd (double rho, double T, int c)
{
   c = c;                                         //OK411
   double phi_a=0,phi_b=0,phi_c=0,phi_d=0;
   int i;
   double delta,tau;
   double thetafn,deltafn,phifn,dphiddelta,d2phiddelta2,dDELTA_deriv,dDELTApowbddelta,d2DELTA_deriv,d2DELTApowbddelta2;

   tau = Tc/T;
   delta = rho/rhoc;

   for (i=0;i<limit[0];i++)
   {
      phi_a = phi_a + (K[0][i]*K[1][i]*(K[1][i]-1)*pow(delta,(K[1][i]-2))*pow(tau,K[2][i]));
   }
   for (i=limit[0];i<limit[1];i++)
   {
      phi_b = phi_b+((K[0][i]*exp(-(pow(delta,K[3][i]))))*(((pow(delta,(K[1][i]-2))*
         pow(tau,K[2][i])))*(((K[1][i]-K[3][i]*pow(delta,K[3][i])))*(K[1][i]-1-K[3][i]*
         pow(delta,K[3][i]))-(pow(K[3][i],2)*pow(delta,K[3][i])))));

   }
   for (i=limit[1];i<limit[2];i++)
   {
      phi_c = phi_c + ((K[0][i]*pow(tau,K[2][i]))*
         exp(-K[10][i]*pow((delta-K[13][i]),2)-K[11][i]*
         pow((tau-K[12][i]),2))*(
         (-2*K[10][i]*pow(delta,K[1][i])+4*pow(K[10][i],2)*pow(delta,K[1][i])*
         pow((delta-K[13][i]),2))+
         (-4*K[1][i]*K[10][i]*pow(delta,(K[1][i]-1))*(delta-K[13][i])+
         K[1][i]*(K[1][i]-1)*pow(delta,(K[1][i]-2)))));

   }

   for (i=limit[2];i<limit[3];i++)
   {
      thetafn        = theta_fn(tau,K[6][i],delta,K[11][i]);
      deltafn        = delta_fn(thetafn,K[7][i],delta,K[4][i]);
      phifn          = phi_fn(K[8][i],delta,K[9][i],tau);
      dphiddelta      = dphi_ddelta(K[8][i],delta,phifn);
      d2phiddelta2   = d2phi_ddelta2(K[8][i],delta,phifn);
      dDELTA_deriv   = dDELTA_ddelta(delta,K[6][i],thetafn,K[11][i],K[7][i],K[4][i]);
      dDELTApowbddelta = dDELTApowb_ddelta(K[5][i],deltafn,dDELTA_deriv);
      dphiddelta     = dphi_ddelta(K[8][i],delta,phifn);
      d2DELTA_deriv = d2DELTA_ddelta2(delta,dDELTA_deriv,K[6][i],thetafn,K[11][i],K[7][i],K[4][i]);
      d2DELTApowbddelta2 = d2DELTApowb_ddelta2(K[5][i],deltafn,d2DELTA_deriv,dDELTA_deriv);

      phi_d = phi_d + (K[0][i]*(pow(deltafn,K[5][i])*(2*dphiddelta+delta*d2phiddelta2)
         +2*dDELTApowbddelta*(phifn+delta*dphiddelta)+d2DELTApowbddelta2*delta*phifn));

   }

   return phi_a+phi_b+phi_c+phi_d;
}


/**********************************************************************
Function for calculating the Pressure of a gas/liquid density on density
and temperature.

Parameters:
         rho  - density
         rhoc - density at the critical point
         T    - temperature
         Tc   - critical temperature
         R    - specific gas constant

Programming: NB
Aug 2008
***********************************************************************/
double pressure (double rho, double T, int c)
{
   CFluidProperties *a;
   double P,R,rhoc;
   a=MFPGet(c);
   rhoc=a->rhoc;
   R=a->Rs;

   P = (1+(rho/rhoc)*a->phi_r_d(rho,T,c))*rho*R*T;

   return P;
}


/**********************************************************************
This function calculates the density depending on pressure and temperature
by iteration. The iteration process may take a long time, so it's not re-
comended to use this function for simulations. Better use the GetMatrixValue
function and a table with rho-p-T data.
This function does not work for every pressure/temperature Range!

Parameters:
         P    - pressure
         rho0 - initial density for iteration
         rhoc - density at the critical point
T    - temperature
Tc   - critical temperature
R    - specific gas constant
prec - precision for iteration criteria

Programming: NB
Aug 2008
last change: NB 4.9.05
***********************************************************************/
double density (double P, double rho0, double T, double prec, int c)
{
   CFluidProperties *a;
   int iterations=0;
   double rho=0.0,p0,rhoc,R;                      //OK411

   a=MFPGet(c);
   rhoc=a->rhoc;
   R=a->Rs;

   p0 = 0;
   while (fabs(P-p0) > prec)                      //change to fabs. 15.09.2008 WW
   {
      rho = P/((1+(rho0/rhoc)*a->phi_r_d(rho0,T,c))*R*T);
      p0  = pressure(rho0,T,c);
      rho0 = rho;
      iterations++;
      if (iterations > 50) return 0;
   }

   return rho;
}


/**********************************************************************
Function for calculating enthalpy depending on density and
Temperature.

Parameters:
         rho  - density
         rhoc - density at the critical point
         T    - temperature
         Tc   - critical temperature
         R    - specific gas constant

Programming: NB
Aug 2008
***********************************************************************/
double enthalpy (double rho, double T, int c)
{
   CFluidProperties *a;
   double h,R;
   double tau, delta;

   a=MFPGet(c);

   tau = a->Tc/T;
   delta = rho/a->rhoc;
   R=a->Rs;

   h = (1+tau*(a->phi_0_t(T,c)+a->phi_r_t(rho,T,c))+delta*a->phi_r_d(rho,T,c))*R*T;

   return h;
}


/**********************************************************************
Function for calculating isochoric heat capacity depending on density and
Temperature.

Parameters:
         rho  - density
         rhoc - density at the critical point
         T    - temperature
         Tc   - critical temperature
         R    - specific gas constant

Programming: NB
Aug 2008
last change: NB 4.9.05
***********************************************************************/
double isochoric_heat_capacity (double rho, double T, int c)
{
   CFluidProperties *a;
   double cv;

   a=MFPGet(c);

   //	thermal_properties (fluid, rhoc, Tc, R);
   cv = -(pow((a->Tc/T),2)*(a->phi_0_tt(T,c)+a->phi_r_tt(rho,T,c)))*a->Rs;

   return cv;
}


/**********************************************************************
Function for calculating isobaric heat capacity depending on density and
Temperature.

Parameters:
         rho  - density
         rhoc - density at the critical point
         T    - temperature
         Tc   - critical temperature
         R    - specific gas constant

Programming: NB
Aug 2008
***********************************************************************/

double isobaric_heat_capacity (double rho, double T, int c)
{
   CFluidProperties *mfp_prop;
   double cp,delta,tau;
   //	double cp,delta,tau;
   mfp_prop = MFPGet(c);

   tau = mfp_prop->Tc/T;
   delta = rho/mfp_prop->rhoc;

   cp = (-pow(tau,2)*(mfp_prop->phi_0_tt(T,c)+mfp_prop->phi_r_tt(rho,T,c))
      +(pow((1+delta*mfp_prop->phi_r_d(rho,T,c)-delta*tau*mfp_prop->phi_r_dt(rho,T,c)),2))
      /((1+2*delta*mfp_prop->phi_r_d(rho,T,c)+pow(delta,2)*mfp_prop->phi_r_dd(rho,T,c))))*mfp_prop->Rs;

   return cp;
}


/**********************************************************************
Function for calculating viscosity of CO2 depending on density and Temperature.
Programming: NB
          Aug 2008
***********************************************************************/
double co2_viscosity (double rho, double T)
{
   double eta, eta_0, d_eta;
   double psi,t_r;
   double a[5],b[8],d[8][5];
   int i,j;

   // coefficients of the representation of the zero-density viscosity of co2
   a[0] =  0.235156;
   a[1] = -0.491266;
   a[2] =  0.05211155;
   a[3] =  0.05347906;
   a[4] = -0.01537102;

   psi = 0;                                       // variable for reduced effective cross-section

   // coefficients of the representation of the excess viscosity of co2

   for (i=0; i<8; i++)
   {
      for (j=0; j<5; j++)
      {
         b[i] = 0;
         d[i][j] = 0;
      }
   }

   d[0][0] = 0.4071119e-02;
   d[1][0] = 0.7198037e-04;
   d[5][3] = 0.2411697e-16;
   d[7][0] = 0.2971072e-22;
   d[7][1] =-0.1627888e-22;

   t_r = T/251.196;                               //reduced temperature

   // deriving the zero-density viscosity eta_0(T)

   for (i=0; i<5; i++)
   {
      psi = psi + a[i] * pow(log (t_r),(i));
   }
   psi = exp (psi);

   eta_0 = 1.00697 * pow (T,0.5) / psi;
   d_eta = 0;

   // deriving the excess viscosity d_eta(rho,T)

   for (i=0; i<8; i++)
   {
      for (j=0; j<5; j++)
      {
         b[i] = b[i] + d[i][j]/pow(t_r,j);
      }
      d_eta = d_eta + b[i] * pow(rho,i+1);
   }

   // deriving dynamic viscosity as sum of eta_o(T) and d_eta(rho,T)

   eta = (eta_0+d_eta)* 1e-06;                    // eta in [Pa s]

   return eta;
}


/**********************************************************************
Function for calculating heat conductivity of CO2 depending on density and Temperature.
   (Vesovic&Wakeham)
Programming: NB 4.8.01
          Nov 2008
***********************************************************************/
double co2_heat_conductivity (double rho, double T)
{
   double b[8],c[6],d[5];
   double G_fn=0,T_r,r,c_int_k,sum_c=0;
   int i;
   double lamda_0,delta_lamda=0,lamda;

   b[0]=0.4226159;
   b[1]=0.6280115;
   b[2]=-0.5387661;
   b[3]=0.6735941;
   b[4]=0;
   b[5]=0;
   b[6]=-0.4362677;
   b[7]=0.2255338;

   c[1]=0.02387869;
   c[2]=4.35079400;
   c[3]=-10.33404000;
   c[4]=7.98159000;
   c[5]=-1.94055800;

   d[1]=2.4471640E-02;
   d[2]=8.7056050E-05;
   d[3]=-6.5479500E-08;
   d[4]=6.5949190E-11;

   for (i=1;i<6;i++)
   {
      sum_c = sum_c + c[i]*pow((T/100),(2-i));
   }

   c_int_k = (1 + exp(-183.5/T))*sum_c;

   r = pow((2*c_int_k/5),(0.5));

   T_r = T/251.196;

   for (i=0;i<8;i++)
   {
      G_fn = G_fn + (b[i]/pow(T_r,i));
   }

   r =

      lamda_0 = (0.475598*pow(T,0.5)*(1+pow(r,2)))/G_fn;

   for (i=1;i<5;i++)
   {
      delta_lamda = delta_lamda + d[i]*pow(rho,i);
   }

   lamda = (lamda_0 + delta_lamda)/1000;

   return lamda;
}


/**********************************************************************
Function for calculating viscosity of CH4 at 295K depending on pressure.
   (Gulik,Mostert,Van den Berg)
Programming: NB 4.8.01
          Nov 2008
***********************************************************************/
double ch4_viscosity_295K (double p)
{
   double h;

   p=p/100000;
   h=(-3.7091411E-14*pow(p,4)+9.1937114E-10*pow(p,3)-6.6099446E-06*pow(p,2)+4.8400147E-02*p+1.0934694E+01)/1.0e6;

   return h;
}


/**********************************************************************
Function for calculating viscosity of pure water at a density rho and
a temperature T.
Programming: NB
          Apr 2009
***********************************************************************/
double h2o_viscosity_IAPWS (double rho, double T)
{
   double my,my_0,my_1;
   double H[4],h[6][7];
   double sum1=0,sum2=0,sum3=0;
   int i,j;

   T = T/647.096;
   rho=rho/322.0;

   H[0]=1.67752;
   H[1]=2.20462;
   H[2]=0.6366564;
   H[3]=-0.241605;

   for (i=0;i<6;i++) for (j=0;j<7;j++)  {h[i][j] = 0;}
   h[0][0]= 0.520094000;  h[1][0]= 0.085089500;  h[2][0]=-1.083740000;  h[3][0]=-0.289555000;  h[0][1]= 0.222531000;
   h[1][1]= 0.999115000;  h[2][1]= 1.887970000;  h[3][1]= 1.266130000;  h[5][1]= 0.120573000;  h[0][2]=-0.281378000;
   h[1][2]=-0.906851000;  h[2][2]=-0.772479000;  h[3][2]=-0.489837000;  h[4][2]=-0.257040000;  h[0][3]= 0.161913000;
   h[1][3]= 0.257399000;  h[0][4]=-0.032537200;  h[3][4]= 0.069845200;  h[4][5]= 0.008721020;  h[3][6]=-0.004356730;
   h[5][6]=-0.000593264;

   for(i=0;i<4;i++)
   {
      sum1=sum1+(H[i]/(pow(T,i)));
   }

   my_0 = 100*pow(T,0.5)/sum1;

   for(i=0;i<6;i++)
   {
      for (j=0;j<7;j++)
      {
         sum3 = sum3 + h[i][j]*pow(rho-1,j);
      }
      sum2= sum2 + (pow(1/T-1,i) * sum3);
      sum3 = 0;
   }

   my_1 = exp(rho*sum2);

   my = (my_0 * my_1)/1e6;
   return my;
}


/**********************************************************************
Viscosity for different Fluids

Programming: NB 4.8.01
          Nov 2008
***********************************************************************/
double Fluid_Viscosity (double rho, double T, double p, int fluid)
{
   p = p;                                         //OK411
   //TODO: make a global function for all properties: Fluid_Property(int c, int property, double T, double P) (NB)

   double h;

   switch (fluid)
   {
      case 0 :                                    // CARBON_DIOXIDE
         h = co2_viscosity (rho,T);
         break;
      case 1 :                                    // WATER
         // h = 1E-3;
         h = h2o_viscosity_IAPWS (rho,T);
         break;
      case 2 :                                    // METHANE
         //h = ch4_viscosity_295K(p);
         h = ch4_viscosity (rho,T);
         break;

      case 3 :                                    // Nitrogen
         h = n2_viscosity (rho,T);
         break;
      default :   h = 1E-3;
   }

   return h;
}


/**********************************************************************
Heat conductivity for different Fluids

Programming: NB 4.8.01
          Nov 2008
***********************************************************************/
double Fluid_Heat_Conductivity (double rho, double T, int fluid)
{
   double h;

   switch (fluid)
   {
      case 0 :                                    // CARBON_DIOXIDE
         h = co2_heat_conductivity (rho,T);
         break;
      case 1 :                                    // WATER
         // h = 0.598; // [W/m/K] at 293K and 1 bar
         h= h2o_heat_conductivity_IAPWS_ind(rho,T);
         break;
      case 2 :                                    // METHANE
         h = ch4_heat_conductivity (rho,T);
         break;
      case 3 :                                    // NITROGEN
         h= n2_heat_conductivity (rho,T);
         break;
      case 4 :                                    // NITROGEN
         h= n2_heat_conductivity (rho,T);
         break;

      default :   h = 0.5;
   }

   //if(caption.find("CARBON_DIOXIDE")!=string::npos)
   //
   //if(caption.find("METHANE")!=string::npos)
   //      h = 0.0338; // [W/m/K] at 298K and 1 bar
   //if(caption.find("WATER")!=string::npos)
   //      h = 0.598; // [W/m/K] at 293K and 1 bar

   return h;
}


/*************************************************************
 * vapour pressure for co2
 * NB, Dec 08
last change: NB 4.9.05
*************************************************************/
double vapour_pressure_co2(double T)
{
   double Tc=304.128;
   double pc=7377300;
   double p=0;
   double a[4],t[4];
   int i;
   a[0]=-7.0602087;  t[0]=1;
   a[1]=1.9391218;       t[1]=1.5;
   a[2]=-1.6463597;  t[2]=2;
   a[3]=-3.2995634;  t[3]=4;

   for (i=0;i<4;i++)
      p=p+a[i]*pow(1-(T/Tc),t[i]);

   p=exp(Tc/T*p)*pc;

   return p;
}


/*************************************************************
 * sublime pressure for co2
 * NB, Dec 08
last change: NB 4.9.05
*************************************************************/
double sublime_pressure_co2(double T,double Tt,double pt)
{
   double p;
   double a[3];

   a[0]=-14.740846;
   a[1]=2.4327015;
   a[2]=-5.3061778;

   p= exp(Tt/T*(a[0]*(1-T/Tt)+a[1]*pow((1-T/Tt),1.9)+a[2]*pow((1-T/Tt),2.9)))*pt;

   return p;
}


/*************************************************************
 * melting pressure for co2
 * NB, Dec 08
last change: NB 4.9.05
*************************************************************/
                                                  //just for CO2
double melting_pressure_co2(double T,double Tt,double pt)
{
   double p;
   double a[2];
   a[0]=1955.539;                                 //CO2
   a[1]=2055.4593;

   p= (1+a[0]*(T/Tt-1)+a[1]*pow((T/Tt-1),2))*pt;

   return p;
}


/*************************************************************
 * Peng&Robinson Equation of State
 * Analytical solving of third grade polynomial
 *
 * Parameters: temperature and pressure
 *             caption defines fluid
 * Programming: NB, Dec 08
 **************************************************************/
double preos(double T, double P, int c)
{
   double z1,z2,z3,h;
   vector<double> roots;

   CFluidProperties *mfp_prop;

   double a,b,Tc,pc,MM,Ru;
   double omega,alpha;

   //int i;
   mfp_prop = MFPGet (c);

   Ru=mfp_prop->Ru;                               //universal gas constant
   MM=mfp_prop->molar_mass;
   Tc=mfp_prop->Tc;                               // critical temperature
   pc=mfp_prop->pc/1000;                          //critical pressure
   omega=mfp_prop->omega;                         // azentric factor

   // Peng Robinson EOS:
   // P= R*T / (V-b) - a*alpha / (V^2+2bV-b^2)   where V = MM/rho

   a=0.457235*pow(Ru,2)*pow(Tc,2)/pc;
   b=0.077796*Ru*Tc/pc;
   P=P/1000;                                      //P in kPa
   alpha=pow((1+(0.37464+1.5422*omega-0.26992*pow(omega,2))*(1-pow((T/Tc),0.5))),2);

   //EOS in the form: 0 = rho^3 + z1*rho^2 + z2*rho + z3

   z1=(MM*a*alpha-3*MM*pow(b,2)*P-2*MM*Ru*T*b)/(b*(P*pow(b,2)+b*Ru*T-a*alpha));
   z2=(pow(MM,2)*(b*P-Ru*T))/(b*(P*pow(b,2)+b*Ru*T-a*alpha));
   z3=(pow(MM,3)*P)/(b*(P*pow(b,2)+b*Ru*T-a*alpha));

   NsPol3(z1,z2,z3,&roots);                       //derives the roots of the polynomial

   h=FindMin(roots);
   return h;                                      //returns the lowest positive root
}


/////****************************************************************************
////* Finds and returns the positive minimum of a vector.
////* Programming: NB Dec 08
////*****************************************************************************/
////double FindMin (vector<double>Vec)
////{
////double x=DBL_MAX;
////int unsigned i;
////
////for(i=0;i<Vec.size();i++) {if ((Vec[i]>=0)&&(Vec[i]<x)) x=Vec[i];}
////
////return x;
////
////}

/*************************************************************
 * Redlich&Kwong Equation of State
 * Analytical solving of third grade polynomial
 *
 * Parameters: temperature and pressure
 *             caption defines fluid
 * Programming: NB, Jan 09
 **************************************************************/
double rkeos(double T, double P, int c)
{
   double z1,z2,z3,h;
   vector<double> roots;

   CFluidProperties *mfp_prop;

   double a,b,Tc,pc,MM;
   double Ru;

   if (P<0) P=100000;                             // set pressure to 1atm if unstable NB
   mfp_prop = MFPGet (c);

   Ru=mfp_prop->Ru*10;                            //universal gas constant [bar cm3/mol/K]
   MM=mfp_prop->molar_mass;
   Tc=mfp_prop->Tc;                               // critical temperature
   pc=mfp_prop->pc/100000;                        //critical pressure

   // Redlich-Kwong EOS:
   // P= R*T (1+y+y^2-y^3)/ v(1-y^3) - a / (T^0.5*v(cv+b)   where V = MM/rho and y = b / (4v)

   a=27*pow(Ru,2)*pow(Tc,2.5)/(64*pc);
   b=0.0866*Ru*Tc/pc;
   P=P/100000;                                    //P in bar

   //EOS in the form: 0 = vm^3 + z1*vm^2 + z2*vm + z3

   z1=-Ru*T/P;
   z2=-(Ru*T*b/P-a/(pow(T,0.5)*P)+pow(b,2));
   z3=-a*b/(pow(T,0.5)*P);

   NsPol3(z1,z2,z3,&roots);                       //derives the roots of the polynomial

   h=FindMax(roots);                              //returns the lowest positive root (molar volume)
   h=MM/h*1000;                                   // density in kg/m3
   return h;
}


/*************************************************************
 * Redlich&Kwong Equation of State
 * Analytical solving of third grade polynomial
 *
 * This form of RKEOS requires adjusting parameters a and b as
input values, also the molecular mass

* Parameters: temperature and pressure
*             adj. parameters a and b
*             molecular mass

* Programming: NB, May 09
**************************************************************/
double rkeos(double T, double P, double MM, double a, double b)
{
   double z1,z2,z3,h;
   vector<double> roots;

   double Ru;
   P=P/100000;                                    //P in bar
   Ru=83.14472;                                   //universal gas constant [bar cm3/mol/K]
   // Redlich-Kwong EOS:
   // 0 = vm^3 + z1*vm^2 + z2*vm + z3

   z1=-Ru*T/P;
   z2=-(Ru*T*b/P-a/(pow(T,0.5)*P)+pow(b,2));
   z3=-a*b/(pow(T,0.5)*P);

   NsPol3(z1,z2,z3,&roots);                       //derives the roots of the polynomial

   h=FindMax(roots);                              //returns the lowest positive root (molar volume)
   h=MM/h*1000;                                   // density in kg/m3
   return h;
}


/**********************************************************************
Function for calculating thermal conductivity of pure water at a density
rho and a temperature T.
IAPWS Formulation 1997 for industrial use
Programming: NB
          Apr 2009
***********************************************************************/
double h2o_heat_conductivity_IAPWS_ind (double rho, double T)
{
   double lamda,lamda_0,lamda_1,lamda_2;
   double sum1=0;
   double S,Q,dT;
   double a[4],b[3],B[2],d[4],C[6];
   int i;

   T = T/647.26;
   rho = rho/317.7;

   a[0] =  0.0102811;
   a[1] =  0.0299621;
   a[2] =  0.0156146;
   a[3] = -0.00422464;

   b[0] = -0.397070;
   b[1] =  0.400302;
   b[2] =  1.060000;

   B[0] = -0.171587;
   B[1] =  2.392190;

   d[0] = 0.0701309;
   d[1] = 0.0118520;
   d[2] = 0.00169937;
   d[3] = -1.0200;

   C[0] = 0.642857;
   C[1] = -4.11717;
   C[2] = -6.17937;
   C[3] = 0.00308976;
   C[4] = 0.0822994;
   C[5] = 10.0932;

   for (i=0;i<4;i++)
   {
      sum1=sum1 + a[i]*pow(T,i);
   }

   lamda_0 = pow(T,0.5)*sum1;
   lamda_1 = b[0] + b[1]*rho + b[2]* exp(B[0]*pow(rho+B[1],2));

   dT = fabs(T-1)+C[3];
   Q = 2 + (C[4]/pow(dT,3./5.));

   if (T>=1)
      S = 1/dT;
   else
      S = C[5]/pow(dT,3./5.);

   lamda_2 = (d[0]/pow(T,10.)+d[1])*pow(rho,9./5.)*exp(C[0]*(1-pow(rho,14./5.)))
      + d[2]*S*pow(rho,Q)*exp((Q/(1.+Q))*(1-pow(rho,(1.+Q))))
      + d[3]*exp(C[1]*pow(T,3./2.)+C[2]/pow(rho,5.));

   lamda = (lamda_0 + lamda_1 + lamda_2);         // lamda in [W/m/K]

   return lamda;
}


/**********************************************************************
Function for calculating viscosity of methane (CH4) depending on density
and Temperature.

see
Friend, Ely, Ingham: The Transport Properties of Methane,
J. Chem. Phys. Ref. Data,1989.

Programming: NB
          Apr 2009
***********************************************************************/
double ch4_viscosity (double rho, double T)
{
   double eta, eta_0, eta_ex;
   double t,delta,tau,Omega=0;
   double C[9],r[11],s[11],g[11];
   double sum1=0,sum2=0;
   int unsigned i;

   rho=rho/16.043;                                //rho in [mol/dm3]
   t = T/174.;

   C[0] =  -3.0328138281;     C[1] =  16.918880086;    C[2] = -37.189364917;     C[3] =  41.288861858;    C[4] = -24.615921140;
   C[5] =   8.9488430959;     C[6] =  -1.8739240542;   C[7] =   0.20966101390;   C[8] =  -9.6570437074e-03;

   r[0] = 1;  r[1] = 1;  r[2] = 2; r[3] = 2;  r[4] = 2;   r[5] = 3;  r[6] = 3; r[7] = 4; r[8] = 4;   r[9] = 1;  r[10] = 1;
   s[0] = 0;  s[1] = 1;  s[2] = 0; s[3] = 1;  s[4] = 1.5; s[5] = 0;  s[6] = 2; s[7] = 0; s[8] = 1;   s[9] = 0;  s[10] = 1;

   g[0] = 0.41250137;  g[1] =-0.14390912;  g[2] = 0.10366993;  g[3] = 0.40287464;  g[4] = -0.24903524;  g[5] = -0.12953131;
   g[6] = 0.06575776;  g[7] = 0.02566628;  g[8] = -0.03716526; g[9] = -0.38798341; g[10] = 0.03533815;

   for (i=0;i<9;i++)
   {
      Omega = Omega + C[i]*pow(t, (i/3.-1));
   }
   Omega = 1/Omega;

   eta_0 = 10.5*pow(t,0.5)/Omega;

   delta = rho/10.139;
   tau = 190.551/T;

   for (i=0;i<9;i++)
   {
      sum1 = sum1 + g[i]*pow(delta,r[i])*pow(tau,s[i]);
   }

   for (i=9;i<11;i++)
   {
      sum2 = sum2 + g[i]*pow(delta,r[i])*pow(tau,s[i]);
   }

   eta_ex = 12.149 * sum1 * pow(1+sum2,-1.0);

   eta = (eta_0 + eta_ex)*1e-6;
   return eta;
}


/**********************************************************************
Function for calculating viscosity of methane (CH4) depending on density
and Temperature.

see
Friend, Ely, Ingham: The Transport Properties of Methane,
J. Chem. Phys. Ref. Data,1989.

Programming: NB
          Apr 2009
***********************************************************************/
double ch4_heat_conductivity (double rho, double T)
{
   double lamda, eta_0, lamda_0, lamda_ex;
   double t,delta,tau,Omega=0,epsilon,beta;
   double C[9],Q[6],r[7],s[7],j[7],f[2],H[6],J[5];
   double sum=0,phi_id_tt,f_int;
   double P_sigma,Z_c,delta_star,T_star,rho_sigma_v;

   int unsigned i;
   const double T_c = 190.551;
   const double rho_c = 10.139;
   const double P_c = 4.599200;
   const double R = 8.314510;                     // [J/mol/K]

   t = T/174.;
   tau = T_c/T;
   T_star = (T_c-T)/T_c; if (T_star < 0) T_star = 0;
   rho=rho/16.043;                                //rho in [mol/dm3]
   delta = rho/rho_c;

   C[0] =  -3.0328138281;     C[1] =  16.918880086;    C[2] = -37.189364917;     C[3] =  41.288861858;    C[4] = -24.615921140;
   C[5] =   8.9488430959;     C[6] =  -1.8739240542;   C[7] =   0.20966101390;   C[8] =  -9.6570437074e-03;
   Q[0] = 2.5998324; Q[1] = -3.3854083; Q[2] = 1.6900979; Q[3] = -0.3911541; Q[4] = 4.7206715; Q[5] = -10.543907;

   f[0] = 1.458850; f[1] = -0.4377162;
   r[0] = 1; r[1] = 3; r[2] = 4; r[3] = 4; r[4] = 5; r[5] = 5; r[6] = 2;
   s[0] = 0; s[1] = 0; s[2] = 0; s[3] = 1; s[4] = 0; s[5] = 1; s[6] = 0;

   j[0] = 2.414927;   j[1] = 0.55166331;   j[2] = -0.52837734;   j[3] = 0.073809553;   j[4] = 0.24465507;   j[5] = -0.047613626;  j[6] = 1.5554612;
   H[0] = -6.589879;  H[1] = 0.6355175;  H[2] = 11.31028;  H[3] = -10.38720;  H[4] = 3.393075;
   J[0] = -0.7377483; J[1] = -1.241532;  J[2] = -1.649972; J[3] = 2.281949;  J[4] = 1.439570;

   beta = 0.355;
   epsilon = 1.9;

   for (i=0;i<9;i++)
   {
      Omega = Omega + C[i]*pow(t, (i/3.-1));
   }
   Omega = 1/Omega;

   eta_0 = 10.5*pow(t,0.5)/Omega;                 // (Eq. 10a)

   phi_id_tt = - Q[0]+ 4*Q[1]/9  * pow(tau,-1./3.) + 10*Q[2]/9 * pow(tau,-2./3.) + 2*Q[3]*pow(tau,-1)
      - Q[4]*pow(Q[5],2)* pow(tau,2)* exp (Q[5]*tau) * pow((exp (Q[5]*tau)-1),-2);

   f_int = f[0] + (f[1]/t);

                                                  //(Eq. 14)
   lamda_0 = 0.51826 * eta_0 * (3.75-f_int*(phi_id_tt+1.5));

                                                  //(Eq. 3)
   P_sigma = P_c * exp(H[0]*T_star/(1-T_star)+H[1]*T_star+H[2]*pow(T_star,epsilon)+H[3]*pow(T_star,2)+H[4]*pow(T_star,3));
   Z_c = P_c/(R*T_c*rho_c);

   rho_sigma_v = P_sigma/(R*T) * (1+P_sigma*pow(tau,8)*(Z_c-1)/P_c*(1+(J[0]*pow(T_star,beta)+J[1]*pow(T_star,2*beta)+J[2]*(T_star+pow(T_star,4))+J[3]*pow(T_star,2))
      /(1+J[4]*T_star)));                         // (Eq. 5)

   if ((T<T_c)&&(rho<rho_c))                      // (Eq. 16)
      delta_star = rho_sigma_v/rho_c;
   else
      delta_star = 11;

   for (i=0;i<7;i++)                              // (Eq. 17)
   {
      sum = sum + j[i]*pow(delta,r[i])*pow(tau,s[i]);
   }

                                                  // (Eq. 17)
   lamda_ex = 6.29638 * ( sum + j[6]*pow(delta,2)/delta_star);

   lamda = (lamda_0 + lamda_ex)/1000;             //lamda in [W/m/K]

   return lamda;
}


/**********************************************************************
Function for calculating viscosity of nitrogen (N2) depending on density
and Temperature.

see
Stephan, Krauss and Laeseke: Viscosity and themal conductivity of fluid
Nitrogen, J. Chem. Phys. Ref. Data,Vol. 16, No. 4, 1987.

Programming: NB
          Apr 2009
***********************************************************************/
double n2_viscosity (double rho, double T)
{
   //OK411 const double M = 28.013;
   //OK411 const double T_c = 126.2; //[K]
   //OK411 const double P_c = 3.4; // [MPa]
   const double rho_c = 314;                      // [kg/m3]
   //OK411 const double T_t = 63.1; //[K]
   const double CVF = 14.058;                     // [1e-3 Pa-s]

   const double sigma = 0.36502496e-09;
   const double k = 1.38062e-23;
   const double eps = 138.08483e-23;
   //OK411 const double N_A = 6.02213E26;
   //OK411 const double pi = 3.14159;
   const double c1 = 0.3125;
   const double c2 = 2.0442e-49;

   double A[5],C[5];

   double sum=0,eta, eta_0,eta_r,T_star,Omega=0;
   int i;

   A[0] = 0.46649;
   A[1] = -0.57015;
   A[2] = 0.19164;
   A[3] = -0.03708;
   A[4] = 0.00241;

   C[0] = -20.09997;
   C[1] = 3.4376416;
   C[2] = -1.4470051;
   C[3] = -0.027766561;
   C[4] = -0.21662362;

   T_star = T*k/eps;
   rho = rho/rho_c;

   for (i=0;i<5;i++)
   {
      Omega = Omega + A[i]*pow(log(T_star),i);
   }

   Omega = exp (Omega);

                                                  //eta in [Pa*s]
   eta_0 = c1 * pow(c2*T,0.5) / (pow(sigma,2)*Omega);

   for (i=2;i<5;i++)
   {
      sum = sum + C[i]*pow(rho,i-1);
   }

                                                  //
   eta_r= CVF *1e-6*(C[0]/(rho-C[1]) + C[0]/C[1] + sum);

   eta =  (eta_0 + eta_r);                        // [Pa*s]

   return eta;
}


/**********************************************************************
Function for calculating thermal conductivity of nitrogen (N2) depending
on density and Temperature.

see
Stephan, Krauss and Laeseke: Viscosity and themal conductivity of fluid
Nitrogen, J. Chem. Phys. Ref. Data,Vol. 16, No. 4, 1987.

Programming: NB
          Apr 2009
***********************************************************************/
double n2_heat_conductivity (double rho, double T)
{
   const double X1 = 0.95185202;
   const double X2 = 1.0205422;

   const double rho_c = 314;                      // [kg/m3]
   const double M = 28.013;
   const double k = 1.38062e-23;
   const double eps = 138.08483e-23;
   const double N_A = 6.02213E26;
   const double R = 8.31434;
   const double CCF = 4.173;                      //mW/m/K

   const double c1 = 0.3125;
   const double c2 = 2.0442e-49;
   const double sigma = 0.36502496e-09;

   double F;
   double A[5],f[9],C[4];
   double sum=0,eta_0,c_v0,T_star,Omega=0;
   double lamda_tr,lamda_in,lamda_r,lamda_0,lamda;

   int i;

   T_star = T*k/eps;
   rho = rho/rho_c;

   A[0] = 0.46649;
   A[1] = -0.57015;
   A[2] = 0.19164;
   A[3] = -0.03708;
   A[4] = 0.00241;

   f[0] =-0.837079888737e3;
   f[1] = 0.37914711487e2;
   f[2] =-0.601737844275;
   f[3] = 0.350418363823e1;
   f[4] =-0.874955653028e-5;
   f[5] = 0.148968607239e-7;
   f[6] =-0.256370354277e-11;
   f[7] = 0.100773735767e1;
   f[8] = 0.335340610e4;

   C[0] =3.3373542;
   C[1] =0.37098251;
   C[2] =0.89913456;
   C[3] =0.16972505;

   // dilute heat conductivity
   for (i=0;i<7;i++)
   {
      sum = sum + f[i]*pow(T,(i-3));
   }
   c_v0 = R * (sum + ((f[7]*pow((f[8]/T),2)*(exp((f[8]/T))))/(pow(exp ((f[8]/T)) -1,2))-1));
   sum = 0;

   double cvint;
   cvint = c_v0*1000/N_A;

   // dilute gas viscosity
   for (i=0;i<5;i++)
   {
      Omega = Omega + A[i]*pow(log(T_star),i);
   }
   Omega = exp (Omega);

                                                  //eta in [Pa*s]
   eta_0 = 1e6*(c1 * pow(c2*T,0.5) / (pow(sigma,2)*Omega));

   F = eta_0 * k * N_A / (M*1000);

   lamda_tr = 2.5 * (1.5-X1);
   lamda_in = X2*(cvint/k+X1);

   lamda_0 = F * (lamda_tr + lamda_in);
   sum = 0;
   for (i=0;i<4;i++)
   {
      sum = sum + C[i]*pow(rho,(i+1));
   }

   lamda_r = sum*CCF;

   lamda = (lamda_0 + lamda_r)/1000;              //lamda in [W/m/K]

   return lamda;
}


/*************************************************************
saturation vapour pressure at a temperature.

Wagner, W. and Kretzschmar, H.-J., 1998: International Steam
Tables. Properties ofWater and Steam Based on the Industrial
Formulation IAPWS-IF97. Second edition. Springer-Verlag.
P. 25

Input: T in K
Output: p in bar

Programming: Dedong Li

*************************************************************/
double vapour_pressure_IF97(double T)
{
   double A,B,C,theta,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10;
   n1 =  0.11670521452767e+04;   n2 = -0.72421316703206e+06;   n3 = -0.17073846940092e+02;
   n4 =  0.12020824702470e+05;   n5 = -0.32325550322333e+07;   n6 =  0.14915108613530e+02;
   n7 = -0.48232657361591e+04;   n8 =  0.40511340542057e+06;   n9 = -0.23855557567849e+00;
   n10=  0.65017534844798e+03;

   theta = T+n9/(T-n10);

   A =    pow(theta,2.0)+n1*theta+n2;
   B = n3*pow(theta,2.0)+n4*theta+n5;
   C = n6*pow(theta,2.0)+n7*theta+n8;
   double h;
   h=pow((2.0*C/(-B+pow((B*B-4.0*A*C),0.5))),4.0)*10.0;
   return h;
}


/*************************************************************
Methane saturation vapour pressure at a temperature, Setzmann,1991.

Input: T in K
Output: p in Pa

Programming: NB
*************************************************************/
double vapour_pressure_ch4(double T)
{
   const double Tc=190.564;
   const double pc=4592200;

   double theta;
   double n1,n2,n3,n4,h;

   n1=-6.036219; n2= 1.409359; n3=-0.4945199; n4=-1.443048;

   theta = (1-T/Tc);
   h=Tc/T*(n1*theta+n2*pow(theta,1.5)+n3*pow(theta,2)+n4*pow(theta,4.5));

   return exp(h)*pc;;
}


/*************************************************************
Water saturation vapour pressure at a temperature, Wagner,2002.

Input: T in K
Output: p in Pa

Programming: NB
*************************************************************/
double vapour_pressure_h2o(double T)
{
   double theta;
   const double Tc=647.096;
   const double pc=22064000;

   double a1,a2,a3,a4,a5,a6,h;
   a1= -7.85951783; a2= 1.84408259; a3=-11.7866497; a4=22.6807411; a5=-15.9618719; a6=1.80122502;

   theta = (1-T/Tc);

   h=Tc/T*(a1*theta+a2*pow(theta,1.5)+a3*pow(theta,3)+a4*pow(theta,3.5)+a5*pow(theta,4)+a6*pow(theta,7));

   return exp(h)*pc;;
}


/*************************************************************
returns min and max interpolation boundaries for zbrent

Programming: Dedong Li
             May 2009
Last modification: NB Jun 2009
*************************************************************/
double dpressure(double TT, double PP,int fluid, double ds)
{

   //return PP-pressure(ds, TT, fluid)*1.0E-5;
   return PP-pressure(ds, TT, fluid);
}


/*************************************************************
Programming: Dedong Li
             May 2009
Last modification:
*************************************************************/
inline float SIGN(const double &a, const float &b)
{
   return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}


/*************************************************************
returns min and max interpolation boundaries for zbrent

Programming: Dedong Li
             May 2009
Last modification: NB Jun 2009
*************************************************************/
void auswahl (double T, double P, int fluid, double* x1, double* x2, double eps)
{
   eps = eps;                                     //OK411
   switch (fluid)
   {
      case 0 :                                    //CO2
         *x1 = 0.0005;
         *x2 = 2000.0;

         if (T<304.128)                           //gas or liquid
         {
            if (P<(vapour_pressure_co2(T)))       // gas
            {
               *x1 = 5e-4;                        //min_density;
               *x2 = vapour_saturation_density_co2(T);
            }
            else                                  // liquid
            {
               *x1 = liquid_saturation_density_co2(T);
               *x2 = 2000;                        //max_density;
            }
         }

         break;
      case 1 :                                    //Water
         if(T<647.096 && P>(vapour_pressure_h2o(T)))
         {
            *x2=1500.0;
            if(T<=440.0)                *x1=9.0e+2;
            if(T> 440.0 && T<573.15)    *x1=7.1e+2;
            if(T>=573.15)               *x1=3.0e+2;
         }
         else
         {
            *x1=1.0e-8;
            *x2=2000;
            if(T<423.0) *x2=2.6;
            if(T<400.0) *x2=1.4;
            if(T<335.0) *x2=0.15;
         }

         break;

      case 2 :                                    //Methane
         *x1=0.0005;                              //min_density
         *x2=600;                                 //max_density

         if (T<190.564)                           //gas or liquid
         {
            if (P<(vapour_pressure_ch4(T)))       // gas
            {
               *x1 = 5e-4;                        //min_density;
               *x2 = vapour_saturation_density_ch4(T);
            } else                                // liquid
            {
               *x1 = liquid_saturation_density_ch4(T);
               *x2 = 600;                         //max_density;
            }
         }
         break;
      case 3 :                                    //Nitrogen

         {
            *x1=0.005;                            //min_density
            *x2=2000;                             //max_density

            if (T<126.192)                        //gas or liquid
            {
               if (P<(vapour_pressure_n2(T)))     // gas
               {
                  *x1 = 5e-4;                     //min_density;
                  *x2 = vapour_saturation_density_n2(T);
               } else                             // liquid
               {
                  *x1 = liquid_saturation_density_n2(T);
                  *x2 = 2000;                     //max_density;
               }
            }
         }

         break;
   }
}


/*************************************************************
Interpolates the Helmholts free energy pressure-temperature-
density relation and returns a density for a temperature and
a pressure

Programming: Dedong Li
             May 2009
Last modification: NB Jun 2009
*************************************************************/
double zbrent(double TT, double PP, int fluid, const double tol)
{
   const int ITMAX=100;
   const double EPS=5.0e-16;                      //numeric_limits<double>::epsilon();
   //double fa=func(a),fb=func(b);
   double fa,fb;
   double fc,p,q,r,s,tol1,xm;
   double x1;
   double x2;
   //	PP=PP/1e5; // P in bar
   auswahl(TT,PP,fluid,&x1,&x2,tol);
   double a=x1,b=x2,c=x2,d=0.0,e=0.0;             //OK411
   fa=dpressure(TT,PP,fluid,x1);
   fb=dpressure(TT,PP,fluid,x2);

   if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
   {                                              //cout << "Error in zbrent, fluid " << fluid << " T: " << TT << " P: " << PP << " b: " << b << endl;
      cout << ".";
   }
   fc=fb;
   for (int iter=0;iter<ITMAX;iter++)
   {
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
      {
         c=a;
         fc=fa;
         e=d=b-a;
      }
      if (fabs(fc) < fabs(fb))
      {
         a=b;
         b=c;
         c=a;
         fa=fb;
         fb=fc;
         fc=fa;
      }
      tol1=2.0*EPS*fabs(b)+0.5*tol;
      xm=0.5*(c-b);
      if (fabs(xm) <= tol1 || fb == 0.0)
      {
         return b;
      }
      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
      {
         s=fb/fa;
         if (a == c)
         {
            p=2.0*xm*s;
            q=1.0-s;
         }
         else
         {
            q=fa/fc;
            r=fb/fc;
            p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
            q=(q-1.0)*(r-1.0)*(s-1.0);
         }
         if (p > 0.0) q = -q;
         p=fabs(p);
         double min1=3.0*xm*q-fabs(tol1*q);
         double min2=fabs(e*q);
         if (2.0*p < (min1 < min2 ? min1 : min2))
         {
            e=d;
            d=p/q;
         }
         else
         {
            d=xm;
            e=d;
         }
      }
      else
      {
         d=xm;
         e=d;
      }
      a=b;
      fa=fb;
      if (fabs(d) > tol1)
         b += d;
      else
         b += (double)SIGN(tol1,(float)xm);       //OK411
      fb=dpressure(TT,PP,fluid,b);
   }
   throw("Maximum number of iterations exceeded in zbrent");
}


/**********************************************************************
Function for the mixing of adjusting parameters for cubic equations of state.
For ternary mixtures.

Fluid List:  0 - carbon dioxide
             1 - water
             2 - methane
             3 - nitrogen
             ...

Arguments: x[i] ... mole fraction of component i
a[i],b[i] ... parameters for pure substance
MM[i]... molar mass of pure substance
ra,rb,MM   ... mixed parameters

Programming: NB
May 2009
***********************************************************************/
double mixing_ternary (double* x,double* a,double* b,double *MM, double *ra, double *rb, double *rMM)
{
   int unsigned i,j,k;
   double a_mix=0,b_mix=0,MM_mix=0;

   for (i=0;i<3;i++)
      for (j=0;j<3;j++)
         for (k=0;k<3;k++)
            a_mix+=x[i]*x[j]*x[k]*pow(a[i]*a[j]*a[k],(1./3.));

   for (i=0;i<3;i++)
      b_mix+=x[i]*b[i];

   for (i=0;i<3;i++)
      MM_mix+=x[i]*MM[i];

   *ra=a_mix;
   *rb=b_mix;
   *rMM=MM_mix;

   return 0;
}


/*************************************************************
Methane liquid saturation density at a temperature, Setzmann,1991.

Input: T in K
Output: p in Pa

Programming: NB
*************************************************************/
double liquid_saturation_density_ch4(double T)
{
   const double Tc=190.564;
   const double rhoc=162.66;

   double theta;
   double n1,n2,n3,h;

   n1=1.9906389; n2= -0.78756197; n3=0.036976723;

   theta = (1-T/Tc);
   h=(n1*pow(theta,0.354)+n2*pow(theta,0.5)+n3*pow(theta,(5./2.)));

   return exp(h)*rhoc;;
}


/*************************************************************
Methane vapour saturation density at a temperature, Setzmann,1991.

Input: T in K
Output: p in Pa

Programming: NB
*************************************************************/
double vapour_saturation_density_ch4(double T)
{
   const double Tc=190.564;
   const double rhoc=162.66;

   double theta;
   double n1,n2,n3,n4,n5,n6,h;

   n1=-1.880284; n2=-2.8526531; n3=-3.000648; n4=-5.251169; n5=-13.191859; n6=-37.553961;

   theta = (1-T/Tc);
   h=(n1*pow(theta,0.354)+n2*pow(theta,(5./6.))+n3*pow(theta,(3./2.))+n4*pow(theta,(5./2.))+n5*pow(theta,(25./6.))+n6*pow(theta,(47./6.)));

   return exp(h)*rhoc;;
}


/*************************************************************
Nitrogen vapour pressure at a temperature, Setzmann,1991.

Input: T in K
Output: p in Pa

Programming: NB
*************************************************************/
double vapour_pressure_n2(double T)
{
   const double Tc=126.192;
   const double pc=3395800.;

   double theta;
   double n1,n2,n3,n4,h;

   n1=-6.12445284; n2=1.26327220; n3=-0.765910082; n4=-1.77570564;

   theta = (1-T/Tc);
   h=(Tc/T)*(n1*theta+n2*pow(theta,(1.5))+n3*pow(theta,(2.5))+n4*pow(theta,(5.)));

   return exp(h)*pc;
}


/*************************************************************
Nitrogen liquid saturation density at a temperature, Setzmann,1991.

Input: T in K
Output: density in kg/m3

Programming: NB
*************************************************************/
double liquid_saturation_density_n2(double T)
{
   const double Tc=126.192;
   const double rhoc=313.3;

   double theta;
   double n1,n2,n3,n4,h;                          //OK411 n5,n6,

   n1=1.48654237; n2=-0.280476066; n3=0.0894143085; n4=-0.119879866;

   theta = (1-T/Tc);
   h=(n1*pow(theta,0.3294)+n2*pow(theta,(2./3.))+n3*pow(theta,(8./3.))+n4*pow(theta,(35./6.)));

   return exp(h)*rhoc;;
}


/*************************************************************
Nitrogen vapour saturation density at a temperature, Setzmann,1991.

Input: T in K
Output: density in kg/m3

Programming: NB
*************************************************************/
double vapour_saturation_density_n2(double T)
{
   const double Tc=126.192;
   const double rhoc=313.3;

   double theta;
   double n1,n2,n3,n4,n5,h;                       //OK411 n6,

   n1=-1.70127164; n2=-3.70402649; n3=1.29859383; n4=-0.561424977; n5=-2.68505381;

   theta = (1-T/Tc);
   h=(Tc/T)*(n1*pow(theta,0.34)+n2*pow(theta,(5./6.))+n3*pow(theta,(7./6.))+n4*pow(theta,(13./6.))+n5*pow(theta,(14./3.)));

   return exp(h)*rhoc;
}


/*************************************************************
Carbon dioxide vapour saturation density at a temperature, Span,1996.

Input: T in K
Output: density in kg/m3

Programming: NB
*************************************************************/
double vapour_saturation_density_co2(double T)
{
   const double Tc=304.128;
   const double rhoc=467.6;

   double n[4],t[4],h=0;
   int i;

   n[0]=1.9245108;  t[0]=0.34;
   n[1]=-0.62385555;t[1]=0.5;
   n[2]=-0.32731127;t[2]=(10./6.);
   n[3]=0.39245142; t[3]=(11./6.);

   for (i=0;i<4;i++) h+=n[i]*pow((1-T/Tc),t[i]);

   return exp(h)*rhoc;
}


/*************************************************************
Carbon dioxide liqiud saturation density at a temperature, Span,1996.

Input: T in K
Output: density in kg/m

Programming: NB
*************************************************************/
double liquid_saturation_density_co2(double T)
{
   const double Tc=304.128;
   const double rhoc=467.6;

   double n[5],t[5],h=0;
   int i;

   n[0]=-1.7074879; t[0]=0.34;
   n[1]=-0.8227467; t[1]=0.5;
   n[2]=-4.6008549; t[2]=1;
   n[3]=-10.111178; t[3]=(7./3.);
   n[4]=-29.742252; t[4]=(14./3.);

   for (i=0;i<4;i++) h+=n[i]*pow((1-T/Tc),t[i]);

   return exp(h)*rhoc;
}


/**********************************************************************
Function thermal_properties (fluid, critical_density, critical_temperature, specific_gas_constant)
returns the thermal properties of a given fluid
Programming: NB Mar09
**********************************************************************/
void CFluidProperties::therm_prop (string caption)
{
   // CFluidProperties fpc;
   char letter;
   int i,j;
   //TODO: Change first letter approach to an integer (NB) - done (fluid_number) NB, 2009-06-02
   letter = caption[0];
   Ru=8.314472;
   switch (letter)
   {
      case 'C' :                                  // CARBON DIOXIDE
      {
         fluid_id = 0;
         rhoc=467.6;                              // critical density [kg/m3]
         Tc=304.1282;                             // critical temperature [K]
         pc=7377300;                              // critical pressure [Pa]
         Tt=216.592;                              // triple point temperature [K]
         pt=517950;                               //  triple point pressure [Pa]
         Rs=188.9241;                             // specific gas constant [J/kg/K]
         molar_mass=44.0099;                      // [g/mol]
         omega=0.22491;                           // azentric factor, see PREOS

         // Limits sums in FHE-derivations

         limit[0]=7;
         limit[1]=34;
         limit[2]=39;
         limit[3]=42;

         // Coefficients for FHE-derivations

         for (i=0;i<2;i++)
            for (j=0;j<8;j++)
               {k[i][j] = 0;}

               for (i=0;i<14;i++)
                  for (j=0;j<56;j++)
                     {K[i][j] = 0;}

                     // ideal gas part
                                                  //a
                     k[0][0]=  8.37304456;   k[0][1]=  -3.70454304;  k[0][2]=  2.5;       k[0][3]=  1.99427042;
         k[0][4]=  0.62105248;   k[0][5]=  0.41195293;   k[0][6]=  1.04028922;   k[0][7]=  0.08327678;
                                                  //theta
         k[1][0]=  0;          k[1][1]=  0;         k[1][2]=  0;        k[1][3]=  3.15163;
         k[1][4]=  6.1119;     k[1][5]=  6.77708;      k[1][6]=  11.32384;    k[1][7]=  27.08792;

         //real gas part

                                                  //n
         K[0][0]=  3.8856823203161E-01; K[0][1]=  2.9385475942740E+00; K[0][2]= -5.5867188534934E+00;
         K[0][3]= -7.6753199592477E-01; K[0][4]=  3.1729005580416E-01; K[0][5]=  5.4803315897767E-01;
         K[0][6]=  1.2279411220335E-01; K[0][7]=  2.1658961543220E+00; K[0][8]=  1.5841735109724E+00;
         K[0][9]= -2.3132705405503E-01; K[0][10]= 5.8116916431436E-02; K[0][11]=-5.5369137205382E-01;
         K[0][12]= 4.8946615909422E-01; K[0][13]=-2.4275739843501E-02; K[0][14]= 6.2494790501678E-02;
         K[0][15]=-1.2175860225246E-01; K[0][16]=-3.7055685270086E-01; K[0][17]=-1.6775879700426E-02;
         K[0][18]=-1.1960736637987E-01; K[0][19]=-4.5619362508778E-02; K[0][20]= 3.5612789270346E-02;
         K[0][21]=-7.4427727132052E-03; K[0][22]=-1.7395704902432E-03; K[0][23]=-2.1810121289527E-02;
         K[0][24]= 2.4332166559236E-02; K[0][25]=-3.7440133423463E-02; K[0][26]= 1.4338715756878E-01;
         K[0][27]=-1.3491969083286E-01; K[0][28]=-2.3151225053480E-02; K[0][29]= 1.2363125492901E-02;
         K[0][30]= 2.1058321972940E-03; K[0][31]=-3.3958519026368E-04; K[0][32]= 5.5993651771592E-03;
         K[0][33]=-3.0335118055646E-04; K[0][34]=-2.1365488688320E+02; K[0][35]= 2.6641569149272E+04;
         K[0][36]=-2.4027212204557E+04; K[0][37]=-2.8341603423999E+02; K[0][38]= 2.1247284400179E+02;
         K[0][39]=-6.6642276540751E-01; K[0][40]= 7.2608632349897E-01; K[0][41]= 5.5068668612842E-02;

                                                  //d
         K[1][0]=  1;   K[1][1]=  1;   K[1][2]=  1;   K[1][3]=  1;   K[1][4]=  2;   K[1][5]=  2;
         K[1][6]=  3;   K[1][7]=  1;   K[1][8]=  2;   K[1][9]=  4;   K[1][10]= 5;   K[1][11]= 5;
         K[1][12]= 5;   K[1][13]= 6;   K[1][14]= 6;   K[1][15]= 6;   K[1][16]= 1;   K[1][17]= 1;
         K[1][18]= 4;   K[1][19]= 4;   K[1][20]= 4;   K[1][21]= 7;   K[1][22]= 8;   K[1][23]= 2;
         K[1][24]= 3;   K[1][25]= 3;   K[1][26]= 5;   K[1][27]= 5;   K[1][28]= 6;   K[1][29]= 7;
         K[1][30]= 8;   K[1][31]= 10;  K[1][32]= 4;   K[1][33]= 8;   K[1][34]= 2;   K[1][35]= 2;
         K[1][36]= 2;   K[1][37]= 3;   K[1][38]= 3;

                                                  //t
         K[2][0]= 0;    K[2][1]= 0.75; K[2][2]= 1;    K[2][3]= 2;    K[2][4]= 0.75; K[2][5]= 2;
         K[2][6]= 0.75; K[2][7]= 1.5;  K[2][8]= 1.5;  K[2][9]= 2.5;  K[2][10]= 0;   K[2][11]= 1.5;
         K[2][12]= 2;   K[2][13]= 0;   K[2][14]= 1;   K[2][15]= 2;   K[2][16]= 3;   K[2][17]= 6;
         K[2][18]= 3;   K[2][19]= 6;   K[2][20]= 8;   K[2][21]= 6;   K[2][22]= 0;   K[2][23]= 7;
         K[2][24]= 12;  K[2][25]= 16;  K[2][26]= 22;  K[2][27]= 24;  K[2][28]= 16;  K[2][29]= 24;
         K[2][30]= 8;   K[2][31]= 2;   K[2][32]= 28;  K[2][33]= 14;  K[2][34]= 1;   K[2][35]= 0;
         K[2][36]= 1;   K[2][37]= 3;   K[2][38]= 3;

                                                  //c
         K[3][7]= 1;    K[3][8]= 1;    K[3][9]= 1;    K[3][10]= 1;   K[3][11]= 1;   K[3][12]= 1;
         K[3][13]= 1;   K[3][14]= 1;   K[3][15]= 1;   K[3][16]= 2;   K[3][17]= 2;   K[3][18]= 2;
         K[3][19]= 2;   K[3][20]= 2;   K[3][21]= 2;   K[3][22]= 2;   K[3][23]= 3;   K[3][24]= 3;
         K[3][25]= 3;   K[3][26]= 4;   K[3][27]= 4;   K[3][28]= 4;   K[3][29]= 4;   K[3][30]= 4;
         K[3][31]= 4;   K[3][32]= 5;   K[3][33]= 6;

                                                  //a
         K[4][39]= 3.5;    K[4][40]= 3.5;    K[4][41]= 3;
                                                  //b
         K[5][39]= 0.875;  K[5][40]= 0.925;  K[5][41]= 0.875;
                                                  //A
         K[6][39]=0.7;     K[6][40]=0.7;     K[6][41]=0.7;
                                                  //B
         K[7][39]=0.3;     K[7][40]=0.3;     K[7][41]=1;
                                                  //C
         K[8][39]= 10;     K[8][40]= 10;     K[8][41]= 12.5;
                                                  //D
         K[9][39]=275;     K[9][40]=275;     K[9][41]=275;

                                                  //alpha
         K[10][34]=25;  K[10][35]=25;     K[10][36]=25;     K[10][37]=15;  K[10][38]=20;
                                                  //beta
         K[11][34]=325; K[11][35]=300;    K[11][36]=300;    K[11][37]=275; K[11][38]=275;
         K[11][39]=0.3; K[11][40]=0.3;    K[11][41]=0.3;
                                                  //gamma
         K[12][34]=1.16; K[12][35]=1.19;     K[12][36]=1.19;      K[12][37]=1.25; K[12][38]=1.22;
                                                  // epsilon
         K[13][34]=1;   K[13][35]=1;      K[13][36]=1;      K[13][37]=1;   K[13][38]=1;

         break;
      }
      case 'W' :                                  // WATER
      {
         fluid_id = 1;
         rhoc=322;                                //[kg/m3]
         Tc=647.096;                              //[K]
         pc=22064000;                             // [Pa]
         Tt=273.16;                               //  [K]
         pt=611.657;                              //  [Pa]
         Rs=461.51805;                            //  [J/kg/K]
         molar_mass=18.01528;                     //  [g/mol]
         omega=0.344;                             // azentric factor, see PREOS

         // Limits for Sums in FHE-derivations

         limit[0]=7;
         limit[1]=51;
         limit[2]=54;
         limit[3]=56;

         // Coefficients for FHE-derivations

         for (i=0;i<2;i++)
            for (j=0;j<8;j++)
               {k[i][j] = 0;}

               for (i=0;i<14;i++)
                  for (j=0;j<56;j++)
                     {K[i][j] = 0;}

                     // ideal gas part

                     k[0][0]=  -8.32044648201;  k[0][1]=   6.6832105268;   k[0][2]=  3.00632;
         k[0][3]=  0.012436;        k[0][4]=  0.97315;         k[0][5]=  1.27950;
         k[0][6]=  0.96956;         k[0][7]=  0.24873;

         k[1][0]=  0;            k[1][1]=  0;            k[1][2]=  0;
         k[1][3]=  1.28728967;      k[1][4]=  3.53734222;      k[1][5]=  7.74073708;
         k[1][6]=  9.24437796;      k[1][7]=  27.5075105;

         //real gas part

         K[0][0]=  1.2533547935523E-02;K[0][1]=  7.8957634722828E+00;K[0][2]= -8.7803203303561E+00;
         K[0][3]=  3.1802509345418E-01;K[0][4]= -2.6145533859358E-01;K[0][5]= -7.8199751687981E-03;
         K[0][6]=  8.8089493102134E-03;K[0][7]= -6.6856572307965E-01;K[0][8]=  2.0433810950965E-01;
         K[0][9]= -6.6212605039687E-05;K[0][10]=-1.9232721156002E-01;K[0][11]=-2.5709043003438E-01;
         K[0][12]= 1.6074868486251E-01;K[0][13]=-4.0092828925807E-02;K[0][14]= 3.9343422603254E-07;

         K[0][15]=-7.5941377088144E-06;K[0][16]= 5.6250979351888E-04;K[0][17]=-1.5608652257135E-05;
         K[0][18]= 1.1537996422951E-09;K[0][19]= 3.6582165144204E-07;K[0][20]=-1.3251180074668E-12;
         K[0][21]=-6.2639586912454E-10;K[0][22]=-1.0793600908932E-01;K[0][23]= 1.7611491008752E-02;
         K[0][24]= 2.2132295167546E-01;K[0][25]=-4.0247669763528E-01;K[0][26]= 5.8083399985759E-01;
         K[0][27]= 4.9969146990806E-03;K[0][28]=-3.1358700712549E-02;K[0][29]=-7.4315929710341E-01;

         K[0][30]= 4.7807329915480E-01;K[0][31]= 2.0527940895948E-02;K[0][32]=-1.3636435110343E-01;
         K[0][33]= 1.4180634400617E-02;K[0][34]= 8.3326504880713E-03;K[0][35]=-2.9052336009585E-02;
         K[0][36]= 3.8615085574206E-02;K[0][37]=-2.0393486513704E-02;K[0][38]=-1.6554050063734E-03;
         K[0][39]= 1.9955571979541E-03;K[0][40]= 1.5870308324157E-04;K[0][41]=-1.6388568342530E-05;
         K[0][42]= 4.3613615723811E-02;K[0][43]= 3.4994005463765E-02;K[0][44]=-7.6788197844621E-02;

         K[0][45]= 2.2446277332006E-02;K[0][46]=-6.2689710414685E-05;K[0][47]=-5.5711118565645E-10;
         K[0][48]=-1.9905718354408E-01;K[0][49]= 3.1777497330738E-01;K[0][50]=-1.1841182425981E-01;
         K[0][51]=-3.1306260323435E+01;K[0][52]= 3.1546140237781E+01;K[0][53]=-2.5213154341695E+03;
         K[0][54]=-1.4874640856724E-01;K[0][55]= 3.1806110878444E-01;

         K[1][0]=1;  K[1][1]=1;  K[1][2]=1;  K[1][3]=2;  K[1][4]=2;  K[1][5]=3;  K[1][6]=4;
         K[1][7]=1;  K[1][8]=1;  K[1][9]=1;  K[1][10]=2; K[1][11]=2; K[1][12]=3; K[1][13]=4;
         K[1][14]=4; K[1][15]=5; K[1][16]=7; K[1][17]=9; K[1][18]=10;K[1][19]=11;K[1][20]=13;
         K[1][21]=15;K[1][22]=1; K[1][23]=2; K[1][24]=2; K[1][25]=2; K[1][26]=3; K[1][27]=4;
         K[1][28]=4; K[1][29]=4; K[1][30]=5; K[1][31]=6; K[1][32]=6; K[1][33]=7; K[1][34]=9;
         K[1][35]=9; K[1][36]=9; K[1][37]=9; K[1][38]=9; K[1][39]=10;K[1][40]=10;K[1][41]=12;
         K[1][42]=3; K[1][43]=4; K[1][44]=4; K[1][45]=5; K[1][46]=14;K[1][47]=3; K[1][48]=6;
         K[1][49]=6; K[1][50]=6; K[1][51]=3; K[1][52]=3; K[1][53]=3;

         K[2][0]=-0.5;     K[2][1]=0.875;    K[2][2]=1;     K[2][3]=0.5;      K[2][4]=0.75;
         K[2][5]=0.375;    K[2][6]=1;        K[2][7]=4;     K[2][8]=6;        K[2][9]=12;
         K[2][10]=1;       K[2][11]=5;       K[2][12]=4;    K[2][13]=2;       K[2][14]=13;
         K[2][15]=9;       K[2][16]=3;       K[2][17]=4;    K[2][18]=11;      K[2][19]=4;
         K[2][20]=13;      K[2][21]=1;       K[2][22]=7;    K[2][23]=1;       K[2][24]=9;
         K[2][25]=10;      K[2][26]=10;      K[2][27]=3;    K[2][28]=7;       K[2][29]=10;
         K[2][30]=10;      K[2][31]=6;       K[2][32]=10;   K[2][33]=10;      K[2][34]=1;
         K[2][35]=2;       K[2][36]=3;       K[2][37]=4;    K[2][38]=8;       K[2][39]=6;
         K[2][40]=9;       K[2][41]=8;       K[2][42]=16;   K[2][43]=22;      K[2][44]=23;
         K[2][45]=23;      K[2][46]=10;      K[2][47]=50;   K[2][48]=44;      K[2][49]=46;
         K[2][50]=50;      K[2][51]=0;       K[2][52]=1;    K[2][53]=4;

         K[3][7]=1;  K[3][8]=1;  K[3][9]=1;  K[3][10]=1; K[3][11]=1; K[3][12]=1; K[3][13]=1; K[3][14]=1;
         K[3][15]=1; K[3][16]=1; K[3][17]=1; K[3][18]=1; K[3][19]=1; K[3][20]=1; K[3][21]=1; K[3][22]=2;
         K[3][23]=2; K[3][24]=2; K[3][25]=2; K[3][26]=2; K[3][27]=2; K[3][28]=2; K[3][29]=2; K[3][30]=2;
         K[3][31]=2; K[3][32]=2; K[3][33]=2; K[3][34]=2; K[3][35]=2; K[3][36]=2; K[3][37]=2; K[3][38]=2;
         K[3][39]=2; K[3][40]=2; K[3][41]=2; K[3][42]=3; K[3][43]=3; K[3][44]=3; K[3][45]=3; K[3][46]=4;
         K[3][47]=6; K[3][48]=6; K[3][49]=6; K[3][50]=6;

         K[4][54] = 3.5;      K[4][55] = 3.5;
         K[5][54] = 0.85;  K[5][55] = 0.95;
         K[6][54] = 0.32;  K[6][55] = 0.32;
         K[7][54] = 0.2;      K[7][55] = 0.2;
         K[8][54] = 28;    K[8][55] = 32;
         K[9][54] = 700;      K[9][55] = 800;

         K[10][51] = 20;      K[10][52] = 20;      K[10][53] = 20;
         K[11][51] = 150;  K[11][52] = 150;  K[11][53] = 250;  K[11][54] = 0.3;  K[11][55] = 0.3;
         K[12][51] = 1.21; K[12][52] = 1.21; K[12][53] = 1.25;
         K[13][51] = 1; K[13][52] = 1; K[13][53] = 1;

         break;
      }
      case 'M' :                                  // METHANE
      {
         fluid_id = 2;
         rhoc=162.66;                             //[kg/m3]
         Tc=190.551;                              //[K]
         pc=4599200;                              // [Pa]
         Tt=90.685;                               //  [K]
         pt=11696;                                //  [Pa]
         Rs=518.3;                                //  [J/kg/K]
         molar_mass=16.04;                        //  [g/mol]
         omega=0.011;                             // azentric factor, see PREOS

         // Limits sums in FHE-derivations

         limit[0]=13;
         limit[1]=36;
         limit[2]=40;
         limit[3]=40;

         // Coefficients for FHE-derivations

         for (i=0;i<2;i++)
            for (j=0;j<8;j++)
               {k[i][j] = 0;}

               for (i=0;i<14;i++)
                  for (j=0;j<56;j++)
                     {K[i][j] = 0;}

                     // ideal gas part

                     k[0][0]=9.91243972;    k[0][1]=-6.33270087;    k[0][2]=3.0016;        k[0][3]=0.008449;            k[0][4]=4.6942;
         k[0][5]=3.4865;    k[0][6]=1.6572;    k[0][7]=1.4115;

         k[1][0]=0; k[1][1]=0;k[1][2]=0;k[1][3]=3.4004324;k[1][4]=10.26951575;k[1][5]=20.43932747;k[1][6]=29.93744884;k[1][7]=79.13351945;

         //real gas part

         K[0][0]=4.368E-02; K[0][1]=6.709E-01; K[0][2]=-1.766E+00; K[0][3]=8.582E-01; K[0][4]=-1.207E+00; K[0][5]=5.120E-01; K[0][6]=-4.000E-04;
         K[0][7]=-1.248E-02; K[0][8]=3.100E-02; K[0][9]=1.755E-03; K[0][10]=-3.172E-06; K[0][11]=-2.240E-06; K[0][12]=2.947E-07; K[0][13]=1.830E-01;
         K[0][14]=1.512E-01; K[0][15]=-4.289E-01; K[0][16]=6.894E-02; K[0][17]=-1.408E-02; K[0][18]=-3.063E-02; K[0][19]=-2.970E-02; K[0][20]=-1.932E-02;
         K[0][21]=-1.106E-01; K[0][22]=9.953E-02; K[0][23]=8.548E-03; K[0][24]=-6.151E-02; K[0][25]=-4.292E-02; K[0][26]=-1.813E-02; K[0][27]=3.446E-02;
         K[0][28]=-2.386E-03;  K[0][29]=-1.159E-02; K[0][30]=6.642E-02; K[0][31]=-2.372E-02; K[0][32]=-3.962E-02; K[0][33]=-1.387E-02; K[0][34]=3.389E-02;
         K[0][35]=-2.927E-03; K[0][36]=9.325E-05; K[0][37]=-6.287E+00; K[0][38]=1.271E+01; K[0][39]=-6.424E+00;

         K[3][13]=1; K[3][14]=1; K[3][15]=1; K[3][16]=1; K[3][17]=1; K[3][18]=1; K[3][19]=1; K[3][20]=2; K[3][21]=2; K[3][22]=2; K[3][23]=2; K[3][24]=2;
         K[3][25]=3; K[3][26]=3; K[3][27]=3; K[3][28]=3; K[3][29]=4; K[3][30]=4; K[3][31]=4; K[3][32]=4; K[3][33]=4; K[3][34]=4; K[3][35]=4;

         K[1][0]=1;K[1][1]=1;K[1][2]=1;K[1][3]=2;K[1][4]=2;K[1][5]=2;K[1][6]=2;K[1][7]=3;K[1][8]=4;K[1][9]=4;K[1][10]=8;K[1][11]=9;K[1][12]=10;K[1][13]=1;
         K[1][14]=1;K[1][15]=1;K[1][16]=2;K[1][17]=4;K[1][18]=5;K[1][19]=6;K[1][20]=1;K[1][21]=2;K[1][22]=3;K[1][23]=4;K[1][24]=4;K[1][25]=3;K[1][26]=5;
         K[1][27]=5;K[1][28]=8;K[1][29]=2;K[1][30]=3;K[1][31]=4;K[1][32]=4;K[1][33]=4;K[1][34]=5;K[1][35]=6;K[1][36]=2;K[1][37]=0;K[1][38]=0;K[1][39]=0;

         K[2][0]=-0.5; K[2][1]=0.5; K[2][2]=1;K[2][3]=0.5;K[2][4]=1;K[2][5]=1.5;K[2][6]=4.5;K[2][7]=0;K[2][8]=1;K[2][9]=3;K[2][10]=1;K[2][11]=3;K[2][12]=3;
         K[2][13]=0;K[2][14]=1;K[2][15]=2;K[2][16]=0;K[2][17]=0;K[2][18]=2;K[2][19]=2;K[2][20]=5;K[2][21]=5;K[2][22]=5;K[2][23]=2;K[2][24]=4;K[2][25]=12;K[2][26]=8;
         K[2][27]=10;K[2][28]=10;K[2][29]=10;K[2][30]=14;K[2][31]=12;K[2][32]=18;K[2][33]=22;K[2][34]=18;K[2][35]=14;K[2][36]=2;K[2][37]=0;K[2][38]=1;K[2][39]=2;
         K[10][36]=20;K[10][37]=40;K[10][38]=40;K[10][39]=40;
         K[11][36]=200;K[11][37]=250;K[11][38]=250;K[11][39]=250;
         K[12][36]=1.07;K[12][37]=1.11;K[12][38]=1.11;K[12][39]=1.11;
         K[13][36]=1;K[13][37]=1;K[13][38]=1;K[13][39]=1;

         break;
      }
      case 'N' :                                  // Nitrogen
      {
         fluid_id = 3;
         rhoc=314.0;                              //[kg/m3]
         Tc=126.20;                               //[K]
         pc=3383000;                              // [Pa]
         Tt=63.148;                               //  [K]
         pt=12500;                                //  [Pa]
         Rs=296.8;                                //  [J/kg/K]
         molar_mass=28.013;                       //  [g/mol]
         omega=0.039;                             // azentric factor, see PREOS

         // Limits sums in FHE-derivations
         limit[0]=6;
         limit[1]=32;
         limit[2]=36;
         limit[3]=36;

         // Coefficients for FHE-derivations
         // ideal gas part
         k[0][0]=  9.912644;
         k[0][1]=  -6.333133;
         k[0][2]=  3.0016;
         k[0][3]=  0.008449;
         k[0][4]=  4.6942;
         k[0][5]=  3.4865;
         k[0][6]=  1.6572;
         k[0][7]=  1.4115;

         k[1][0] = 0;  k[1][1] = 0;  k[1][2] = 0; k[1][3] = 3.400664; k[1][4] = 10.27022;
         k[1][5] = 20.44072; k[1][6] = 29.93949; k[1][7] = 79.13892;
         //real gas part
         K[0][0]=  0.924803575275;        K[0][1]= -0.492448489428;        K[0][2]=  0.661883336938;        K[0][3]= -0.192902649201e1;
         K[0][4]= -0.622469309629e-1;     K[0][5]=  0.349943957581;        K[0][6]=  0.564857472498;        K[0][7]= -0.161720005987e1;
         K[0][8]= -0.481395031883 ;          K[0][9]=  0.421150636384;        K[0][10]= -0.161962230825e-1;       K[0][11]= 0.172100994165;
         K[0][12]= 0.735448924933e-2;     K[0][13]= 0.168077305479e-1;     K[0][14]= -0.107626664179e-2;       K[0][15]= -0.137318088513e-1;
         K[0][16]= 0.635466899859e-3;     K[0][17]= 0.304432279419e-2;     K[0][18]= -0.435762336045e-1;       K[0][19]= -0.723174889316e-1;
         K[0][20]= 0.389644315272e-1;     K[0][21]= -0.212201363910e-1;       K[0][22]= 0.408822981509e-2;     K[0][23]= -0.551990017984e-4;
         K[0][24]= -0.462016716479e-1;       K[0][25]= -0.300311716011e-2;       K[0][26]= 0.368825891208e-1;     K[0][27]= -0.255856846220e-2;
         K[0][28]= 0.896915264558e-2;     K[0][29]= -0.441513370350e-2;       K[0][30]= 0.133722924858e-2;     K[0][31]= 0.264832491957e-3;
         K[0][32]= 0.196688194015e2;         K[0][33]= -0.209115600730e2;     K[0][34]= 0.167788306989e-1;     K[0][35]= 0.262767566274e4;

         K[1][0] = 1;         K[1][1] = 1;         K[1][2] = 2;         K[1][3] = 2;         K[1][4] = 3;         K[1][5] = 3;
         K[1][6] = 1;         K[1][7] = 1;         K[1][8] = 1;         K[1][9] = 3;         K[1][10] =3;         K[1][11] =4;
         K[1][12] =6;         K[1][13] =6;         K[1][14] =7;         K[1][15] =7;         K[1][16] =8;         K[1][17] =8;
         K[1][18] =1;         K[1][19] =2;         K[1][20] =3;         K[1][21] =4;         K[1][22] =5;         K[1][23] =8;
         K[1][24] =4;         K[1][25] =5;         K[1][26] =5;         K[1][27] =8;         K[1][28] =3;         K[1][29] =5;
         K[1][30] =6;          K[1][31] =9;         K[1][32] =1;         K[1][33] =1;         K[1][34] =3;         K[1][35] =2;

         K[2][00]=0.25; K[2][01]=0.875; K[2][02]=0.5; K[2][03]=0.875; K[2][04]=0.375; K[2][05]=0.75; K[2][6]=0.5;
         K[2][7]=0.75;  K[2][8]=2;      K[2][9]=1.25; K[2][10]=3.5;   K[2][11]=1;

         K[2][12]=0.5;          K[2][13]=3;           K[2][14]=0;            K[2][15]=2.75;         K[2][16]=0.75;            K[2][17]=2.5;
         K[2][18]=4;            K[2][19]=6;           K[2][20]=6;            K[2][21]=3;            K[2][22]=3;               K[2][23]=6;
         K[2][24]=16;           K[2][25]=11;          K[2][26]=15;           K[2][27]=12;           K[2][28]=12;              K[2][29]=7;
         K[2][30]=4;            K[2][31]=16;          K[2][32]=0;            K[2][33]=1;            K[2][34]=2;               K[2][35]=3;

         K[3][00] = 0;            K[3][01] = 0;            K[3][02] = 0;            K[3][03] = 0;            K[3][04] = 0;            K[3][05] = 0;
         K[3][06] = 1;            K[3][07] = 1;            K[3][8] = 1;             K[3][9] = 1;             K[3][10] = 1;            K[3][11] = 1;
         K[3][12] = 1;            K[3][13] = 1;            K[3][14] = 1;            K[3][15] = 1;            K[3][16] = 1;            K[3][17] = 1;
         K[3][18] = 2;            K[3][19] = 2;            K[3][20] = 2;            K[3][21] = 2;            K[3][22] = 2;            K[3][23] = 2;
         K[3][24] = 3;            K[3][25] = 3;            K[3][26] = 3;            K[3][27] = 3;            K[3][28] = 4;            K[3][29] = 4;
         K[3][30] = 4;            K[3][31] = 4;            K[3][32] = 2;            K[3][33] = 2;            K[3][34] = 2;            K[3][35] = 2;

         K[10][32] = 20;     K[10][33] = 20;     K[10][34] = 15;    K[10][35] = 25;
         K[11][32] = 325;    K[11][33] = 325;    K[11][34] = 300;    K[11][35] = 275;
         K[12][32] = 1.16;   K[12][33] = 1.16;   K[12][34] = 1.13;   K[12][35] = 1.25;
         K[13][32] = 1;      K[13][33] = 1;      K[13][34] = 1;      K[13][35] = 1;
         break;
      }
      default: cout << "Error in eos.cpp: no fluid name specified!" << endl;
      break;
   }
}
