//#include <math.h>
//#include <fstream>
//#include <iostream>
//#include <stdio.h>
//#include <string>

double pressure (double rho, double T, int fluid);
// calculates the pressure depending on density and temperature
//double density (double P, double rho0, double T, double prec, string c);
// calculates the density iteratively depending on pressure and temperature
double enthalpy (double rho, double T, std::string c);
// calculates the density iteratively depending on pressure and temperature
double isochoric_heat_capacity (double rho, double T, std::string c);
// calculates the isochoric heat capacity depending on density and temperature
double isobaric_heat_capacity (double rho, double T, std::string c);
// calculates the isobaric heat capacity depending on density and temperature
double co2_viscosity (double rho, double T);
// calculates the viscosity depending on density and temperature !ONLY for CO2!!!
double co2_heat_conductivity (double, double);
// calculates the heat conductivity of co2 depending on density and temperature
double ch4_viscosity_295K (double);
// calculates the viscosity of CH4 at 25 �C depending on pressure
double Fluid_Viscosity (double rho, double T, double p, int fluid);
//Viscosity for several fluids
double Fluid_Heat_Conductivity (double rho , double T, int fluid);
//Heat conductivity for several fluids
double rkeos(double T, double P, double MM, double a, double b);
double melting_pressure_co2(double T,double Tt,double pt);
double sublime_pressure_co2(double T,double Tt,double pt);
double vapour_pressure_co2(double T);
double rkeos(double T, double P, int fluid);
double preos(double T, double P, int fluid);
double h2o_viscosity_IAPWS (double rho, double T);
double h2o_heat_conductivity_IAPWS_ind (double rho, double T);
double ch4_viscosity (double rho, double T);
double ch4_heat_conductivity (double rho, double T);
double n2_viscosity (double rho, double T);
double n2_heat_conductivity (double rho, double T);
double mixing_ternary (double* x,double* a,double* b,double *MM, double *ra, double *rb, double *rMM);

double dpressure(double TT, double PP,std::string cs, double ds);
double zbrent(double TT, double PP, int fluid, const double tol);

double vapour_saturation_density_ch4(double T);
double liquid_saturation_density_ch4(double T);

double vapour_pressure_n2(double T);
double liquid_saturation_density_n2(double T);
double vapour_saturation_density_n2(double T);
double vapour_saturation_density_co2(double T);
double liquid_saturation_density_co2(double T);
