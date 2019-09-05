#ifndef DAMP_H_INCLUDED
#define DAMP_H_INCLUDED

#include "planet.h"
#include <iostream>

//double damp(double M_star, double plan_radius, double semi_major, double ecc, double M_plan, double mean_motion, double gdot, double mu, double alpha, double tau_M, double tau_A);
double damp(Planet *planet, double t);
double H(int p, double ecc);
double K_c(Planet *planet, double chi);
//double K_c(double chi, double radius, double M_plan, double mu, double alpha, double tau_M, double tau_A);


#endif // DAMP_H_INCLUDED
