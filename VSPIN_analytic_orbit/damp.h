#ifndef DAMP_H_INCLUDED
#define DAMP_H_INCLUDED

#include "planet.h"
#include <iostream>

double damp(Planet *planet, double t);
double H(int p, double ecc);
double K_c(Planet *planet, double chi);
double tau_a(double chi);


#endif // DAMP_H_INCLUDED
