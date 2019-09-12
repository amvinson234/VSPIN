#ifndef ATMOSPHERE_H_INCLUDED
#define ATMOSPHERE_H_INCLUDED

#include "constants.h"

class Planet;

class Atmosphere
{
public:
    Atmosphere();
    Atmosphere(Planet *planet);
    Atmosphere(Planet *planet, double q_0, double omega_0);

    double damp(Planet *planet);

private:
    Planet *_planet;

    double _mean_density;

    //Following taken from Table 1, Leconte et al. 2015
    double _q_0; //amplitude of the atmospheric quadrupole. pressure, SI except time in years.
    double _omega_0; //intrinsic thermal frequency of the atmosphere
    double K_a;

    double b_a(double freq);
};

#endif // ATMOSPHERE_H_INCLUDED
