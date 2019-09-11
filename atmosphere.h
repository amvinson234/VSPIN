#ifndef ATMOSPHERE_H_INCLUDED
#define ATMOSPHERE_H_INCLUDED

#include "planet.h"
#include "constants.h"

class Atmosphere
{
public:
    Atmosphere();
    Atmosphere(Planet *planet);

    double damp();

private:
    Planet *_planet;

    double _mean_density;

    //Following taken from Table 1, Leconte et al. 2015
    double _q_0 = (365. * 24. * 3600.) * (365. * 24. * 3600.) * 1000.0; //amplitude of the atmospheric quadrupole. pressure, SI except time in years.
    double _omega_0 = 71.0; //intrinsic thermal frequency of the atmosphere
    double b_a(double freq);
};

#endif // ATMOSPHERE_H_INCLUDED
