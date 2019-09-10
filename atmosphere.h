#ifndef ATMOSPHERE_H_INCLUDED
#define ATMOSPHERE_H_INCLUDED

#include "planet.h"

class Atmosphere
{
public:
    Atmosphere();
    Atmosphere(Planet planet);

    double damp();


private:
    Planet _planet;
};

#endif // ATMOSPHERE_H_INCLUDED
