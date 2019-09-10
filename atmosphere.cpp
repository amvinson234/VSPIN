#include "atmosphere.h"

Atmosphere::Atmosphere()
{

}

Atmosphere::Atmosphere(Planet planet)
{
    _planet = planet;
}

double Atmosphere::damp()
{
    double mean_motion
    double K_a; //still needs to be initialized
    double b_a; //still needs to be initialized

    return -3.0/2.0 * 3.0 * K_a * b_a * 2 * (omega - mean_motion);
}

