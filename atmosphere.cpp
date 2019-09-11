#include "atmosphere.h"
#include "planet.h"

Atmosphere::Atmosphere()
{

}

Atmosphere::Atmosphere(Planet *planet)
{
    _planet = planet;
    _mean_density = planet->get_mass() / (4.0 / 3.0 * PI * std::pow(planet->get_radius(), 3));
    K_a = 3 * planet->get_stellar_mass() * std::pow(planet->get_radius(), 3) / (5.0 * _mean_density * std::pow(planet->get_semi_major(), 3));
}

double Atmosphere::damp(Planet *planet)
{
    double sigma = 2.0 * planet->get_gamma_dot();

    return -3.0/2.0 * 3.0 * K_a * b_a(sigma);
}

double Atmosphere::b_a(double freq)
{
    double b = - std::sqrt(10.0 / 3 / PI) * _q_0 * freq / 2.0 / _omega_0 / (1.0 + freq * freq / 4.0 /_omega_0 /_omega_0);
    return b;
}

