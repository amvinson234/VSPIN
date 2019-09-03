#include "planet.h"
#include <cmath>


Planet::Planet(std::string name, double radius, double mass, double B_A_C, double mu, double alpha, double tau_M, double tau_A, double mass_star)
{
    _name = name;
    _radius = radius;
    _mass = mass;
    _B_A_C = B_A_C;
    _mu = mu;
    _alpha = alpha;
    _tau_M = tau_M;
    _tau_A = tau_A;
    _mass_star = mass_star;

    gamma = 0.;
    gamma_dot = 0.;

    //initialize to Earth values if other planet details not specified
    semi_major = AEARTH;
    _mean_motion = 2 * PI / 365. / 24. / 3600.;
    ecc = 0.0167;

    time = 0.;
    time_step = 365. * 24. * 3600.;

    _min_dt = 0.;
    _max_dt = INFINITY;
}

std::string Planet::name()
{
    return _name;
}

double Planet::get_gamma()
{
    return gamma;
}
double Planet::get_gamma_dot()
{
    return gamma_dot;
}

double Planet::get_ecc()
{
    return ecc;
}
double Planet::get_mean_motion()
{
    return _mean_motion;
}
double Planet::get_semi_major()
{
    return semi_major;
}

double Planet::get_radius()
{
    return _radius;
}
double Planet::get_mass()
{
    return _mass;
}
double Planet::get_mu()
{
    return _mu;
}
double Planet::get_alpha()
{
    return _alpha;
}
double Planet::get_tau_M()
{
    return _tau_M;
}
double Planet::get_tau_A()
{
    return _tau_A;
}
double Planet::get_time()
{
    return time;
}
