#include <iostream>

#include "damp.h"
#include "constants.h"


int main()
{

    double damping;

    double M_star = MSUN;
    double plan_radius = REARTH;
    double semi_major = AEARTH;
    double ecc = 0.02;
    double M_plan = MEARTH;
    double mean_motion = 2 * PI / (365*24*3600.);
    double gdot = 0.0;
    double mu = 0.8 * pow(10,11);
    double alpha = 0.2;
    double tau_M = 500 * 365 * 24 * 3600.;
    double tau_A = tau_M;
    double B_A_C = 10e-5;
    double mass_star = 0.08 * MSUN;

    std::string name = "test planet";
    Planet planet(name, plan_radius, M_plan, B_A_C, mu, alpha, tau_M, tau_A, mass_star);
    planet.solve();

    return 0;
}
