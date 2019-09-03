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
    double gdot = 1.0 / (365*24*3600.);
    double mu = 0.8 * pow(10,11);
    double alpha = 0.2;
    double tau_M = 500 * 365 * 24 * 3600.;
    double tau_A = tau_M;

    std::string name = "test planet";
    Planet planet(name, plan_radius, M_plan, mu, alpha, tau_M, tau_A);

    return 0;
}
