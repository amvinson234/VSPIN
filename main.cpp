#include <iostream>
#include <iomanip>
#include "damp.h"
#include "constants.h"


int main()
{

    double plan_radius = REARTH;
    double M_plan = MEARTH;
    double mu = 0.8 * pow(10,11);
    double alpha = 0.2;
    double tau_M = 500;
    double tau_A = tau_M;
    double B_A_C = 10.e-5;
    double mass_star = 0.08 * MSUN;

    std::string name = "trappf";
    Planet planet(name, plan_radius, M_plan, B_A_C, mu, alpha, tau_M, tau_A, mass_star);

    std::string output_name = "test_output.txt";
    std::ofstream output;
    output.open(output_name.c_str());
    output << "time" << '\t' << "gamma" << '\t' << "g-dot" << '\t' << "mean motion" << '\t' << "n-dot" << '\t' << "eccentricity" << std::endl;

    while(planet.get_time() < 1000)
    {
        planet.solve();
        output << std::setprecision(8) << planet.get_time() << '\t' << std::setprecision(6) << planet.get_gamma() << '\t' << planet.get_gamma_dot()  << '\t'
        << '\t' << planet.mean_motion(planet.get_time()) << '\t' << planet.eccentricity(planet.get_time()) << std::endl;
    }



    return 0;
}
