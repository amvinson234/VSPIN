#include <iostream>

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
    double B_A_C = 10e-5;
    double mass_star = 0.08 * MSUN;

    std::string name = "trappf";
    Planet planet(name, plan_radius, M_plan, B_A_C, mu, alpha, tau_M, tau_A, mass_star);
    std::cerr << planet.get_time() << '\t' << planet.get_gamma() << '\t' << planet.get_gamma_dot() << std::endl;
    planet.solve();
    std::cerr << planet.get_time() << '\t' << planet.get_gamma() << '\t' << planet.get_gamma_dot() << std::endl;
    planet.solve();
    std::cerr << planet.get_time() << '\t' << planet.get_gamma() << '\t' << planet.get_gamma_dot() << std::endl;
    planet.solve();
    std::cerr << planet.get_time() << '\t' << planet.get_gamma() << '\t' << planet.get_gamma_dot() << std::endl;
    planet.solve();
    std::cerr << planet.get_time() << '\t' << planet.get_gamma() << '\t' << planet.get_gamma_dot() << std::endl;
    planet.solve();
    std::cerr << planet.get_time() << '\t' << planet.get_gamma() << '\t' << planet.get_gamma_dot() << std::endl;
    planet.solve();
    std::cerr << planet.get_time() << '\t' << planet.get_gamma() << '\t' << planet.get_gamma_dot() << std::endl;
    planet.solve();
    std::cerr << planet.get_time() << '\t' << planet.get_gamma() << '\t' << planet.get_gamma_dot() << std::endl;

    return 0;
}
