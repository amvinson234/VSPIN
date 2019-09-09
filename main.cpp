#include <iostream>
#include <iomanip>
#include "damp.h"
#include "constants.h"


int main()
{
    std::string input_file_name = "parameters.txt";
    std::ifstream input_params(input_file_name.c_str());
    std::string header_trash;
    std::getline(input_params, header_trash);

    std::string planet_name;
    input_params >> planet_name;
    double input;
    std::vector<double> inputs;
    while(input_params >> input)
    {
        inputs.push_back(input);
    }

    double mass_star = inputs[2] * MSUN;

    double plan_radius, M_plan;
    if(planet_name == "d")
    {
        plan_radius = 0.784 * REARTH;
        M_plan = 0.297 * MEARTH;
    }
    else if (planet_name == "e")
    {
        plan_radius = 0.91 * REARTH;
        M_plan = 0.772 * MEARTH;
    }
    else if (planet_name == "f")
    {
        plan_radius = 1.045 * REARTH;
        M_plan = 0.68 * MEARTH;
    }
    else
    {
        std::cout << "Input planet mass in Earth mass units: ";
        std::cin >> M_plan;
        M_plan *= MEARTH;
        std::cout << std::endl << "Input planet radius in Earth radius units: ";
        std::cin >> plan_radius;
        plan_radius *= REARTH;
    }

    //std::string name = "trapp" + planet_name;
    std::string name = planet_name;
    Planet planet(name, M_plan, plan_radius, inputs);

    std::string output_name = "output_" + planet_name + ".txt";
    std::ofstream output;
    output.open(output_name.c_str());
    output << "time" << '\t' << "gamma" << '\t' << "g-dot" << '\t' << "mean motion" << '\t' << "n-dot" << '\t' << "eccentricity" << std::endl;

    double run_time = 9000; //in simulation years

    while(planet.get_time() < run_time)
    {
        planet.solve();
        output << std::setprecision(8) << planet.get_time() << '\t' << std::setprecision(6) << planet.get_gamma() << '\t' << planet.get_gamma_dot()  << '\t'
        << '\t' << planet.mean_motion(planet.get_time()) << '\t' << planet.eccentricity(planet.get_time()) << std::endl;
    }



    return 0;
}
