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

    std::string output_name = "output_" + planet_name + ".csv";
    std::ofstream output;
    output.open(output_name.c_str());
    output << "time" << ',' << "gamma" << ',' << "g-dot" << ',' << "theta" << ',' << "theta-dot" << ',' << "mean-anomaly" << ','
        << "mean-motion" << ',' << "n-dot" << ',' << "eccentricity" << std::endl;

    double run_time = 9000; //in simulation years

    double last_printed_time = 0.0;
    double dt = 2*PI / planet.mean_motion(0.0) / 2.0;
    while(planet.get_time() < run_time)
    {
        planet.solve();
        if(planet.get_time() - last_printed_time >= dt)
        {
            output << std::setprecision(8) << planet.get_time() << ',' << std::setprecision(6) << planet.get_gamma() << ',' << planet.get_gamma_dot()  << ','
                << planet.get_theta() << ',' << planet.get_theta_dot() << ',' << planet.get_mean_anomaly() << ','
                << planet.mean_motion(planet.get_time()) << ',' << planet.mean_motion_dot(planet.get_time()) << ',' << planet.eccentricity(planet.get_time())
                << std::endl;
            last_printed_time = planet.get_time();
        }

    }


    return 0;
}
