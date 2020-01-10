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

    double input;
    std::vector<double> inputs;
    while(input_params >> input)
    {
        inputs.push_back(input);
    }


    input_file_name = "orbital_parameters.txt";
    std::ifstream orbital_params(input_file_name.c_str());
    std::getline(orbital_params, header_trash);

    std::vector<double> orbit_inputs;
    while(orbital_params >> input)
    {
        orbit_inputs.push_back(input);
    }

    //std::string name = "trapp" + planet_name;
    Planet planet(inputs, orbit_inputs);

    std::string output_name = "output.csv";
    std::ofstream output;
    output.open(output_name.c_str());
    output << "time" << ',' << "gamma" << ',' << "g-dot" << ',' << "mean-motion" << ',' << "n-dot" << ',' << "eccentricity" << ',' << "phi" << ',' << "phi_dot" << std::endl;

    double run_time = 500000; //in simulation years

    while(planet.get_time() < run_time)
    {
        planet.solve();
        output << std::setprecision(8) << planet.get_time() << ',' << std::setprecision(6) << planet.get_gamma() << ',' << planet.get_gamma_dot()  << ','
         << planet.mean_motion(planet.get_time()) << ','  << planet.mean_motion_dot(planet.get_time()) << ',' << planet.eccentricity(planet.get_time())
         << ',' << planet.get_phi() << ',' << planet.get_phi_dot() << std::endl;
    }



    return 0;
}
