#include "simulate.h"

void Simulate::Start(std::string run, double run_time)
{
    _run_time = run_time;
    run_name = run;
    setup();

    while(!exiting)
    {
        sim_loop();
    }
    output.close();
}

void Simulate::sim_loop()
{
    planet.solve();
    if(planet.get_time() - last_printed_time >= dt_sampling)
    {
        output << std::setprecision(8) << planet.get_time() << ',' << std::setprecision(6) << planet.get_gamma() << ',' << planet.get_gamma_dot()  << ','
            << planet.mean_motion(planet.get_time()) << ',' << planet.mean_motion_dot(planet.get_time()) << ',' << planet.eccentricity(planet.get_time()) << ','
            << planet.get_damp(planet.get_time()) << ',' << planet.get_atmospheric_damp()
            << std::endl;
        last_printed_time = planet.get_time();
    }
    if (planet.get_time() > _run_time) exiting = true;
}

void Simulate::setup()
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


    planet = Planet(inputs, orbit_inputs);

    std::string output_name = "output.csv";

    output.open(output_name.c_str());
    output << "time" << ',' << "gamma" << ',' << "g-dot" << ','
        << "mean-motion" << ',' << "n-dot" << ',' << "eccentricity"
        << ',' << "tidal_damp" << ',' << "atmospheric_damp" << std::endl;


    dt_sampling = 2*PI / planet.mean_motion(0.0) / 2.0;
}

double Simulate::_run_time = 0.0;
bool Simulate::exiting = false;
std::ofstream Simulate::output;
Planet Simulate::planet;
double Simulate::last_printed_time = 0.0;
double Simulate::dt_sampling = 0.0;
std::string Simulate::run_name = "";
