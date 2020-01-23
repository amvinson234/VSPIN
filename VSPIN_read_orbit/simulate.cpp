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
            << planet.get_theta() << ',' << planet.get_theta_dot() << ',' << planet.get_mean_anomaly() << ','
            << planet.mean_motion(planet.get_time()) << ',' << planet.mean_motion_dot(planet.get_time()) << ',' << planet.eccentricity(planet.get_time())
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
    planet = Planet(name, M_plan, plan_radius, inputs);

    std::string output_name = "output_" + planet_name + run_name + ".csv";

    output.open(output_name.c_str());
    output << "time" << ',' << "gamma" << ',' << "g-dot" << ',' << "theta" << ',' << "theta-dot" << ',' << "mean-anomaly" << ','
        << "mean-motion" << ',' << "n-dot" << ',' << "eccentricity" << std::endl;

    double last_printed_time = 0.0;
    double dt_sampling = 2*PI / planet.mean_motion(0.0) / 2.0;
}

double Simulate::_run_time = 0.0;
bool Simulate::exiting = false;
std::ofstream Simulate::output;
Planet Simulate::planet;
double Simulate::last_printed_time = 0.0;
double Simulate::dt_sampling = 0.0;
std::string Simulate::run_name = "";
