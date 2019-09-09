#include "planet.h"
#include <cmath>


Planet::Planet(std::string name, double mass, double radius, std::vector<double> inputs)
{

    if(inputs.size() != 9) //return error

    gamma = inputs[0];
    gamma_dot = inputs[1];
    _mass_star = inputs[2];
    _B_A_C = inputs[3];
    _moi_coeff = inputs[4];
    _tau_M = inputs[5];
    _tau_A = inputs[6];
    _mu = inputs[7];
    _alpha = inputs[8];

    _mass = mass;
    _radius = radius;

    read_orbit("tamayo_runs/tamayo_grimm_" + name + ".txt");

    //initialize to Earth values if other planet details not specified
    _semi_major = std::pow(4 * PI * PI / std::pow(mean_motion(0),2) * _mass_star / MSUN, 1.0/3.0) * AEARTH;

    time = 0.; //years
    time_step = 1.; //years

    _min_dt = 2*PI/mean_motion(0);
    _max_dt = INFINITY;

}

std::string Planet::name()
{
    return _name;
}

double Planet::get_stellar_mass()
{
    return _mass_star;
}
double Planet::get_gamma()
{
    return gamma;
}
double Planet::get_gamma_dot()
{
    return gamma_dot;
}

double Planet::get_semi_major()
{
    return _semi_major;
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

void Planet::read_orbit(std::string input_file_name)
{
    std::ifstream input(input_file_name.c_str());
    std::string header_trash;
    std::getline(input, header_trash);
    int begin_line = sim_loop(input_file_name);

    for(int ctr = 0; ctr < begin_line - 1; ctr++)
    {
        std::string line;
        std::getline(input,line);
    }

    double t, mm, e, p;
    std::vector<double> input_time, input_mm, input_ecc;//, input_peri;
    while(input >> t)
    {
        input >> mm >> e;// >> p;
        input_time.push_back(t);
        input_mm.push_back(mm);
        input_ecc.push_back(e);
        //input_peri.push_back(p);
    }

    double initial_time = input_time[0];
    for(int i = 0; i < input_time.size(); i++) input_time[i] = input_time[i] - initial_time; //set initial time to 0.

    //Following lines for for smooth stitch back to front. ensures that back of vectors equal its front.
    /*
    input_time.push_back(t+input_time[1]-input_time[0]);
    input_mm.push_back(input_mm[1]);
    input_ecc.push_back(input_ecc[1]);
    input_peri.push_back(input_peri[1]);
    */

    input.close();

    spline_time = input_time;
    spline_a_mm = input_mm;
    spline_a_e = input_ecc;
    spline_delta_t = (spline_time[spline_time.size()-1] - spline_time[0]) / spline_time.size();

    spline(spline_time, spline_a_mm, spline_b_mm, spline_c_mm, spline_d_mm);
    spline(spline_time, spline_a_e, spline_b_e, spline_c_e, spline_d_e);

}
