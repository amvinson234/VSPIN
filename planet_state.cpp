#include "planet.h"
#include <cmath>


Planet::Planet(std::string name, double radius, double mass, double B_A_C, double mu, double alpha, double tau_M, double tau_A, double mass_star)
{
    _name = name;
    _radius = radius;
    _mass = mass;
    _B_A_C = B_A_C;
    _mu = mu;
    _alpha = alpha;
    _tau_M = tau_M;
    _tau_A = tau_A;
    _mass_star = mass_star;
    _moi_coeff = 2.0/5.0;

    gamma = 0.;
    gamma_dot = 2 * PI; //rad per year

    /*************
    Read in orbit
    *************/

    read_orbit(name + "_avg_longg_multisim13.txt");


    //initialize to Earth values if other planet details not specified
    semi_major = AEARTH;
    _mean_motion = 2 * PI; //rad per year
    ecc = 0.0167;

    time = 0.; //years
    time_step = 1.; //years

    _min_dt = 0.;
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

double Planet::get_ecc()
{
    return ecc;
}
double Planet::get_mean_motion()
{
    return _mean_motion;
}
double Planet::get_semi_major()
{
    return semi_major;
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
    std::vector<double> input_time, input_mm, input_ecc, input_peri;
    while(input >> t)
    {
        input >> mm >> e >> p;
        input_time.push_back(t);
        input_mm.push_back(mm);
        input_ecc.push_back(e);
        input_peri.push_back(p);
    }

    for(int i = 0; i < input_time.size(); i++) input_time[i] = input_time[i] - input_time[0]; //set initial time to 0.

    //Following lines for for smooth stitch back to front. ensures that back of vectors equal its front.

    input_time.push_back(t+input_time[1]-input_time[0]);
    input_mm.push_back(input_mm[1]);
    input_ecc.push_back(input_ecc[1]);
    input_peri.push_back(input_peri[1]);

    input.close();

    spline_time = input_time;
    spline_mm = input_mm;
    spline_delta_t = (spline_time[spline_time.size()-1] - spline_time[0]) / spline_time.size();

    spline(spline_time, spline_mm, spline_b, spline_c, spline_d);

}
