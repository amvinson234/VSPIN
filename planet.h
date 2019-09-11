#ifndef PLANET_H_INCLUDED
#define PLANET_H_INCLUDED

#include <iostream>
#include <string>
#include "constants.h"
#include "spline.h"
#include <fstream>

class Atmosphere;

class Planet
{
public:
    //possibly collapse all params into a vector later
    Planet();
    Planet(std::string name,  double mass, double radius, std::vector<double> inputs);
    std::string name();
    double get_gamma();
    double get_gamma_dot();

    double get_ecc();
    double get_mean_motion();
    double get_semi_major();

    double get_radius();
    double get_mass();
    double get_mu();
    double get_alpha();
    double get_tau_M();
    double get_tau_A();
    double get_time();
    double get_stellar_mass();

    double mean_motion(double t); //returns mean motion at arbitrary time
    double mean_motion_dot(double t); //returns time derivative of mean motion at arbitrary time
    double eccentricity(double t); //returns eccentricity at arbitrary time
    double omega_s(double t); //returns natural libration frequency, omega_s, at arbitrary time.

    void solve(); //updates gamma and gamma_dot by performing integration. goes forward, adjusting time_step until convergence crit is reached, then stops.

private:
    std::string _name;
    double  _radius;
    double  _mass;
    double  _mu;
    double  _alpha;
    double  _tau_M;
    double  _tau_A;
   // double _mean_motion;
    double _B_A_C; //triaxiality
    double _mass_star;
    double _moi_coeff; //moment of inertia coefficient

    double  _min_dt;
    double  _max_dt;
    double _epsilon = 1.0e-6; //convergence criterion.


    double _semi_major;
   // double ecc; //may want to change this to a function, ecc(t).

    double gamma;
    double gamma_dot;

    double time;
    double time_step;

    //Following are spline data, calculated in read_orbit, needed to calculate mean_motion and eccentricity at arbitrary times
    std::vector<double> spline_b_mm;
    std::vector<double> spline_c_mm;
    std::vector<double> spline_d_mm;
    std::vector<double> spline_a_mm;
    std::vector<double> spline_b_e;
    std::vector<double> spline_c_e;
    std::vector<double> spline_d_e;
    std::vector<double> spline_a_e;
    std::vector<double> spline_time;
    double spline_delta_t;


    double func_gdd(double t, double x, double x_dot);

    std::pair<double,double> integrate(double t, double h, double y, double y_dot); //integrates forward from time, t, to one timestep, t + h. returns pair: first = gamma, second = gamma_dot
    void read_orbit(std::string file_name);

    Atmosphere *atmosphere;

};

#endif // PLANET_H_INCLUDED
