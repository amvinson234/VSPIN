#ifndef PLANET_H_INCLUDED
#define PLANET_H_INCLUDED

#include <iostream>
#include <string>
#include "constants.h"

class Planet
{
public:
    //possibly collapse all params into a vector later
    Planet(std::string name, double radius, double mass, double B_A_C, double mu, double alpha, double tau_M, double tau_A, double mass_star);

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

    double mean_motion(double t); //returns mean motion at arbitrary time
    double mean_motion_dot(double t); //returns time derivative of mean motion at arbitrary time
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
    double _mean_motion;
    double _B_A_C; //triaxiality
    double _mass_star;

    double  _min_dt;
    double  _max_dt;


    double semi_major;
    double ecc; //may want to change this to a function, ecc(t).

    double gamma;
    double gamma_dot;

    double time;
    double time_step;

    double func_gdd(double t, double x, double x_dot);

    std::pair<double,double> integrate(double t, double h, double y, double y_dot); //integrates forward from time, t, to one timestep, t + h. returns pair: first = gamma, second = gamma_dot
    void update_orbit();
};

#endif // PLANET_H_INCLUDED
