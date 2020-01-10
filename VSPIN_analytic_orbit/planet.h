#ifndef PLANET_H_INCLUDED
#define PLANET_H_INCLUDED

#include <iostream>
#include <string>
#include "constants.h"
#include <fstream>
#include <vector>

class Atmosphere;

class Planet
{
public:
    Planet();
    Planet(std::vector<double> inputs, std::vector<double> orbit_inputs);

    double get_gamma();
    double get_gamma_dot();
    double get_phi();
    double get_phi_dot();

    double get_ecc();
    double get_mean_motion();
    double mean_motion(double t);
    double eccentricity(double t);
    double mean_motion_dot(double t);
    double get_semi_major();

    double get_radius();
    double get_mass();
    double get_mu();
    double get_alpha();
    double get_tau_M();
    double get_tau_A();
    double get_time();
    double get_stellar_mass();

    double get_semi_major_sec();
    double get_mass_sec();

    double omega_s(double t); //returns natural libration frequency, omega_s, at arbitrary time.
    double omega_m();

    void solve(); //updates gamma and gamma_dot by performing integration. goes forward, adjusting time_step until convergence crit is reached, then stops.

private:
    double  _radius;
    double  _mass;
    double _mass_sec;
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
    //double _epsilon = 1.0e-6; //convergence criterion.


    double _semi_major;
    double _semi_major_sec;
    double _mean_motion;
    double _ecc;

    int _j1;
    int _j2;

   // double ecc; //may want to change this to a function, ecc(t).

    double gamma;
    double gamma_dot;

    double phi;
    double phi_dot;

    double time;
    double time_step;

    double _omega_m;
    double _omega_s;

    double f_d(int j, double q, double alpha, double ecc);

    double func_gdd(double t, double x, double x_dot);
    double func_phi(double t, double x, double x_dot);

    std::pair<double,double> integrate(double t, double h, double y, double y_dot); //integrates forward from time, t, to one timestep, t + h. returns pair: first = gamma, second = gamma_dot
    std::pair<double,double> integrate_phi(double t, double h, double y, double y_dot);

    Atmosphere *atmosphere;

    bool driving_on;
    bool atmosphere_on;
    bool damping_on;

};

#endif // PLANET_H_INCLUDED
