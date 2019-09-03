#ifndef PLANET_H_INCLUDED
#define PLANET_H_INCLUDED

#include <string>
#include "constants.h"

class Planet
{
public:
    Planet(std::string name, double radius, double mass, double mu, double alpha, double tau_M, double tau_A);

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

    double mean_motion(double t);
    double mean_motion_dot(double t);

    void solve(); //updates gamma and gamma_dot by performing integration. also calls update_orbit()

private:
    std::string _name;
    double  _radius;
    double  _mass;
    double  _mu;
    double  _alpha;
    double  _tau_M;
    double  _tau_A;
    double  _min_dt;
    double  _max_dt;
    double _mean_motion; //

    double semi_major;
    double ecc;

    double gamma;
    double gamma_dot;

    double time;
    double time_step;

    double func_gdd(double time, double x, double x_dot, double n_dot);

    std::pair<double,double> integrate(double t, double h, double y, double y_dot);
    void update_orbit();
    void update_time_step();
};

#endif // PLANET_H_INCLUDED
