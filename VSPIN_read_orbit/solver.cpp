#include "planet.h"
#include "damp.h"
#include "spline.h"
#include "atmosphere.h"
#include <cmath>


//Runge-Kutta integration
std::pair<double,double> Planet::integrate(double t, double h, double y, double y_dot)
{
    double k1,k2,k3,k4,l1,l2,l3,l4;

    k1 = y_dot;
    l1 = func_gdd(t, y, y_dot);
    k2 = y_dot + h / 2. * l1;
    l2 = func_gdd(t + h/2., y + h/2. * k1 , y_dot + h/2. * l1);
    k3 = y_dot + h / 2. * l2;
    l3 = func_gdd(t + h/2., y + h/2. * k2 , y_dot + h/2. * l2);
    k4 = y_dot + h * l3;
    l4 = func_gdd(t + h, y + h * k3 , y_dot + h * l3);

    y = y + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
    y_dot = y_dot + h / 6. * (l1 + 2. * l2 + 2. * l3 + l4);

    return std::make_pair(y,y_dot);
}


void Planet::solve()
{
    double gamma_init = gamma;
    double gamma_dot_init = gamma_dot;

    double delta_gamma;
    double delta_gamma_dot;

    do
    {
        std::pair<double,double> candidate_b = integrate(time,time_step * 2,gamma,gamma_dot);
        double gamma_b = candidate_b.first;
        double gamma_dot_b = candidate_b.second;

        std::pair<double,double> candidate_s = integrate(time,time_step,gamma,gamma_dot);
        gamma = candidate_s.first;
        gamma_dot = candidate_s.second;
        candidate_s = integrate(time+time_step,time_step,gamma,gamma_dot);
        gamma = candidate_s.first;
        gamma_dot = candidate_s.second;

        /*******************************************************************************
         * Stepsize adjustments                                                        *
         *******************************************************************************/

        double delta_gamma = gamma - gamma_b;
        double delta_gamma_dot = gamma_dot - gamma_dot_b;

        if(std::abs(delta_gamma) > std::abs((1.0e-6) * gamma) && time_step > _min_dt) //convergence criterion not reached
        {
            time_step = time_step / 2.0;
            gamma = gamma_init;
            gamma_dot = gamma_dot_init;
        }
        else //convergence criterion reached. updated new time and double next time step.
        {
            time_step *= 2.;
            time += time_step;
            break; //break out of while loop once convergence is reached
        }

    }while(1);

    //Refine estimates:
    gamma = gamma + delta_gamma / 15.;
    gamma_dot = gamma_dot + delta_gamma_dot / 15.;

    //Set radians to [-pi,pi]
    int num_g = gamma/2/PI;
    gamma = gamma - num_g*2*PI;
    if(gamma < 0) gamma = 2*PI+gamma;
    if(gamma > PI) gamma = gamma - 2*PI;
}

//functional form of double derivative of gamma.
double Planet::func_gdd(double t, double x, double x_dot)
{

    double atmosphere_damp;
    double driving;
    double damping;

    if(atmosphere_on) atmosphere_damp = atmosphere->damp(this); else atmosphere_damp = 0;
    if(driving_on) driving = mean_motion_dot(t); else driving = 0;
    if(damping_on) damping = damp(this,t); else damping = 0;

    return -1/2. * omega_s(t) * omega_s(t) * std::sin(2*x) - driving + (damping + atmosphere_damp) / _moi_coeff / _mass / _radius / _radius; //x and x_dot represent gamma and gamma_dot, respectively
}

double Planet::mean_motion(double t)
{
    /********************************************************************
     * Following does not yet support looping over input orbital sims   *
     ********************************************************************/
    int i = int(t / spline_delta_t);
    if(driving_on) return spline_mm.spline_interpolate(t,i);
    else return spline_mm.y_data[0];

}

double Planet::mean_motion_dot(double t)
{
    /********************************************************************
     * Following does not yet support looping over input orbital sims   *
     ********************************************************************/
    int i = int(t / spline_delta_t);
    if(driving_on) return spline_mm.spline_interpolate_derivative(t,i);
    else return 0.0;

    //return 0;
}

double Planet::eccentricity(double t)
{
    /********************************************************************
     * Following does not yet support looping over input orbital sims   *
     ********************************************************************/
    int i = int(t / spline_delta_t);
    if(driving_on) return spline_ecc.spline_interpolate(t,i);
    else return spline_ecc.y_data[0];

}

double Planet::omega_s(double t)
{
    return mean_motion(t) * std::sqrt(3*std::abs(H(0,eccentricity(t)) * _B_A_C));
}
