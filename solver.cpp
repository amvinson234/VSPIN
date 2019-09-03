#include "planet.h"
#include "damp.h"
#include <cmath>

double epsilon = 0.000001; //error criterion

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

        /********************************************************************************
        * Stepsize adjustments
        ********************************************************************************/

        double delta_gamma = gamma - gamma_b;
        double delta_gamma_dot = gamma_dot - gamma_dot_b;

        if(abs(delta_gamma) > abs(epsilon * gamma)) //convergence criterion not reached
        {
            time_step = time_step / 2.0;
            gamma = gamma_init;
            gamma_dot = gamma_dot_init;
        }
        else //convergence criterion reached. updated new time and double next time step.
        {
            time_step *= 2.;
            time += time_step;
        }


    }while(abs(delta_gamma) > abs(epsilon * gamma) && time_step > _min_dt);

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
    double damping = damp(this,_mass_star);
    return -1/2. * omega_s(t) * omega_s(t) * sin(2*x) - mean_motion_dot(t) - damping * x_dot; //x and x_dot represent gamma and gamma_dot, respectively
}

double Planet::mean_motion(double t)
{
    /********************************************************************
     * Need to put stuff in here. Placeholder for now.                  *
     ********************************************************************/
    return 2*PI / 365/24/3600;
}

double Planet::mean_motion_dot(double t)
{
    /********************************************************************
     * Need to put stuff in here. Placeholder for now.                  *
     ********************************************************************/
    return 0;
}

double Planet::omega_s(double t)
{
    return mean_motion(t) * sqrt(3*abs(H(0,ecc) * _B_A_C));
}
