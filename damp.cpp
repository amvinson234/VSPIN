#include "damp.h"
#include "constants.h"
#include <cmath>
#include <iostream>

//"simplified" Efroimsky Torque
//Taken from Makarov 2012, eqn 10
double damp(Planet *planet, double M_star)
//double damp(double M_star, double plan_radius, double semi_major, double ecc, double M_plan, double mean_motion, double gdot, double mu, double alpha, double tau_M, double tau_A)
{
    double mean_motion = planet->get_mean_motion();
    double gdot = planet->get_gamma_dot();
    double ecc = planet->get_ecc();

    double sum = 0;
    for(int q = -1; q < 5; q++)
    {
        double omega_220q = (2+q)*mean_motion - 2*(gdot + mean_motion);
        int sign_omega_220q = omega_220q / abs(omega_220q);
        double chi_220q = abs(omega_220q);
        double G_20q = H(q, ecc);
        //Need to give values for mu, alpha, taus...
        //perhaps taken from a param file, and/or as an input in the damp() function
        //sum += pow(G_20q,2) * K_c(chi_220q,plan_radius,M_plan,mu,alpha,tau_M,tau_A) * sign_omega_220q;
        sum += pow(G_20q,2) * K_c(planet,chi_220q) * sign_omega_220q;
    }

    double plan_radius = planet->get_radius();
    double semi_major = planet->get_semi_major();

    return (3/2.0 * G * M_star * M_star * pow(plan_radius,5) / pow(semi_major,6)) * sum;
    //consider making function just return the sum part, so don't have to keep re-calculating the constant coefficient in front.
}

double H(int q, double ecc)
{
    if(q == -1)
    {
        return -1/2 * ecc + 1/16.0 * pow(ecc, 3);
    }
    else if(q == 0)
    {
        return 1 - 5/2.0 * pow(ecc,2) + 13/16.0 * pow(ecc,4);
    }
    else if(q == 1)
    {
        return 7/2.0 * ecc - 123/16.0 * pow(ecc,3);
    }
    else if(q == 2)
    {
        return 17/2.0 * pow(ecc,2) - 115/6.0 * pow(ecc,4);
    }
    else if (q == 3)
    {
        return 845 / 48.0 * pow(ecc,3);
    }
    else if(q == 4)
    {
        return 533 / 16.0 * pow(ecc,4);
    }

    //else error...
}


//double K_c(double chi, double radius, double M_plan, double mu, double alpha, double tau_M, double tau_A)
//mu = "unrelaxed rigidity modulus" = 0.8 \times 10^11 kg m^-1 s^-2 for Mercury
//alpha = "tidal responsivity" = 0.2 for Mercury
//tau_M = "maxwell time" = 500 yrs for Mercury
//tau_A = some other freq-dependent time, but can be set equal to tau_M for simplicity, or:
//      tau_A = 50000 exp(-chi/0.2) + 500
//This is all taken from Makarov 2012, eqns 6 - 9

double K_c(Planet *planet, double chi)
{
    double radius = planet->get_radius();
    double M_plan = planet->get_mass();
    double alpha = planet->get_alpha();
    double tau_A = planet->get_tau_A();
    double tau_M = planet->get_tau_M();
    double mu = planet->get_mu();

    double lambda_2 = 4 * PI * (2 * 2 * 2 + 4 * 2 + 3) * pow(radius, 4) * mu / ( 3 * 2 * G * M_plan *M_plan );
    double real_part = chi + pow(chi,1 - alpha) * pow(tau_A, -alpha) * cos(alpha * PI / 2) * tgamma(1+alpha);
    double imaginary_part = -pow(tau_M,-1) - pow(chi,1-alpha) * pow(tau_A,-alpha) * sin(alpha*PI/2) * tgamma(1+alpha);

    double K = -3 / 2.0 * lambda_2 * chi * imaginary_part / (pow(real_part + lambda_2 * chi, 2) + pow(imaginary_part,2) );
    return K;
}


