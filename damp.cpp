#include "damp.h"
#include "constants.h"
#include <cmath>
#include <iostream>

//"simplified" Efroimsky Torque
//Taken from Makarov 2012, eqn 10
double damp(Planet *planet)
//double damp(double M_star, double plan_radius, double semi_major, double ecc, double M_plan, double mean_motion, double gdot, double mu, double alpha, double tau_M, double tau_A)
{
    double mean_motion = planet->get_mean_motion();
    double gdot = planet->get_gamma_dot();
    double ecc = planet->get_ecc();
    double M_star = planet->get_stellar_mass();

    double sum = 0;
    for(int q = -1; q < 5; q++)
    {
        double omega_220q = (2+q)*mean_motion - 2*(gdot + mean_motion);
        int sign_omega_220q = omega_220q / std::abs(omega_220q);
        double chi_220q = std::abs(omega_220q);
        double G_20q = H(q, ecc);
        //Need to give values for mu, alpha, taus...
        //perhaps taken from a param file, and/or as an input in the damp() function
        //sum += pow(G_20q,2) * K_c(chi_220q,plan_radius,M_plan,mu,alpha,tau_M,tau_A) * sign_omega_220q;
        sum += std::pow(G_20q,2) * K_c(planet,chi_220q) * sign_omega_220q;
    }

    double plan_radius = planet->get_radius();
    double semi_major = planet->get_semi_major();

    double sec_per_year = 365. * 24. * 3600.; //for conversion of G to time units of years

    return (3/2.0 * G * sec_per_year * sec_per_year * M_star * M_star * std::pow(plan_radius,5) / std::pow(semi_major,6)) * sum;
    //consider making function just return the sum part, so don't have to keep re-calculating the constant coefficient in front.
}

double H(int q, double ecc)
{
    if(q == -1)
    {
        return -1.0/2.0 * ecc + 1.0/16.0 * std::pow(ecc, 3);
    }
    else if(q == 0)
    {
        return 1.0 - 5/2.0 * std::pow(ecc,2) + 13.0/16.0 * std::pow(ecc,4);
    }
    else if(q == 1)
    {
        return 7.0/2.0 * ecc - 123.0/16.0 * std::pow(ecc,3);
    }
    else if(q == 2)
    {
        return 17.0/2.0 * std::pow(ecc,2) - 115.0/6.0 * std::pow(ecc,4);
    }
    else if (q == 3)
    {
        return 845.0 / 48.0 * std::pow(ecc,3);
    }
    else if(q == 4)
    {
        return 533.0 / 16.0 * std::pow(ecc,4);
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

    double lambda_2 = 4 * PI * (2 * 2 * 2 + 4 * 2 + 3) * std::pow(radius, 4) * mu / ( 3 * 2 * G * M_plan *M_plan );
    double real_part = chi + std::pow(chi,1 - alpha) * std::pow(tau_A, -alpha) * std::cos(alpha * PI / 2) * std::tgamma(1+alpha);
    double imaginary_part = - std::pow(tau_M,-1) - std::pow(chi,1-alpha) * std::pow(tau_A,-alpha) * std::sin(alpha*PI/2) * std::tgamma(1+alpha);

    double K = -3 / 2.0 * lambda_2 * chi * imaginary_part / (std::pow(real_part + lambda_2 * chi, 2) + std::pow(imaginary_part,2) );
    return K;
}


