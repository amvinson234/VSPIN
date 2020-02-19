#include "planet.h"
#include "atmosphere.h"
#include <cmath>
#include "damp.h"

Planet::Planet()
{

}

Planet::Planet(std::vector<double> inputs, std::vector<double> orbit_inputs)
{

    if(inputs.size() != 9) //return error

    //parameters and initial conds
    gamma = inputs[0];
    gamma_dot = inputs[1];
    _mass_star = inputs[2] * MSUN;
    _B_A_C = inputs[3];
    _moi_coeff = inputs[4];
    _tau_M = inputs[5];
    _tau_A = inputs[6];
    _mu = inputs[7] * SEC_PER_YEAR * SEC_PER_YEAR;
    _alpha = inputs[8];

    //features on/off
    damping_on = inputs[9];
    atmosphere_on = inputs[10];
    driving_on = inputs[11];

    //mass and radius params

    _j1 = orbit_inputs[0];
    _j2 = orbit_inputs[1];
    _semi_major = orbit_inputs[2] * AEARTH;
    _semi_major_sec = orbit_inputs[3] * AEARTH;
    _mass = orbit_inputs[4] * MEARTH;
    _mass_sec = orbit_inputs[5] * MEARTH;
    _ecc = orbit_inputs[6];
    _radius = orbit_inputs[7] * REARTH;
    phi = orbit_inputs[8];
    phi_dot = orbit_inputs[9];

    //initialize to Earth values if other planet details not specified
    _mean_motion = 63.0;
    std::cerr << "remember to generalize mean motion later!!!" << std::endl;

    omega_m();

    std::cerr << "omega_m: " <<  _omega_m << std::endl;

    time = 0.; //years
    time_step = 1.; //years

    _min_dt = 2*PI/mean_motion(0)/32.;
    _max_dt = INFINITY;

    atmosphere = new Atmosphere(this);
}


double Planet::omega_m()
{
	int q = _j1 - _j2; //MMR order
    double alpha = std::min(get_semi_major(),get_semi_major_sec()) / std::max(get_semi_major(), get_semi_major_sec());
    double C_r = get_mass_sec() / get_stellar_mass() * alpha * mean_motion(0) * f_d(_j1,q,alpha,eccentricity(0));
    _omega_m = std::sqrt(std::abs(3 * _j2 * _j2 * C_r * mean_motion(0)));
}
double Planet::f_d(int j, double q, double alpha, double ecc)
// integer j = j_1 from M&D99 Disturbing Function Expansion
// q = order of resonant angle
// alpha is the ratio of inner planet semi-major axis to outer
// ecc is eccentricity
{
    if(q == 1)
    {
        if (j == 2)
        {
            return -0.749967 / alpha  * ecc;
        }
        else if (j == 3)
        {
            return -1.54553 / alpha  * ecc;
        }
        else if (j == 4)
        {
            return -2.34472 / alpha  * ecc;
        }
        else if (j == 5)
        {
            return -3.14515 / alpha  * ecc;
        }
        else if (j == 6)
        {
            return -3.94613 / alpha  * ecc;
        }
    }
    else if(q == 2)
    {
        if (j == 3)
        {
            return 0.287852 / alpha * ecc * ecc;
        }
        else if (j == 5)
        {
            return 2.32892 / alpha * ecc * ecc;
        }
        else if (j == 7)
        {
            return 6.28903 / alpha * ecc * ecc;
        }
        else if (j == 9)
        {
            return 12.1673 / alpha * ecc * ecc;
        }
        else if (j == 11)
        {
            return 19.9639 / alpha * ecc * ecc;
        }
    }
    else if (q == 3 || q == 4)
    {
        double s = 1/2.;
        double coeff1 = 2;
        for (int i = 1; i < j - 1; i++)
        {
            coeff1 = coeff1 * (s+i-1.) / i;
        }

        double coeff2 = coeff1 * s * (s+j) / (j+1);

        double coeff3 = coeff1 * s * s * (s+1) * (s+j) * (s+j+1) / (2 * (j+1) * (j+2));

        double A = coeff1 * std::pow(alpha,j) + coeff2 * std::pow(alpha,j+2) + coeff3 * std::pow(alpha,4+j);

        double DA = coeff1 * j * std::pow(alpha,(j-1));
                + coeff2 * (j+2) * std::pow(alpha,(j+1))
                + coeff3 * (4+j) * std::pow(alpha,(3+j));
        double DDA = coeff1 * j * (j-1) * std::pow(alpha,(j-2))
                + coeff2 * (j+2) * (j+1) * std::pow(alpha,(j))
                + coeff3 * (4+j) * (3+j) * std::pow(alpha,(2+j));
        double DDDA = coeff1 * j * (j-1) * (j-2) * std::pow(alpha,(j-3))
                + coeff2 * (j+2) * (j+1) * j * std::pow(alpha,(j-1))
                + coeff3 * (4+j) * (3+j) * (2+j) * std::pow(alpha,(1+j));
        double DDDDA = coeff1 * j * (j-1) * (j-2) * (j-3) * std::pow(alpha,(j-4))
                + coeff2 * (j+2) * (j+1) * j * (j-1) * std::pow(alpha,(j-2))
                + coeff3 * (4+j) * (3+j) * (2+j) * (1+j) * std::pow(alpha,(j));

        if (q == 3)
        {
            double f82 = 1/48. * (A * (-26 * j + 30 * j*j - 8 * j*j*j)
                        + DA * (-9 * alpha + 27 * j*alpha -12 * j*j*alpha)
                        + DDA * (6 * alpha*alpha - 6 * j*alpha*alpha)
                        + DDDA * alpha*alpha*alpha);
            return ecc*ecc*ecc*f82;
        }
        else if (q == 4)
        {
            double f90 = 1/384. * (A * (-206. * j + 283 * std::pow(j,2) -120 * std::pow(j,3) + 16 * std::pow(j,4))
                    + DA * (-64*alpha + 236*alpha*j* - 168*std::pow(j,2) *alpha + 32 * std::pow(j,3) * alpha)
                    + DDA *(48 * std::pow(alpha,2) - 78 *j * std::pow(alpha,2) + 24 * std::pow(j,2) * std::pow(alpha,2))
                    + DDDA * (-12 * std::pow(alpha,3) + 8 * j * std::pow(alpha,3))
                    + DDDDA * std::pow(alpha,4));
            return pow(ecc,4) * f90; //return f_d
        }
    }

}


double Planet::mean_motion(double t)
{
    return _mean_motion;
}
double Planet::eccentricity(double t)
{
    return _ecc;
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
double Planet::get_phi()
{
    return phi;
}
double Planet::get_phi_dot()
{
    return phi_dot;
}

double Planet::get_semi_major()
{
    return _semi_major;
}

double Planet::get_semi_major_sec()
{
    return _semi_major_sec;
}

double Planet::get_ecc()
{
    return _ecc;
}

double Planet::get_radius()
{
    return _radius;
}
double Planet::get_mass()
{
    return _mass;
}

double Planet::get_mass_sec()
{
    return _mass_sec;
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

double Planet::get_damp(double t)
{

    return damp(this,t) / _moi_coeff / _mass / _radius / _radius;
}

double Planet::get_atmospheric_damp()
{
    return atmosphere->damp(this) / _moi_coeff / _mass / _radius / _radius;
}
