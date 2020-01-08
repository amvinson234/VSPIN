double f_d(int j, double q, double s, double alpha, double ecc)
// integer j = j_1 from M&D99 Disturbing Function Expansion
// q = order of resonant angle
// s is an integer factor of 1/2
// alpha is the ratio of inner planet semi-major axis to outer
// ecc is eccentricity
{
    //First, define A
    double coeff1 = 2
    for (int i = 1; i < j - 1; i++)
    {
        coeff1 = coeff1 * (s+i-1.) / i;
    }

    double coeff2 = coeff1 * s * (s+j) / (j+1);

    double coeff3 = coeff1 * s * s * (s+1) * (s+j) * (s+j+1) / (2 * (j+1) * (j+2));

    double A = coeff1 * std::pow(alpha,j) + coeff2 * std::pow(alpha,j+2) + coeff3 * std::pow(alpha,4+j);

    if(q == 1)
    {

    }
    else if(q == 2)
    {

    }
    else if (q == 3)
    {

    }
    else if (q == 4)
    {
        double DA = coeff1 * j * alpha**(j-1)
                + coeff2 * (j+2) * alpha**(j+1)
                + coeff3 * (4+j) * alpha**(3+j);
        double DDA = coeff1 * j * (j-1) * alpha**(j-2)
                + coeff2 * (j+2) * (j+1) * alpha**(j)
                + coeff3 * (4+j) * (3+j) * alpha**(2+j);
        double DDDA = coeff1 * j * (j-1) * (j-2) * alpha**(j-3)
                + coeff2 * (j+2) * (j+1) * j * alpha**(j-1)
                + coeff3 * (4+j) * (3+j) * (2+j) * alpha**(1+j);
        double DDDDA = coeff1 * j * (j-1) * (j-2) * (j-3) * alpha**(j-4)
                + coeff2 * (j+2) * (j+1) * j * (j-1) * alpha**(j-2)
                + coeff3 * (4+j) * (3+j) * (2+j) * (1+j) * alpha**(j);

        double f90 = 1/384. * (A * (-206. * j + 283 * j**2 -120 * j**3 + 16 * j**4)
                + DA * (-64*alpha + 236*alpha*j* - 168*j**2 *alpha + 32 * j**3 * alpha)
                + DDA *(48 * alpha**2 - 78 *j * alpha**2 + 24 * j**2 * alpha**2)
                + DDDA * (-12 *alpha**3 + 8 * j * alpha**3)
                + DDDDA * alpha**4);
        return pow(ecc,4) * f90; //return f_d
    }
    else
    {
        //throw an error here
    }
}
