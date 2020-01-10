double f_d(int j, double q, double alpha, double ecc)
// integer j = j_1 from M&D99 Disturbing Function Expansion
// q = order of resonant angle
// alpha is the ratio of inner planet semi-major axis to outer
// ecc is eccentricity
{
    //First, define A
    double s = 1/2.;
    double coeff1 = 2;
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

double omega_m(double mean_motion, double eccentricity, double a, double a_sec, double m_sec, double m_star, int j, int j2)
{
	int q = j - j2; //MMR order
    double alpha = std::min(a,a_sec) / std::max(a, a_sec);
    double C_r = m_sec / m_star * alpha * mean_motion * f_d(j,q,alpha,eccentricity);
    _omega_m = std::sqrt(std::abs(3 * j2 * j2 * C_r * mean_motion * eccentricity));
}

std::pair<double,double> integrate_phi(double t, double h, double y, double y_dot)
{
    double k1,k2,k3,k4,l1,l2,l3,l4;

    k1 = y_dot;
    l1 = func_phi(t, y, y_dot);
    k2 = y_dot + h / 2. * l1;
    l2 = func_phi(t + h/2., y + h/2. * k1 , y_dot + h/2. * l1);
    k3 = y_dot + h / 2. * l2;
    l3 = func_phi(t + h/2., y + h/2. * k2 , y_dot + h/2. * l2);
    k4 = y_dot + h * l3;
    l4 = func_phi(t + h, y + h * k3 , y_dot + h * l3);

    y = y + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
    y_dot = y_dot + h / 6. * (l1 + 2. * l2 + 2. * l3 + l4);

    return std::make_pair(y,y_dot);
}

double func_phi(double t, double x, double x_dot)
{
    return -1/j * _omega_m * _omega_m * sin(x); //x and x_dot represent gamma and gamma_dot, respectively
}
