double f_d(int j, double q, double s, double alpha)
// integer j = j_1 from M&D99 Disturbing Function Expansion
// q = order of resonant angle
// s is an integer factor of 1/2
// alpha is the ratio of inner planet semi-major axis to outer
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

    }
    else
    {
        //throw an error here
    }
}
