#include "spline.h"
#include <iostream>

Spline::Spline()
{

}

Spline::Spline(std::vector<double> x_dat, std::vector<double> y_dat)
{
    y_data = y_dat;
    x_data = x_dat;

    int n = x_dat.size();
    if(n != y_dat.size()) std::cerr << "number of x data points does not match number of y data points for spline interpolation." << std::endl;

    //initialize ddy
    for(int i = 0; i < n; i++) ddy.push_back(0.0);

    std::vector<double> u;
    for(int i = 0; i < n-1; i++) u.push_back(0.0);

    for(int i = 1; i < n-1; i++)
    {
        double sig = (x_dat[i] - x_dat[i-1]) / (x_dat[i+1] - x_dat[i-1]);
        double p = sig * ddy[i-1] + 2.0;
        ddy[i] = (sig-1.0) / p;
        u[i] = (y_dat[i+1] - y_dat[i]) / (x_dat[i+1] - x_dat[i]) - (y_dat[i] - y_dat[i-1]) / (x_dat[i] - x_dat[i-1]);
        u[i] = (6.0 * u[i] / (x_dat[i+1] - x_dat[i-1]) - sig * u[i-1]) / p;
    }

    for (int i = n-2; i >= 0; i--)
    {
        ddy[i] = ddy[i] * ddy[i+1] + u[i];
    }
}

double Spline::spline_interpolate(double x)
{
    int index = 0;
    while(x < x_data[index])
    {
        index++;
    }

    double A = (x_data[index + 1] - x) / (x_data[index + 1] - x_data[index]);
    double B = 1.0 - A;
    double C = 1/6.0 * (A*A*A - A) * std::pow(x_data[index+1] - x_data[index],2);
    double D = 1/6.0 * (B*B*B - B) * std::pow(x_data[index+1] - x_data[index],2);


    return A * y_data[index] + B * y_data[index+1] + C*ddy[index] + D*ddy[index+1];
}

double Spline::spline_interpolate(double x, int index)
{

    double A = (x_data[index + 1] - x) / (x_data[index + 1] - x_data[index]);
    double B = 1.0 - A;
    double C = 1/6.0 * (A*A*A - A) * std::pow(x_data[index+1] - x_data[index],2);
    double D = 1/6.0 * (B*B*B - B) * std::pow(x_data[index+1] - x_data[index],2);


    return A * y_data[index] + B * y_data[index+1] + C*ddy[index] + D*ddy[index+1];
}

double Spline::spline_interpolate_derivative(double x)
{
    int index = 0;
    while(x < x_data[index])
    {
        index++;
    }

    double A = (x_data[index + 1] - x) / (x_data[index + 1] - x_data[index]);
    double B = 1.0 - A;

    double dA = 1.0 / (x_data[index] - x_data[index+1]);
    double dB = - dA;
    double dC = 1/6.0 * std::pow(x_data[index+1] - x_data[index],2) * (3.0*A*A - 1.0) * dA;
    double dD = 1/6.0 * std::pow(x_data[index+1] - x_data[index],2) * (3.0*B*B - 1.0) * dB;

    return dA * y_data[index] + dB * y_data[index+1] + dC*ddy[index] + dD*ddy[index+1];
}

double Spline::spline_interpolate_derivative(double x, int index)
{
    double A = (x_data[index + 1] - x) / (x_data[index + 1] - x_data[index]);
    double B = 1.0 - A;

    double dA = 1.0 / (x_data[index] - x_data[index+1]);
    double dB = - dA;
    double dC = 1/6.0 * std::pow(x_data[index+1] - x_data[index],2) * (3.0*A*A - 1.0) * dA;
    double dD = 1/6.0 * std::pow(x_data[index+1] - x_data[index],2) * (3.0*B*B - 1.0) * dB;

    return dA * y_data[index] + dB * y_data[index+1] + dC*ddy[index] + dD*ddy[index+1];
}



