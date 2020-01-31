#include "spline.h"
#include <iostream>
#include <float.h>
#include <cmath>

Spline::Spline()
{

}

Spline::Spline(std::vector<double> x_dat, std::vector<double> y_dat)
{
    y_data = y_dat;
    x_data = x_dat;

    delta_t = (x_data[x_data.size()-1] - x_data[0]) / x_data.size();

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

    end_index = x_data.size() - 1; //initialize to end of data
    begin_end_index(); //set new end index
}

void Spline::begin_end_index()
{
    int test_index = y_data.size() - y_data.size() / 4; //start a quarter of the way towards the end

    double y_begin = y_data[0];
    double y_deriv_begin = spline_interpolate_derivative(x_data[0],0);

    double score = INFINITY;
    int best_index = test_index;

    for(test_index; test_index < y_data.size(); test_index++)
    {
        //difference between values at beginning and at current test index
        double delta_y = y_data[test_index] - y_begin;
        double delta_y_deriv = spline_interpolate_derivative(x_data[test_index],test_index) - y_deriv_begin;

        double test_score = std::sqrt(delta_y * delta_y / std::pow(y_begin + DBL_EPSILON, 2)
                        + delta_y_deriv * delta_y_deriv / std::pow(y_deriv_begin + DBL_EPSILON, 2));
        if(test_score < score)
        {
            score = test_score;
            best_index = test_index;
        }
    }

}

double Spline::spline_interpolate(double x)
{
    //int index = 0;
    x = std::fmod(x - x_data[0], x_data[end_index] - x_data[0]) + x_data[0];

    int index = int((x-x_data[0]) / delta_t) - 2;
    if(index < 0) index = 0;

    while(x >= x_data[index])
    {
        index++;
    }
    index--;

    double A = (x_data[index + 1] - x) / (x_data[index + 1] - x_data[index]);
    double B = 1.0 - A;
    double C = 1/6.0 * (A*A*A - A) * std::pow(x_data[index+1] - x_data[index],2);
    double D = 1/6.0 * (B*B*B - B) * std::pow(x_data[index+1] - x_data[index],2);


    return A * y_data[index] + B * y_data[index+1] + C*ddy[index] + D*ddy[index+1];
}

double Spline::spline_interpolate(double x, int index)
{

    if(index >= end_index)
    {
        index = index % end_index;
        x = std::fmod(x - x_data[0], x_data[end_index] - x_data[0]) + x_data[0];
    }

    double A = (x_data[index + 1] - x) / (x_data[index + 1] - x_data[index]);
    double B = 1.0 - A;
    double C = 1/6.0 * (A*A*A - A) * std::pow(x_data[index+1] - x_data[index],2);
    double D = 1/6.0 * (B*B*B - B) * std::pow(x_data[index+1] - x_data[index],2);


    return A * y_data[index] + B * y_data[index+1] + C*ddy[index] + D*ddy[index+1];
}

double Spline::spline_interpolate_derivative(double x)
{
   // int index = 0;
    x = std::fmod(x - x_data[0], x_data[end_index] - x_data[0]) + x_data[0];

    //leapfrog index close to, but less than (or equal) to actual value, so that while loop doesn't slow us down

    int index = int((x-x_data[0]) / delta_t) - 2;
    if(index < 0) index = 0;


    while(x >= x_data[index])
    {
        index++;
    }
    index--;

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

    if(index >= end_index)
    {
        index = index % end_index;
        x = std::fmod(x - x_data[0], x_data[end_index] - x_data[0]) + x_data[0];
    }

    double A = (x_data[index + 1] - x) / (x_data[index + 1] - x_data[index]);
    double B = 1.0 - A;

    double dA = 1.0 / (x_data[index] - x_data[index+1]);
    double dB = - dA;
    double dC = 1/6.0 * std::pow(x_data[index+1] - x_data[index],2) * (3.0*A*A - 1.0) * dA;
    double dD = 1/6.0 * std::pow(x_data[index+1] - x_data[index],2) * (3.0*B*B - 1.0) * dB;

    return dA * y_data[index] + dB * y_data[index+1] + dC*ddy[index] + dD*ddy[index+1];
}



