#ifndef SPLINE_H_INCLUDED
#define SPLINE_H_INCLUDED

#include <vector>
#include <cmath>
#include <string>
#include <fstream>


struct Spline
{
    Spline();
    Spline(std::vector<double> x_input, std::vector<double> y_input);

    double spline_interpolate(double x);
    double spline_interpolate(double x, int index);
    double spline_interpolate_derivative(double x);
    double spline_interpolate_derivative(double x, int index);

    std::vector<double> ddy; //double derivative of y, to be calculated
    std::vector<double> x_data;
    std::vector<double> y_data;
};

#endif // SPLINE_H_INCLUDED
