#ifndef SPLINE_H_INCLUDED
#define SPLINE_H_INCLUDED

#include <vector>
#include <cmath>
#include <string>
#include <fstream>

/*
void spline(std::vector<double> &time, std::vector<double> &data, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d);
double spline_inter(double a, double b, double c, double d, double t, double ti);
double spline_inter_deriv(double b, double c, double d, double t, double ti);

int sim_loop(std::string file_name);
*/

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
