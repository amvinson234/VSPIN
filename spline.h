#ifndef SPLINE_H_INCLUDED
#define SPLINE_H_INCLUDED

#include <vector>
#include <cmath>
#include <string>
#include <fstream>

void spline(std::vector<double> &time, std::vector<double> &data, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d);
double spline_inter(double a, double b, double c, double d, double t, double ti);
double spline_inter_deriv(double b, double c, double d, double t, double ti);

int sim_loop(std::string file_name);

#endif // SPLINE_H_INCLUDED
