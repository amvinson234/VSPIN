#include "spline.h"
using namespace std;

void spline(vector<double> &time, vector<double> &data, vector<double> &b, vector<double> &c, vector<double> &d)
{
    vector<double> h, alpha, l, mu, z;
    for (int i = 0; i < time.size()-1; i++)
    {
        h.push_back(time[i+1]-time[i]);
        if(i>0)
        alpha.push_back(3/h[i]*(data[i+1]-data[i])-3/h[i-1]*(data[i]-data[i-1]));
    }
    l.push_back(1.0);
    mu.push_back(0.0);
    z.push_back(0.0);
    for (int i = 1; i < time.size()-1; i++)
    {
        l.push_back(2*(time[i+1]-time[i-1])-h[i-1]*mu[i-1]);
        mu.push_back(h[i]/l[i]);
        z.push_back((alpha[i]-h[i-1]*z[i-1])/l[i]);
    }
    l.push_back(1.0);
    z.push_back(0.0);
    mu.push_back(0.0);

    c.resize(time.size()-1);
    b.resize(time.size()-1);
    d.resize(time.size()-1);

    c[time.size()-2] = 0.0;

    for (int j=time.size()-3;j>-1;j--)
    {
        c[j]=z[j]-mu[j]*c[j+1];
        b[j]=(data[j+1]-data[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3.;
        d[j] = (c[j+1]-c[j])/3./h[j];
    }
}

double spline_inter(double a, double b, double c, double d, double t, double ti)
{
    return a + b*(t-ti) + c*(t-ti)*(t-ti) + d*(t-ti)*(t-ti)*(t-ti);
}

double spline_inter_deriv(double b, double c, double d, double t, double ti)
{
    return b + 2*c*(t-ti) + 3*d*(t-ti)*(t-ti);
}
