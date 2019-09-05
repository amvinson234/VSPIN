#include "spline.h"
#include <iostream>

void spline(std::vector<double> &time, std::vector<double> &data, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d)
{
    std::vector<double> h, alpha, l, mu, z;
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

int sim_loop(std::string file_name)
{
    /*********************************************************************************************************************
    *RETURNS JUST LINE 0 FOR NOW! MAY NEED TO CHANGE THIS LATER IF I WANT TO LOOP OVER SIMS (BUT MAY NOT BE NECESSARY)   *
    *********************************************************************************************************************/
    return 0;

    std::ifstream input(file_name.c_str());
    std::string last_line, line;
    getline(input,line);
    int ctr = 0;
    double time, n, ecc, peri;
    double n_prev_last, e_prev_last, peri_prev_last, time_prev_last,h;
    input >> time;
    do
    {
        ctr++;
        h = abs(time_prev_last - time);
        time_prev_last = time;
        n_prev_last = n;
        e_prev_last = ecc;
        peri_prev_last = peri;
        input >> n >> ecc >> peri;
    }while(input >> time);
    input.close();

    size_t sz;
    double time_last = time;
    double n_last = n;
    double e_last = ecc;
    double peri_last = peri;
    double d_n_dt_last = (n_last - n_prev_last)/h;
    double d_e_dt_last = (e_last - e_prev_last)/h;
    double d_peri_dt_last = (peri_last - peri_prev_last)/h;

    input.open(file_name.c_str());
    getline(input,line);
    double n_prev, e_prev, peri_prev, time_prev;
    double score = 1000.;
    double line_best = 0;

    for (int i = 1; i < ctr/2; i++)
    {

        n_prev = n;
        e_prev = ecc;
        peri_prev = peri;
        input >> time >> n >> ecc >> peri;
        if(i > 0)
        {
            double d_n_dt = (n - n_prev)/h;
            double d_e_dt = (ecc - e_prev)/h;
            double d_peri_dt = (peri - peri_prev)/h;
            double score_test = std::pow(std::abs((d_n_dt-d_n_dt_last)/d_n_dt_last),2) + std::pow(std::abs((d_e_dt-d_e_dt_last)/d_e_dt_last),2) + std::pow(std::abs((d_peri_dt - d_peri_dt_last)/d_peri_dt_last),2)+
            std::pow(std::abs((n_last-n)/n_last),2) + std::pow(std::abs((ecc-e_last)/e_last),2) + std::pow(std::abs((peri-peri_last)/peri_last),2);
            //if(i % 5000 == 0)
            if (score_test < score)
            {
                score = score_test;
                line_best = i;
            }

        }
    }
    input.close();
    return line_best;


    /*psuedo:
    Find last line in output file.
    Look at n, n-dot, (peri, peri-dot, ecc, ecc-dot??)
    Loop through first quarter or so of the file.
    Find line which matches best with the values of the last line.
        idea: find min(n/n_end + n-dot/n-dot_end + (... other variables))
    return line number.
    */
}

