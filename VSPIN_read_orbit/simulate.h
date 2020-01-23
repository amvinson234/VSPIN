#ifndef SIMULATE_H_INCLUDED
#define SIMULATE_H_INCLUDED

#include "planet.h"
#include "constants.h"
#include <iostream>
#include <iomanip>

class Simulate
{
public:

    static void Start(double run_time);

private:

    static void setup();
    static void sim_loop();

    static Planet planet;
    static std::ofstream output;

    static double _run_time;
    static bool exiting;

    static double last_printed_time;
    static double dt_sampling; //time sampling minimum spacing


};

#endif // SIMULATE_H_INCLUDED
