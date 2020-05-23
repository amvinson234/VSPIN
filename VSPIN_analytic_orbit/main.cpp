#include <iostream>
#include <iomanip>
#include "damp.h"
#include "constants.h"
#include "simulate.h"

int main()
{

    //std::stringstream runID;
    //runID << std::getenv("SGE_TASK_ID");
    //std::string run = runID.str();
    std::string run = "";

    Simulate::Start(run, 3000000.);

    return 0;
}
