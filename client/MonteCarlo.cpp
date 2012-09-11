//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "MonteCarlo.h"
#include "Log.h"

using namespace helper_functions;

MonteCarlo::MonteCarlo(Matter const *matterIn, Parameters *params)
{
    parameters    = params;
    matter = new Matter(parameters);
    *matter = *matterIn;
}

MonteCarlo::~MonteCarlo()
{
}

void MonteCarlo::run(int numSteps, double temperature, double stepSize)
{

}
