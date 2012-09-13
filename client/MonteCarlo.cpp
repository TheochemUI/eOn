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
#include <stdio.h>

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

    remove("movie.con");
    AtomMatrix p1=matter->getPositions();
    AtomMatrix store;
    double T=temperature;
    T = 1.0;
    numSteps=1000;
    stepSize=1.0;
    int accepts=0;
    for(int steps=0; steps<numSteps; steps++){
        store=matter->getPositions();
        double e1, e2;
        e1 = matter->getPotentialEnergy();
        matter->matter2con("movie.con", true);
        for (int i=0; i<p1.rows(); i++){
            for(int j=0; j<3; j++){
                double d = gaussRandom(0.0, stepSize);
                p1(i,j)+=d;
            }
        }

        matter->setPositions(p1);
        e2 = matter->getPotentialEnergy();
        double de= e2-e1;
        if(de <= 0.0){
            printf("%i: accept de <= 0.0\n", steps);
            accepts++;
            continue;
        }
        double r=randomDouble();
        double kb= 8.6173324e-5;
        double arg = -de/(kb*T);
        printf("arg: %f\n", arg);
        if (arg < -50.0) {
            matter->setPositions(store);
            printf("%i: reject\n", steps);
            continue;
        }
            
        double p=exp(arg);
        if(r<p){
            accepts++;
            printf("%i: accept r<p\n", steps);
            matter->setPositions(p1);
        }else{
            matter->setPositions(store);
            printf("%i: reject\n", steps);
        }

    }
    cout<<accepts<<"\n";

    //double e1, e2;
    //e1 = matter->getPotentialEnergy();
    //p1 = matter->getPositions();
    //p1 += d;
    ////d(0,0) = gaussRandom(0.0, stepSize);
    //matter->setPositions(p1);
    //e2 = matter->getPotentialEnergy();
    //

    //matter->matter2con(filename, append);
    //matter->matter2con("movie.con", true);

    //double de = e2 - e1;
    //r = randomDouble();
    //p = exp(-de/(kB*T));
    //if (r<p) { accept; }

}
