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

    matter->matter2con("movie.con");
    AtomMatrix current;
    AtomMatrix trial;
    double T=temperature;
    T = 100.0;
    numSteps=100;
    stepSize=0.1;
    int accepts=0;
    for(int steps=0; steps<numSteps; steps++){
        current=matter->getPositions();
	trial=current;
        double ecurrent, etrial;
        ecurrent = matter->getPotentialEnergy();
        matter->matter2con("movie.con", true);
        for (int i=0; i<trial.rows(); i++){
            for(int j=0; j<3; j++){
                double d = gaussRandom(0.0, stepSize);
                trial(i,j)+=d;
            }
        }
        matter->setPositions(trial);
        etrial = matter->getPotentialEnergy();
        double de= ecurrent-etrial;
	cout<<"de="<<de<<"\n";
        if(de <= 0.0){
            printf("%i: accept de <= 0.0\n", steps);
	    accepts++;
	    continue;
        }
        double r=randomDouble();
        double kB = parameters->kB;
        double arg = -de/(kB*T);
        printf("arg: %f\n", arg);
        if (arg < -50.0) {
	    matter->setPositions(current);
            printf("%i: reject small arg\n", steps);
	    continue;	  
        }
            
        double p=exp(arg);
        if(r<p){
            printf("%i: accept r<p\n", steps);
	    accepts++;
	    continue;
        }else{
	    matter->setPositions(current);	    
            printf("%i: reject\n", steps);
        }

    }
    cout<<accepts<<"\n";


}
