#include "Matter.h"
#include "Constants.h"
#include "Dynamics.h"
#include "ParallelReplica.h"

ParallelReplica::ParallelReplica(Parameters *params)
{
    parameters = params;
}

ParallelReplica::~ParallelReplica(){ }

void ParallelReplica::run(int bundleNumber)
{
    char buff[STRING_SIZE];
    string reactant_passed("reactant_passed");
    string reactant_output("reactant");

    if (bundleNumber < 0) {
        reactant_passed += ".con";
        reactant_output += ".con";
    }else{
        snprintf(buff, STRING_SIZE, "_%i.con", bundleNumber);
        reactant_passed += buff;
        reactant_output += buff;
    }

    Matter *reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);
     
	printf("Hello,I am ParralelReplica\n");
	printf("mdSteps=%ld\n",parameters->mdSteps);
	printf("mdTemperture=%lf\n",parameters->mdTemperture);
	printf("mdTimeStep=%lf\n",parameters->mdTimeStep);

	Dynamics mdstep(reactant,parameters);
	mdstep.fullSteps();
	//mdstep.onestep();


    printf("Saving result to %s\n", reactant_output.c_str());
    reactant->matter2con(reactant_output);
	
}
