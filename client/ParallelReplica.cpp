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
    string reactant_output("product");

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
     
	printf("Now running Parralel Replica Dynamics\n");

	Dynamics mdstep(reactant,parameters);
	mdstep.fullSteps();

    printf("Saving result to %s\n", reactant_output.c_str());
    reactant->matter2con(reactant_output);
	
}
