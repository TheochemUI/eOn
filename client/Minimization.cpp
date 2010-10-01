#include "Minimization.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Constants.h"

Minimization::Minimization(Parameters *params)
{
    parameters = params;
}

Minimization::~Minimization(){ }

void Minimization::run(int bundleNumber)
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

    printf("\nBeginning minimization of %s\n", reactant_passed.c_str());

    ConjugateGradients cgMin(reactant, parameters);
    cgMin.fullRelax();

    if (cgMin.isItConverged(parameters->getConverged_Relax())) {
        printf("Minimization converged within tolerence\n");
    }else{
        printf("Minimization did not converge to tolerence!\n"
               "Maybe try to increase maximum_iterations?\n");
    }

    printf("Saving result to %s\n", reactant_output.c_str());
    reactant->matter2con(reactant_output);
}
