#include "HessianJob.h"
#include "Matter.h"
#include "Hessian.h"
#include "Potentials.h"

HessianJob::HessianJob(Parameters *params)
{
    parameters = params;
}

HessianJob::~HessianJob()
{
}

void HessianJob::run(int bundleNumber)
{
    char buff[STRING_SIZE];
    string reactant_passed("reactant_passed");
    string saddle_passed("saddle_passed");
    string product_passed("product_passed");

    if (bundleNumber < 0) {
        reactant_passed += ".con";
        saddle_passed += ".con";
        product_passed += ".con";
    }else{
        snprintf(buff, STRING_SIZE, "_%i.con", bundleNumber);
        reactant_passed += buff;
        saddle_passed += buff;
        product_passed += buff;
    }

    Matter *reactant = new Matter(parameters);
    Matter *saddle = new Matter(parameters);
    Matter *product = new Matter(parameters);

    reactant->con2matter(reactant_passed);
    saddle->con2matter(saddle_passed);
    product->con2matter(product_passed);

    Hessian hessian(reactant, saddle, product, parameters);
    double modeProduct = hessian.getModeProduct(parameters->hessianKind);

    bool failed = modeProduct<0;

    FILE *fileResults;

    char filename[STRING_SIZE];

    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "results_%i.dat", bundleNumber);
    }else{
        strncpy(filename, "results.dat", STRING_SIZE);
    }
	fileResults = fopen(filename, "wb");

	fprintf(fileResults, "%s good\n", failed ? "false" : "true");
	fprintf(fileResults, "%d force_calls\n", Potentials::fcalls);
	if(!failed)
    {
	    fprintf(fileResults, "%f mode_product\n", modeProduct);
    }

    delete reactant;
    delete product;
    delete saddle;
}
    
    

    

