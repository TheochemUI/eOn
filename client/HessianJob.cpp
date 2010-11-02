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
    FILE *fileMode;

    char resultsname[STRING_SIZE];
    char modename[STRING_SIZE];

    if (bundleNumber != -1) {
        snprintf(resultsname, STRING_SIZE, "results_%i.dat", bundleNumber);
        snprintf(modename, STRING_SIZE, "mode_%i.dat", bundleNumber);
    }else{
        strncpy(resultsname, "results.dat", STRING_SIZE);
        strncpy(modename, "mode.dat", STRING_SIZE);
    }
	fileResults = fopen(resultsname, "wb");
	fileMode = fopen(modename, "wb");

	fprintf(fileResults, "%s good\n", failed ? "false" : "true");
	fprintf(fileResults, "%d force_calls\n", Potentials::fcalls);
    fprintf(fileResults, "%d hessian_size\n", hessian.getHessian(parameters->hessianKind).rows());
	if(!failed)
    {
	    fprintf(fileResults, "%f mode_product\n", modeProduct);
        
        VectorXd modes = hessian.getModes(parameters->hessianKind);
        for(int i=0; i<modes.size(); i++)
        {
            fprintf(fileMode, "%f\n", modes[i]);
        }
    }

    delete reactant;
    delete product;
    delete saddle;
}
    
    

    

