//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "PrefactorJob.h"
#include "Prefactor.h"
#include "Matter.h"
#include "Hessian.h"
#include "Potential.h"

const char PrefactorJob::PREFACTOR_REACTANT[] = "reactant";
const char PrefactorJob::PREFACTOR_SADDLE[]   = "saddle";
const char PrefactorJob::PREFACTOR_PRODUCT[]  = "product";

PrefactorJob::PrefactorJob(Parameters *params)
{
    parameters = params;
}

PrefactorJob::~PrefactorJob()
{
}

std::vector<std::string> PrefactorJob::run(void)
{
    std::vector<std::string> returnFiles;
    VectorXd freqs;

    string reactant_passed("reactant_passed.con");
    string saddle_passed("saddle_passed.con");
    string product_passed("product_passed.con");

    Matter *reactant = new Matter(parameters);
    Matter *saddle = new Matter(parameters);
    Matter *product = new Matter(parameters);

    
    VectorXi atoms;
    if (parameters->prefactorAllFreeAtoms)
    {
        // it is sufficient to pass the configuration 
        // for which the frequencies should be determined
        string matter_passed;
        if (parameters->prefactorConfiguration == 
            PrefactorJob::PREFACTOR_REACTANT)
        {
            matter_passed = reactant_passed;
        }
        else if (parameters->prefactorConfiguration == 
                 PrefactorJob::PREFACTOR_SADDLE)
        {
            matter_passed = saddle_passed;
        }
        else if (parameters->prefactorConfiguration == 
                 PrefactorJob::PREFACTOR_PRODUCT)
        {
            matter_passed = product_passed;
        }
        reactant->con2matter(matter_passed);
        saddle->con2matter(matter_passed);
        product->con2matter(matter_passed);

        // account for all free atoms
        atoms = allFreeAtoms(reactant);        
    }
    else
    {
        reactant->con2matter(reactant_passed);
        saddle->con2matter(saddle_passed);
        product->con2matter(product_passed);

        // determine which atoms moved in the process
        atoms = movedAtoms(parameters, reactant, saddle, product);
    }
    assert(3*atoms.rows() > 0);
    
    // calculate frequencies
    if (parameters->prefactorConfiguration == 
        PrefactorJob::PREFACTOR_REACTANT)
    {
        Hessian hessian(parameters, reactant);    
        freqs = hessian.getFreqs(reactant, atoms);
    }
    else if (parameters->prefactorConfiguration == 
             PrefactorJob::PREFACTOR_SADDLE)
    {
        Hessian hessian(parameters, saddle);    
        freqs = hessian.getFreqs(saddle, atoms);
    }    
    else if (parameters->prefactorConfiguration == 
             PrefactorJob::PREFACTOR_PRODUCT)
    {
        Hessian hessian(parameters, product);    
        freqs = hessian.getFreqs(product, atoms);
    }

    bool failed = freqs.size() != 3*atoms.rows();

    FILE *fileResults;
    FILE *fileFreq;

    std::string results_file("results.dat");
    std::string freq_file("freq.dat");

    returnFiles.push_back(results_file);
    returnFiles.push_back(freq_file);

    fileResults = fopen(results_file.c_str(), "wb");
    fileFreq = fopen(freq_file.c_str(), "wb");

    fprintf(fileResults, "%s good\n", failed ? "false" : "true");
    fprintf(fileResults, "%d force_calls\n", Potential::fcalls);

    if(!failed)
    {
        for(int i=0; i<freqs.size(); i++)
        {
            if ( 0. < freqs[i] )
            {
                fprintf(fileFreq, "%f\n", sqrt(freqs[i])/(2*M_PI*10.18e-15));
            }
            else
            {
                fprintf(fileFreq, "%f\n", -sqrt(-freqs[i])/(2*M_PI*10.18e-15));                
            }
        }
    }

    delete reactant;
    delete product;
    delete saddle;

    return returnFiles;
}

