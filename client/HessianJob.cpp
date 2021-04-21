#include "HessianJob.h"
#include "Matter.h"
#include "Hessian.h"
#include "Potential.h"

HessianJob::HessianJob(Parameters *params)
{
    parameters = params;
}

HessianJob::~HessianJob()
{
}

std::vector<std::string> HessianJob::run(void)
{
    string matter_in("pos.con");

    std::vector<std::string> returnFiles;

    Matter *matter = new Matter(parameters);

    matter->con2matter(matter_in);

    Hessian hessian(parameters, matter);
    long nAtoms = matter->numberOfAtoms();
    
    VectorXi moved(nAtoms);
    moved.setConstant(-1);
    
    int nMoved = 0;
    for(int i=0; i<nAtoms; i++)
    {
        if(!matter->getFixed(i))
        {
            moved[nMoved] = i;
            nMoved++;
        }
    }
    moved = moved.head(nMoved);
    hessian.getFreqs(matter, moved);

    FILE *fileResults;
//    FILE *fileMode;

    std::string results_file("results.dat");

    returnFiles.push_back(results_file);

    fileResults = fopen(results_file.c_str(), "wb");

    fprintf(fileResults, "%d force_calls\n", Potential::fcalls);
    fclose(fileResults);

    delete matter;

    return returnFiles;
}

