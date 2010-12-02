//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------
#include "DimerDrJob.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Constants.h"

DimerDrJob::DimerDrJob(Parameters *params)
{
    parameters = params;
}

DimerDrJob::~DimerDrJob(){ }

void DimerDrJob::run(int bundleNumber)
{
    // No bundling for this job, so bundleNumber is ignored.
    
    // Load the displacement con file and get the position.
    Matter *displacement = new Matter(parameters);
    displacement->con2matter("displacement_passed.con");
    Matrix<double, Eigen::Dynamic, 3> dimer1;    
    dimer1 = displacement->getPositions();        

    // Load and normalize the mode.
    FILE *modeFile = fopen("mode_passed.dat", "r");
    Matrix<double, Eigen::Dynamic, 3> mode;    
    long nall = 0, nfree = 0;
    fscanf(modeFile, "%ld %ld", &nall, &nfree);
    mode.resize(nall/3, 3);
    mode.setZero();
    for (int i = 0; i < nall / 3; i++) 
    {
        for(int j = 0; j < 3; j++)
        {
            fscanf(modeFile, "%lf", &mode(i, j));
        }
    }
    fclose(modeFile);
    mode.normalize();

    double projectedForceA, projectedForceB, curvature;
    double dRs[] = {0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 1.0, 1.5, 2.0};
    int ndRs = 16;
    Matrix<double, Eigen::Dynamic, 3> dimer2;    
    Matrix<double, Eigen::Dynamic, 3> forceA;    
    Matrix<double, Eigen::Dynamic, 3> forceB;    

    forceA = displacement->getForces();
    
    // Loop over values of dimer dR and print the output to results.dat.
    FILE *results = fopen("results.dat", "w");
    fprintf(results, "%14s    %14s\n", "dR", "curvature");
    for (int dRi = 0; dRi < ndRs; dRi++)
    {
        dimer2 = dimer1 + mode * dRs[dRi];
        displacement->setPositions(dimer2);
        forceB = displacement->getForces();
        projectedForceA = (mode.cwise() * forceA).sum();
        projectedForceB = (mode.cwise() * forceB).sum();
        curvature = (projectedForceB - projectedForceA) / (2.0 * dRs[dRi]);        
        fprintf(results, "%14.8f    %14.8f\n", dRs[dRi], curvature);
    }
    fclose(results);
}
