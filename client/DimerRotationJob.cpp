//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
#include "DimerRotationJob.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Constants.h"
#include "HelperFunctions.h"

using namespace helper_functions;

DimerRotationJob::DimerRotationJob(Parameters *params)
{
    parameters = params;
}

DimerRotationJob::~DimerRotationJob(){ }

void DimerRotationJob::run(int bundleNumber)
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

    // Create the rotationPlaneNorm.
    Matrix<double, Eigen::Dynamic, 3> rotationPlaneNorm;        
    rotationPlaneNorm = mode;
    for(int i = 0; i < nall / 3; i ++)
    {
        for(int j = 0; j < 3; j++)
        {
            rotationPlaneNorm(i, j) = randomDouble(10.0);
        }
    }
    rotationPlaneNorm.normalize(); 
    rotationPlaneNorm = makeOrthogonal(rotationPlaneNorm, mode);
    rotationPlaneNorm.normalize();

    double projectedForceA, projectedForceB, curvature1;
    double dRots[] = {0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 1.0, 1.5, 2.0};
    int ndRots = 16;
    Matrix<double, Eigen::Dynamic, 3> dimer2;    
    Matrix<double, Eigen::Dynamic, 3> forceA;    
    Matrix<double, Eigen::Dynamic, 3> forceB;    

    forceA = displacement->getForces();
    dimer2 = dimer1 + mode * parameters->dimerSeparation;
    displacement->setPositions(dimer2);
    forceB = displacement->getForces();
    projectedForceA = (mode.cwise() * forceA).sum();
    projectedForceB = (mode.cwise() * forceB).sum();
    curvature1 = (projectedForceB - projectedForceA) / (2.0 * parameters->dimerSeparation);    
    
    // Loop over values of dimer dRot and print the output to results.dat.
    FILE *results = fopen("results.dat", "w");
    fprintf(results, "%14s    %15s\n", "dRot", "dCurvature/dRot");
    Matrix<double, Eigen::Dynamic, 3> mode2;    
    double cosAngle, sinAngle, curvature2;
    for (int dRoti = 0; dRoti < ndRots; dRoti++)
    {
        cosAngle = cos(dRots[dRoti]);
        sinAngle = sin(dRots[dRoti]);
        mode2 = mode;
        mode2 = mode2 * cosAngle + rotationPlaneNorm * sinAngle;
        mode2.normalize();
        dimer2 = dimer1 + mode2 * parameters->dimerSeparation;
        displacement->setPositions(dimer2);
        forceB = displacement->getForces();
        projectedForceA = (mode2.cwise() * forceA).sum();
        projectedForceB = (mode2.cwise() * forceB).sum();
        curvature2 = (projectedForceB - projectedForceA) / (2.0 * parameters->dimerSeparation);        
        fprintf(results, "%14.8f    %15.8f\n", dRots[dRoti], (curvature1 - curvature2) / dRots[dRoti]);
    }
    fclose(results);
    
}
