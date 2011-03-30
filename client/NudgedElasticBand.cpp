//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "NudgedElasticBand.h"
#include <cassert>

using namespace helper_functions;

NEB::NEB(Matter const *matterInitial, Matter const *matterFinal, Parameters *params)
{
    parameters = params;
    images = parameters -> nebImages;
    Matter *neb[images+2];
    for(long i=0; i<images+2; i++){
        neb[i] = new Matter(parameters);
    }

    matterInitial = new Matter(parameters);
    matterFinal = new Matter(parameters);
    *matterInitial = *matterA;
    *matterFinal = *matterB;
    nAtoms = matterInitial->numberOfAtoms();
    assert(nAtoms == matterFinal->numberOfAtoms());


    direction.resize(nAtoms, 3);
    rotationalPlane.resize(nAtoms, 3);
    direction.setZero();
    rotationalPlane.setZero();
    totalForceCalls = 0;
}

NEB::~NEB()
{
    delete matterInitial;
    delete matterFinal;
}

void NEB::compute(void)
{
    long rotations = 0;
    long forceCallsCenter;
    long forceCallsDimer;
    double rotationalForce1;
    double rotationalForce2;
    double curvature, rotationalForceChange, forceDimer, rotationAngle;
    double lengthRotationalForceOld;
    double torque = 0;
    bool doneRotating = false;

    *matterCenter = *matter;
    rotationalForceChange = forceDimer = rotationAngle = curvature = 0;
    rotationalForce1 = 0;
    rotationalForce2 = 0;
    Matrix<double, Eigen::Dynamic, 3> rotationalForce(nAtoms,3);
    Matrix<double, Eigen::Dynamic, 3> rotationalForceOld(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> rotationalPlaneOld(nAtoms, 3);
    rotationalForce.setZero();
    rotationalForceOld.setZero();
    rotationalPlaneOld.setZero();
    initialDirection.normalize();
    direction = initialDirection;

    statsAngle = 0;
    lengthRotationalForceOld = 0;
    forceCallsCenter = matterCenter->getForceCalls();
    forceCallsDimer = matterDimer->getForceCalls();

    return;
}

