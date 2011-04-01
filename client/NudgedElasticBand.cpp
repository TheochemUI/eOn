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

NEB::NEB(Matter const *matterA, Matter const *matterB, Parameters *params)
{
    parameters = params;
    images = parameters -> nebImages;
    Matter *neb[images+2];
    AtomMatrix tangent[images+2];
    for(long i=0; i<images+2; i++){
        neb[i] = new Matter(parameters);
    }
    matterInitial = new Matter(parameters);
    matterFinal = new Matter(parameters);
    *matterInitial = *matterA;
    *matterFinal = *matterB;
    nAtoms = matterInitial->numberOfAtoms();
    assert(nAtoms == matterFinal->numberOfAtoms());

    imageSep = matterFinal->getPosition()
    for(long i=0; i<images+2; i++){
        tangent[i].resize(nAtoms,3);

}

NEB::~NEB()
{
    delete matterInitial;
    delete matterFinal;
    for(long i=0; i<images+2; i++){
        delete neb[i];
    }
}

void NEB::compute(void)
{
    return;
}

