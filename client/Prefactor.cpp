//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <math.h>
#include "Prefactor.h"
#include "Hessian.h"
#include "Log.h"

int Prefactor::getPrefactors(Parameters* parameters, Matter *min1, Matter *saddle, 
                             Matter *min2, double &pref1, double &pref2)
{
    VectorXd min1Freqs, saddleFreqs, min2Freqs;

    // determine which atoms moved in the process
    VectorXi atoms;
    
    if (parameters->prefactorFilterScheme == Prefactor::FILTER_FRACTION) {
        atoms = movedAtomsPct(parameters, min1, saddle, min2);
    }else{
        atoms = movedAtoms(parameters, min1, saddle, min2);
    }

    
    int size = 3*atoms.rows();
    assert(size > 0);

    // calculate min1 frequencies
    Hessian hessian(parameters, min1);
    min1Freqs = hessian.getFreqs(min1, atoms);
    if(min1Freqs.size() == 0)
    {
        log("[Prefactor] Bad hessian: min1\n");
         return -1;
    }
    // remove zero modes
    if(parameters->checkRotation){
        min1Freqs = hessian.removeZeroFreqs(min1Freqs);
    }

    // calculate saddle frequencies
    saddleFreqs = hessian.getFreqs(saddle, atoms);
    if(saddleFreqs.size() == 0)
    {
        log("[Prefactor] Bad hessian: saddle\n");
        return -1;
    }
    // remove zero modes
    if(parameters->checkRotation){
        saddleFreqs = hessian.removeZeroFreqs(saddleFreqs);
    }

    // calculate min2 frequencies
    min2Freqs = hessian.getFreqs(min2, atoms);
    if(min2Freqs.size() == 0)
    {
        if(!parameters->quiet) {
            log("[Prefactor] Bad hessian: min2\n");
        }
        return -1;
    }
    // remove zero modes
    if(parameters->checkRotation){
        min2Freqs = hessian.removeZeroFreqs(min2Freqs);
    }

    // check Hessian sizes
    if((min1Freqs.size() != saddleFreqs.size()) || (min1Freqs.size() != saddleFreqs.size())) {
        if(!parameters->quiet) {
            log("[Prefactor] Bad prefactor: Hessian sizes do not match\n");
        }
        return -1;
    }

    logFreqs(min1Freqs, (char*) "minimum 1");
    logFreqs(saddleFreqs, (char*) "saddle");
    logFreqs(min2Freqs, (char*) "minimum 2");

    // check for correct number of negative modes
    int i, numNegFreq = 0;
    for(i=0; i<size; i++)
    {
        if(min1Freqs(i) < 0) { 
            log("[Prefactor] min1 had negative mode of %e\n", min1Freqs(i));
            numNegFreq++;
        }
    }
    if(numNegFreq != 0)
    {
        log("[Prefactor] Error: %i negative modes at min1\n", numNegFreq);
        return -1;
    }

    numNegFreq = 0;
    for(i=0; i<size; i++)
    {
        if(saddleFreqs(i) < 0) { numNegFreq++; }
    }
    if(numNegFreq != 1)
    {
        log("Error: %i negative modes at saddle\n", numNegFreq);
        return -1;
    }

    numNegFreq = 0;
    for(i=0; i<size; i++)
    {
        if(min2Freqs(i) < 0) { numNegFreq++; }
    }
    if(numNegFreq != 0)
    {
        log("Error: %i negative modes at min2\n", numNegFreq);
        return -1;
    }

    // calculate the prefactors
    pref1 = 1.0;
    pref2 = 1.0;
    
    if (parameters->prefactorRate == Prefactor::RATE_HTST){

        // products are calculated this way in order to avoid overflow
        for(int i=0; i<saddleFreqs.size(); i++)
        {
            pref1 *= min1Freqs[i];
            pref2 *= min2Freqs[i];
            if(saddleFreqs[i]>0)
            {
                pref1 /= saddleFreqs[i];
                pref2 /= saddleFreqs[i];
            }
        }
        pref1 = sqrt(pref1)/(2*M_PI*10.18e-15);
        pref2 = sqrt(pref2)/(2*M_PI*10.18e-15);
    }
    else if (parameters->prefactorRate == Prefactor::RATE_QQHTST){
        float kB_T = parameters->temperature * 8.617332e-5; // eV
        float h_bar = 6.582119e-16; // eV*s
        float h = 4.135667e-15; // eV*s
        float temp = (h_bar / ( 2. * kB_T));
        
        for(int i=0; i<min1Freqs.size(); i++)
        {
            pref1 = pref1 * (sinh (temp * (sqrt(min1Freqs[i]) / 10.18e-15)));
            pref2 = pref2 * (sinh (temp * (sqrt(min2Freqs[i]) / 10.18e-15)));
            
            if(saddleFreqs[i]>0)
            {
                pref1 = pref1 / (sinh ( temp * (sqrt(saddleFreqs[i]) / 10.18e-15)));
                pref2 = pref2 / (sinh ( temp * (sqrt(saddleFreqs[i]) / 10.18e-15)));
            }
        }
        pref1 = 2. * kB_T / (h) * pref1;
        pref2 = 2. * kB_T / (h) * pref2;
    }
    log("[Prefactor] reactant to product prefactor: %.3e\n", pref1);
    log("[Prefactor] product to reactant prefactor: %.3e\n", pref2);
    return 0;
}

void Prefactor::logFreqs(VectorXd freqs, char *name)
{
    log_file("[Prefactor] Frequencies at %s\n", name);
    int i;
    for (i=0;i<freqs.size();i++) {
        log_file("%10.6f ", freqs(i));
        if ((i+1)%5 == 0) {
            log_file("\n");
        }
    }
    log_file("\n");
}

VectorXi Prefactor::movedAtoms(Parameters* parameters, Matter *min1, Matter *saddle, Matter *min2)
{
    long nAtoms = saddle->numberOfAtoms();

    VectorXi moved(nAtoms);
    moved.setConstant(-1);

    AtomMatrix diffMin1 = saddle->pbc(saddle->getPositions() - min1->getPositions());
    AtomMatrix diffMin2 = saddle->pbc(saddle->getPositions() - min2->getPositions());

    diffMin1.cwise() *= saddle->getFree();
    diffMin2.cwise() *= saddle->getFree();

    int nMoved = 0;
    for(int i=0; i<nAtoms; i++)
    {
        if( (diffMin1.row(i).norm() > parameters->prefactorMinDisplacement) || 
            (diffMin2.row(i).norm() > parameters->prefactorMinDisplacement) )
        {
            if(!(moved.cwise() == i).any())
            {
                moved[nMoved] = i;
                nMoved++;
            }
            for(int j=0; j<nAtoms; j++)
            {
                double diffRSaddle = saddle->distance(i,j);

                if(diffRSaddle<parameters->prefactorWithinRadius
                   && (!saddle->getFixed(j)))
                {
                    if(!(moved.cwise() == j).any())
                    {
                        moved[nMoved] = j;
                        nMoved++;
                    }
                }
            }
        }
    }
    return (VectorXi) moved.block(0,0,nMoved,1);
}

VectorXi Prefactor::movedAtomsPct(Parameters* parameters, Matter *min1, Matter *saddle, Matter *min2)
{
    long nAtoms = saddle->numberOfAtoms();
    long nFree = saddle->numberOfFreeAtoms();

    VectorXi moved(nAtoms);
    moved.setConstant(-1);

    AtomMatrix diffMin1 = saddle->pbc(saddle->getPositions() - min1->getPositions());
    AtomMatrix diffMin2 = saddle->pbc(saddle->getPositions() - min2->getPositions());

    diffMin1.cwise() *= saddle->getFree();
    diffMin2.cwise() *= saddle->getFree();

    VectorXd diff(nAtoms);
    diff.setConstant(0.0);
    
    log("[Prefactor] including all atoms that make up %.3f%% of the motion\n", 
            parameters->prefactorFilterFraction*100);
    double sum = 0.0;
    int mini = 0;
    for (int i = 0; i < nAtoms; i++) {
        diff[i] = max(diffMin1.row(i).norm(), diffMin2.row(i).norm());
        sum += diff[i];
        if (diff[i] <= diff[mini]) {
            mini = i;
        }

    }

    log("[Prefactor] sum of atom distances moved %.4f\n", sum);
    log("[Prefactor] max moved atom distance: %.4f\n", diff.maxCoeff());

    int nMoved = 0;
    double d = 0.0;
    while (d/sum <= parameters->prefactorFilterFraction && nMoved < nFree) {
        int maxi = mini;
        for (int i = 0; i < nAtoms; i++) {
            if (diff[i] >= diff[maxi]) {
                if (!(moved.cwise() == i).any()) {
                    maxi = i;
                }
            }
        }
        moved[nMoved] = maxi;
        nMoved++;
        d += diff[maxi];
    }

    int totalAtoms = nMoved;
    for (int i=0;i<nMoved;i++) {
        for(int j=0;j<nAtoms;j++)
        {
            if (moved[i] == j) continue;

            double diffRSaddle = saddle->distance(moved[i],j);
            
            if(diffRSaddle<parameters->prefactorWithinRadius
               && (!saddle->getFixed(j))) {

                if(!(moved.cwise() == j).any()) {
                    moved[totalAtoms++] = j;
                }

            }
        }
    }
    log("[Prefactor] including %i atoms in the hessian (%i moved + %i neighbors)\n", totalAtoms, nMoved, totalAtoms-nMoved);    
    return (VectorXi) moved.block(0,0,totalAtoms,1);
}


VectorXi Prefactor::allFreeAtoms(Matter *matter)
{
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
    return moved.head(nMoved);
}
