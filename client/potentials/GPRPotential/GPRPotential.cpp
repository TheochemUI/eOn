//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include"GPRPotential.h"
#include"../../subprojects/gprdimer/structures/Structures.h"
#include"../../subprojects/gprdimer/gpr/auxiliary/AdditionalFunctionality.h"

namespace {

    const char *elementArray[] = {"Unknown", "H","He","Li","Be","B","C","N","O",
           "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
           "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
           "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag",
           "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
           "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
           "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
           "Fr","Ra","Ac","Th","Pa","U", NULL};

    // guess the atom type from the atomic mass,
    std::string mass2atom(double atomicmass) {
        return elementArray[int(atomicmass+.5)];
    }


    int symbol2atomicNumber(char const * symbol)
    {
        int i=0;

        while (elementArray[i] != NULL) {
            if (strcmp(symbol, elementArray[i]) == 0) {
                return i;
            }
            i++;
        }
        // invalid symbol
        return -1;
    }

    char const *atomicNumber2symbol(int n)
    {
        return elementArray[n];
    }
}

GPRPotential::GPRPotential(Parameters *p){
    gpr_model = nullptr;
}

void GPRPotential::registerGPRObject(gpr::GaussianProcessRegression *_gpr_model){
    gpr_model = _gpr_model;
}

void GPRPotential::initialize(void){
}

void GPRPotential::cleanMemory(void){
}

std::pair<double, AtomMatrix> GPRPotential::force(AtomMatrix positions, Eigen::VectorXi atomicNrs,
                                    Matrix3d box, int nImages){
    const int nAtoms = positions.rows();
    gpr::Observation obs;
    // TODO: Be better with the number of images
    obs.clear();
    obs.R.resize(positions.rows(), positions.cols());
    obs.G.resize(positions.rows(), positions.cols());
    obs.E.resize(1); // should be nImages
    obs.R.assignFromEigenMatrix(positions);

    // TODO: Benchmark this, see Potential.cpp

    // See GPRTrainTest.cpp for the functions to be called before this
    this->gpr_model->calculatePotential(obs);


    return std::make_pair(obs.E.extractEigenMatrix()(0),
                          obs.G.extractEigenMatrix() * -1);
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void GPRPotential::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box, int nImages){
    throw std::runtime_error("whoops, you called into the wrong force call");
}
