/** @file
Wrapper for Eon
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/

#ifndef WATER_FOR_EON_HPP
#define WATER_FOR_EON_HPP
#include "tip4p_ccl.hpp"
#include "spce_ccl.hpp"
#include "../../Potential.h"


class Tip4p : public Potential, private forcefields::Tip4p {
public:
    Tip4p();
    // Functions
    // constructor and destructor
    
    // To satify interface
    void initialize(void) {}
    void cleanMemory(void) {}
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};

class SpceCcl : public Potential, private forcefields::SpceCcl {
public:
    SpceCcl();
    // Functions
    // constructor and destructor
    
    // To satisfy interface
    void initialize(void) {}
    void cleanMemory(void) {}
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};

#endif

