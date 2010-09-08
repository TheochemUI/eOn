/** @file
Wrapper for Eon
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/

#ifndef FORCEFIELDS_ZHU_PHILPOTT_FOR_EON_HPP
#define FORCEFIELDS_ZHU_PHILPOTT_FOR_EON_HPP
#include "zhu_philpott.hpp"
#include "tip4p.hpp"
#include "../../PotentialsInterface.h"

class ZpIce : public PotentialsInterface, private forcefields::ZhuPhilpott<> {
public:
    ZpIce();
    // Functions
    // constructor and destructor
    
    // To satify interface
    void initialize(void) {}
    void cleanMemory(void) {}
    void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);
};

class Tip4p : public PotentialsInterface, private forcefields::Tip4p {
public:
    Tip4p();
    // Functions
    // constructor and destructor
    
    // To satify interface
    void initialize(void) {}
    void cleanMemory(void) {}
    void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);
};
#endif

