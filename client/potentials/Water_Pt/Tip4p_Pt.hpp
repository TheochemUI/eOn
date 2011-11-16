/** @file
Wrapper for Eon
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/

#ifndef TIP4P_PT
#define TIP4P_PT
#include "zhu_philpott.hpp"
#include "../../Potential.h"


class Tip4p_Pt : public Potential, private forcefields::ZhuPhilpott<> {
public:
    Tip4p_Pt();
    // Functions
    // constructor and destructor
    
    // To satify interface
    void initialize(void) {}
    void cleanMemory(void) {}
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
};


#endif

