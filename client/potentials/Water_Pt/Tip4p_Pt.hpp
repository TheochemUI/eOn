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
    Tip4p_Pt(Parameters* params) : Potential(params), forcefields::ZhuPhilpott<>(8.5, 1.0){};
    // Functions
    // constructor and destructor
    
    // Interface
    void initialize(void) {}
    void cleanMemory(void) {}
    std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                         const VectorXi atmnrs,
                                         const Matrix3d m_box) override;
};


#endif

