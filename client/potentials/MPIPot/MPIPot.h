//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef __MPI_POTENTIAL__
#define __MPI_POTENTIAL__

#include "../../Potential.h"
#include "../../Parameters.h"

class MPIPot : public Potential
{

    public:
        MPIPot(Parameters *p);
		~MPIPot();
        void initialize(){};
        void cleanMemory(void);    
        std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                             const VectorXi atmnrs,
                                             const Matrix3d m_box) override;

    private:
        int potentialRank;
        double poll_period;
};

#endif
