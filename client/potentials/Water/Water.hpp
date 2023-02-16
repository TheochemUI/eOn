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
    Tip4p(Parameters* params) : Potential(params), forcefields::Tip4p(8.5, 1.0){};
    // Functions
    // constructor and destructor
    
    // To satisfy interface
    void cleanMemory(void) {}
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
    std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                         const VectorXi atmnrs,
                                         const Matrix3d box) override {
      double energy{std::numeric_limits<double>::infinity()};
      long nAtoms{pos.rows()};
      AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
      this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy,
                  box.data());
      return std::make_pair(energy, forces);
    };
};

class SpceCcl : public Potential, private forcefields::SpceCcl {
public:
    SpceCcl(Parameters* params) : Potential(params), forcefields::SpceCcl(8.5, 1.0){}
    // Functions
    // constructor and destructor
    
    // To satisfy interface
    void cleanMemory(void) {}
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
    std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                         const VectorXi atmnrs,
                                         const Matrix3d box) override {
      double energy{std::numeric_limits<double>::infinity()};
      long nAtoms{pos.rows()};
      AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
      this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy,
                  box.data());
      return std::make_pair(energy, forces);
    };
};

#endif

