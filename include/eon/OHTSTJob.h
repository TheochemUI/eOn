/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once
#include "eon/Eigen.h"
#include "eon/Job.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"

namespace eonc {

/**
 * @file
 * @ingroup Jobs
 *
 * \brief Optimized hyperplanar transition state theory (OH-TST).
 *
 * Implements the reversible-work variational TST of Johannesson and
 * Jonsson, J. Chem. Phys. 115, 9644 (2001), doi:10.1063/1.1415499.
 * A hyperplanar dividing surface in the 3N-dimensional configuration
 * space, parameterized by a point on a reactant->product guideline and
 * a unit normal, is advanced by damped (velocity-zeroing) Verlet
 * dynamics on the plane position (Eqs 5-7) and orientation (Eqs 8-10,
 * Appendix A), driven by canonical-ensemble averages of the force
 * sampled with thermostatted dynamics constrained to the plane. The
 * free-energy barrier follows from the translational and rotational
 * reversible work (Eqs 18-19), the flux prefactor from the
 * direction-dependent effective mass (Eqs 23-24), and the
 * configuration-integral ratio from crossing statistics of an
 * unconstrained reactant trajectory (Eq 22).
 */
class OHTSTJob : public Job {

public:
  OHTSTJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~OHTSTJob(void) = default;
  std::vector<std::string> run(void);

private:
  //! Canonical averages over one thermostatted, plane-constrained
  //! sampling block.
  struct PlaneAverages {
    double fn{0.0};      //!< <F.n>
    VectorXd rotNorm;    //!< <(n.F) R / (alpha |R|^2)>, drives rotation
    VectorXd rotRaw;     //!< <(n.F) R>, integrand of Eq 19
    VectorXd pos;        //!< <r>, anchors Eq 18 and the Eq 12 restart
  };

  PlaneAverages samplePlane(Matter &matter, const VectorXd &gamma,
                            const VectorXd &normal);

  //! Eq 22: Q^ZR/Q^R from crossing statistics of an unconstrained
  //! reactant-basin trajectory through the plane (gammaR, normal).
  double reactantQRatio(Matter &matter, const VectorXd &gammaR,
                        const VectorXd &normal);

  //! Maxwell-Boltzmann draw on the free DOF, then projection onto the
  //! plane (n.v = 0) when a normal is supplied.
  void drawThermalVelocities(VectorXd &vel, const VectorXd *normal);

  //! Eqs 14-17 symmetry restriction: if the configuration is closer
  //! to another equivalent product half-line than to the primary one,
  //! revert the position and reflect the velocity about the mirror
  //! that maps the primary direction onto the offending one. Returns
  //! true when a reflection was applied.
  bool symmetryReflect(const VectorXd &xR, VectorXd &x, VectorXd &v,
                       const VectorXd &xOld, const VectorXd *normal);

  std::vector<VectorXd> m_symDirs; //!< p-hat_i, index 0 = primary
  VectorXd m_symXR;                //!< reactant anchor R of the half-lines

  VectorXd m_masses3N;        //!< per-DOF masses of the free atoms (amu)
  double m_dt{0.0};           //!< integration step, internal units
  double m_kbt{0.0};          //!< k_B T (eV)
  double m_andersenProb{0.0}; //!< per-step Andersen collision probability
  long m_seedState{12345};    //!< LCG state for the thermostat draws
  double uniformDraw();       //!< uniform (0,1)
  double gaussDraw();         //!< standard normal (Box-Muller)
  bool m_gaussHave{false};
  double m_gaussSpare{0.0};
};

} // namespace eonc

using eonc::OHTSTJob;
