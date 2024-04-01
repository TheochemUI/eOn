// An interface to the GPDimer library

#include "AtomicGPDimer.h"

const char AtomicGPDimer::OPT_SCG[] = "scg";
const char AtomicGPDimer::OPT_LBFGS[] = "lbfgs";

void AtomicGPDimer::compute(std::shared_ptr<Matter> a_matter,
                            AtomMatrix a_initialDirectionAtomMatrix) {
  m_atoms_config = helpers::gproptim::input::eon_matter_to_atmconf(a_matter);
  *m_dimer_center = *a_matter;
  m_R_init.resize(1, m_dimer_center->getPositionsFree().rows() *
                         m_dimer_center->getPositionsFree().cols());
  int counter = 0;
  for (int i = 0; i < m_dimer_center->getPositionsFree().rows(); ++i) {
    for (int j = 0; j < m_dimer_center->getPositionsFree().cols(); ++j) {
      m_R_init[counter++] = m_dimer_center->getPositionsFree()(i, j);
    }
  }
  m_init_middle_point.clear();
  m_init_middle_point.R = m_R_init;
  // TODO: This is now conditional! Only if there's no existing data.
  // HACK: Obviously this can be improved, add it as a parameter, check for existence etc.
  // m_init_observations.clear();
  // m_problem_setup.cutOffEnergy(E_level, middle_point.E);

  m_problem_setup.activateFrozenAtoms(m_R_init, params->gprActiveRadius,
                                      m_atoms_config);
  m_orient_init.clear();
  m_orient_init.resize(m_dimer_center->getPositionsFree().rows(),
                       m_dimer_center->getPositionsFree().cols());
  AtomMatrix freeOrient(m_dimer_center->numberOfFreeAtoms(), 3);
  int i, j = 0;
  for (i = 0; i < m_dimer_center->numberOfAtoms(); i++) {
    if (!m_dimer_center->getFixed(i)) {
      freeOrient.row(j) = a_initialDirectionAtomMatrix.row(i);
      j++;
      if (j == m_dimer_center->numberOfFixedAtoms()) {
        break;
      }
    }
  }
  m_orient_init.resize(1, freeOrient.rows() * freeOrient.cols());
  counter = 0;
  for (int i = 0; i < freeOrient.rows(); ++i) {
    for (int j = 0; j < freeOrient.cols(); ++j) {
      m_orient_init[counter++] = freeOrient(i, j);
    }
  }
  m_atomic_dimer.initialize(m_gp_params, m_init_observations,
                            m_init_middle_point, m_orient_init, m_atoms_config);

  auto potential = helper_functions::makePotential(params->potential, params);
  m_atomic_dimer.execute(*potential);
  // Forcefully set the right positions
  a_matter->setPositionsFreeV(m_atomic_dimer.getFinalCoordOfMidPoint());
  this->totalIterations = m_atomic_dimer.getIterations();
  this->totalForceCalls = m_atomic_dimer.getTotalForceCalls();
  return;
}

double AtomicGPDimer::getEigenvalue() {
  return m_atomic_dimer.getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector() {
  return m_atomic_dimer.getFinalOrientation()->extractEigenMatrix();
}
