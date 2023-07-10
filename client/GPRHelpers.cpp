#include "GPRHelpers.h"

#include <map>
#include <set>
#include <unordered_map>

namespace helpers::gproptim::input {
gpr::InputParameters eon_parameters_to_gpr(shared_ptr<Parameters> a_params) {
  gpr::InputParameters gppi;
  // Problem parameters
  gppi.actdist_fro.value = a_params->gprActiveRadius;
  gppi.dimer_sep.value = a_params->gprDimerSep;
  gppi.method_rot.value = a_params->gprDimerRotOptMethod;
  gppi.method_trans.value = a_params->gprDimerTransOptMethod;
  gppi.param_trans.value[0] = a_params->gprDimerConvStep;
  gppi.param_trans.value[1] = a_params->gprDimerMaxStep;
  gppi.T_dimer.value = a_params->gprForceThreshold;
  gppi.initrot_nogp.value = a_params->gprDimerInitRotGP;
  gppi.T_anglerot_init.value = a_params->gprDimerConvergedAngle;
  gppi.num_iter_initrot.value = a_params->gprDimerInitRotationsMax;
  gppi.inittrans_nogp.value = a_params->gprDimerInitTransGP;
  gppi.T_anglerot_gp.value = a_params->gprDimerRelaxConvAngle;
  gppi.num_iter_rot_gp.value = a_params->gprDimerRelaxRotationsMax;
  gppi.divisor_T_dimer_gp.value = a_params->gprDimerDivisorTdimerGP;
  gppi.disp_max.value = a_params->gprDimerMidpointMaxDisp;
  gppi.ratio_at_limit.value = a_params->gprDimerRatioAtLimit;
  gppi.num_bigiter.value = a_params->gprDimerMaxOuterIterations;
  gppi.num_iter.value = a_params->gprDimerMaxInnerIterations;
  gppi.islarge_num_iter.value = a_params->gprDimerManyIterations;
  // GPR Parameters
  gppi.gp_sigma2.value = a_params->gprDimerSigma2;
  gppi.jitter_sigma2.value = a_params->gprDimerJitterSigma2;
  gppi.sigma2.value = a_params->gprDimerNoiseSigma2;
  gppi.prior_mu.value = a_params->gprDimerPriorMu;
  gppi.prior_nu.value = a_params->gprDimerPriorNu;
  gppi.prior_s2.value = a_params->gprDimerPriorSigma2;
  gppi.check_derivative.value = a_params->gprOptCheckDerivatives;
  gppi.max_iter.value = a_params->gprOptMaxIterations;
  gppi.tolerance_func.value = a_params->gprOptTolFunc;
  gppi.tolerance_sol.value = a_params->gprOptTolSol;
  gppi.lambda_limit.value = a_params->gprOptLambdaLimit;
  gppi.lambda.value = a_params->gprOptLambdaInit;
  // Prune
  gppi.use_prune.value = a_params->gprUsePrune;
  gppi.start_prune_at.value = a_params->gprPruneBegin;
  gppi.nprune_vals.value = a_params->gprPruneNVals;
  gppi.prune_threshold.value = a_params->gprPruneThreshold;
  // Debugging
  gppi.report_level.value = a_params->gprReportLevel;
  gppi.debug_level.value = a_params->gprDebugLevel;
  gppi.debug_output_dir.value = a_params->gprDebugOutDir;
  gppi.debug_output_file_R.value = a_params->gprDebugPosFile;
  gppi.debug_output_file_E.value = a_params->gprDebugEnergyFile;
  gppi.debug_output_file_G.value = a_params->gprDebugGradFile;
  gppi.debug_output_file_extension.value = a_params->gprDebugOutExt;
  gppi.debug_offset_from_mid_point.value = a_params->gprDebugOffsetMidPoint;
  gppi.debug_dy.value = a_params->gprDebugDy;
  gppi.debug_dz.value = a_params->gprDebugDz;
  return gppi;
}
// FIXME: Take in the active / inactive pairs / atomtypes
gpr::AtomsConfiguration eon_matter_to_atmconf(shared_ptr<Matter> a_matter) {
  //   AtomsConfiguration a;
  //   aux::ProblemSetUp problem_setup;
  //   std::vector<int> atomnrs;
  //   a.positions.resize(a_matter->getPositions().rows(),
  //                      a_matter->getPositions().cols());
  //   a.is_frozen.resize(a_matter->numberOfAtoms());
  //   a.id.resize(a_matter->numberOfAtoms());
  //   a.positions.assignFromEigenMatrix(a_matter->getPositions());
  //   for (auto i = 0; i < a_matter->numberOfAtoms(); i++) {
  //     atomnrs.push_back(a_matter->getAtomicNr(i));
  //     a.is_frozen[i] = a_matter->getFixed(i);
  //     a.id[i] = i + 1;
  //   }
  //   a.atoms_froz_active.clear();
  //   a.atoms_mov.resize(a_matter->numberOfFreeAtoms());
  //   // FIXME: Might have more than one kind of freely moving atom
  //   a.atoms_mov.type.set(0); // Corresponds to H in the CuH example, 0 for Pt
  //   a.atoms_froz_active.clear();
  //   // Atomtypes
  //   Index_t n_at = std::set<int>(atomnrs.begin(), atomnrs.end()).size();
  //   a.pairtype.resize(n_at, n_at);
  //   std::cout << "n_at: " << n_at << "\n";
  //   a.pairtype.set(EMPTY);
  //   a.n_pt = 0;
  //   problem_setup.setPairtypeForMovingAtoms(a.atoms_mov.type, a.n_pt,
  //   a.pairtype); a.atoms_froz_inactive.resize(3 *
  //   a_matter->numberOfFixedAtoms()); for (auto i = 0; i <
  //   a_matter->numberOfFixedAtoms(); ++i) {
  //     a.atoms_froz_inactive.positions.set(0, i,
  //                               {a_matter->getPosition(i, 0),
  //                                a_matter->getPosition(i, 1),
  //                                a_matter->getPosition(i, 2)});
  //   }
  //   // FIXME: Might have more than one kind OR the SAME KIND (Pt)
  //   a.atoms_froz_inactive.type.set(0); // 1 for Cu in the example, 0 for Pt
  //   return a;

  gpr::AtomsConfiguration atoms_config;
  aux::ProblemSetUp problem_setup;
  gpr::Index_t number_of_mov_atoms;
  gpr::Index_t number_of_fro_atoms;
  std::set<int> unique_atomtypes;
  gpr::Index_t n_at;
  vector<int> atomnrs;
  unordered_map<int, int>
      atype_to_gprd_atype; //!> Remember that the atom type in EON is the real
                           //! atomic number, while in GPR Dimer it is a set of
                           //! values from 0 to n-1 so this is EON
  int fake_atype;          //!> False "atomtype" for GPR Dimer

  atoms_config.clear();
  atoms_config.is_frozen.resize(a_matter->numberOfAtoms());
  atoms_config.id.resize(a_matter->numberOfAtoms());
  atoms_config.atomicNrs.resize(a_matter->numberOfAtoms());
  for (auto i = 0; i < a_matter->numberOfAtoms(); i++) {
    atomnrs.push_back(a_matter->getAtomicNr(i));
    atoms_config.atomicNrs[i] = a_matter->getAtomicNr(i);
    atoms_config.is_frozen[i] = a_matter->getFixed(i);
    atoms_config.id[i] = i + 1;
  }

  gpr::Field<double> falseConDat = generateAtomsConfigField(*a_matter);
  atoms_config.assignFromField(falseConDat);

  unique_atomtypes = std::set<int>(atomnrs.begin(), atomnrs.end());
  n_at = unique_atomtypes.size();
  fake_atype = 0;
  for (auto uatom : unique_atomtypes) {
    atype_to_gprd_atype.insert(
        std::pair<int, int>(static_cast<int>(uatom), fake_atype));
    fake_atype++;
  }
  // for (const auto& pair : atype_to_gprd_atype) {
  //     fmt::print("{}: {}\n", pair.first, pair.second);
  // }

  number_of_mov_atoms = atoms_config.countMovingAtoms();
  number_of_fro_atoms = atoms_config.is_frozen.getSize() - number_of_mov_atoms;

  // fmt::print("Number of moving atoms: {}\n", number_of_mov_atoms);
  // fmt::print("Number of frozen atoms: {}\n", number_of_fro_atoms);

  if (number_of_fro_atoms > 0 && number_of_mov_atoms > 0) {
    // Resize structures for moving and frozen atoms
    atoms_config.atoms_mov.resize(number_of_mov_atoms);
    atoms_config.atoms_froz_inactive.resize(number_of_fro_atoms);

    //!> Does a horrible hack to ensure that this is filled correctly.
    //! Essentially we
    //! use the Map of <EON atomtype, GPR faketype> to generate the fully filled
    //! vectors for moving and frozen_inactive
    //! FIXME: We should really just use the EON atomtype everywhere
    if (atype_to_gprd_atype.size() > 1) {
      int mov_counter = 0;
      int froz_inactive_counter = 0;
      for (auto i = 0; i < a_matter->numberOfAtoms(); i++) {
        if (a_matter->getFixed(i)) {
          //!> Is a fixed atom
          //!> Use EON's atomtype as a key for the GPR's fake atomtype
          atoms_config.atoms_froz_inactive.type[froz_inactive_counter] =
              atype_to_gprd_atype.at(atomnrs[i]);
          froz_inactive_counter++;
        } else {
          //!> Is moving
          atoms_config.atoms_mov.type[mov_counter] =
              atype_to_gprd_atype.at(atomnrs[i]);
          mov_counter++;
        }
      }
    }
    //!> Special case when there's only one atom type, we can now just use the
    //!`set` function of the `Field`. Essentially now we only have one atom type
    else if (atype_to_gprd_atype.size() == 1) {
      atoms_config.atoms_mov.type.set(atype_to_gprd_atype.at(atomnrs[0]));
      atoms_config.atoms_froz_inactive.type.set(
          atype_to_gprd_atype.at(atomnrs[0]));
    }
    // Assign moving and frozen atoms and list all frozen atoms as inactive
    gpr::Index_t counter_f = 0, counter_m = 0;
    for (gpr::Index_t n = 0; n < atoms_config.is_frozen.getSize(); ++n) {
      if (atoms_config.is_frozen[n] == MOVING_ATOM)
        atoms_config.atoms_mov.positions.set(0, counter_m++,
                                             atoms_config.positions.at(n));
      else
        atoms_config.atoms_froz_inactive.positions.set(
            0, counter_f++, atoms_config.positions.at(n));
    }
    //!> End case where we have both nonzero moving and nonzero frozen atoms
  } else {
    if (number_of_mov_atoms == 0) {
      //!> Sometimes, nothing happens
      SPDLOG_CRITICAL(
          " You need to have atoms move!!!\nIn stillness there is only "
          "death\n");
      ::exit(1);
    }
    //!> Now we will consider the case when everything is moving
    //! Everything is almost exactly the same, only we don't have frozen atoms
    atoms_config.atoms_mov.resize(number_of_mov_atoms);

    //!> FIXME: Same caveats as documented above
    if (atype_to_gprd_atype.size() > 1) {
      int mov_counter = 0;
      for (auto i = 0; i < a_matter->numberOfAtoms(); i++) {
        atoms_config.atoms_mov.type[mov_counter] =
            atype_to_gprd_atype.at(atomnrs[i]);
        mov_counter++;
      }
    }
    //!> Special case when there's only one atom type, we can now just use the
    //!`set` function of the `Field`. Essentially now we only have one atom type
    else if (atype_to_gprd_atype.size() == 1) {
      atoms_config.atoms_mov.type.set(atype_to_gprd_atype.at(atomnrs[0]));
    }
  }
  // Pairtype indices for pairs of atomtypes (n_at x n_at)
  // Active pairtypes are indexed as 0,1,...,n_pt-1. Inactive pairtypes are
  // given index EMPTY.
  atoms_config.n_pt = n_at;
  atoms_config.pairtype.resize(n_at, n_at);
  atoms_config.pairtype.set(EMPTY);

  // fmt::print("We have the following moving atom types");
  // atoms_config.atoms_mov.type.print();

  // Set pairtype indices for moving+moving atom pairs (and update number of
  // active pairtypes)
  problem_setup.setPairtypeForMovingAtoms(
      atoms_config.atoms_mov.type, atoms_config.n_pt, atoms_config.pairtype);

  // atoms_config.pairtype(0, 0) = 0;
  // atoms_config.pairtype(0, 1) = 1;
  // atoms_config.pairtype(1, 0) = 1;
  // atoms_config.pairtype(1, 1) = -1;
  // fmt::print("We have the following pairtypes indices");
  // atoms_config.pairtype.print();

  // // Activate frozen atoms within activation distance
  //   // This is now in AtomicGPDimer.cpp
  // problem_setup.activateFrozenAtoms(R_init, parameters.actdist_fro.value,
  //                                 atoms_config);

  return atoms_config;
}

gpr::Field<double> generateAtomsConfigField(const Matter &mat) {
  const size_t natoms = mat.numberOfAtoms();
  gpr::Field<double> conf_prim; // always takes a 5 membered field
  conf_prim.resize(natoms, 5);
  gpr::EigenMatrix positions = mat.getPositions();
  positions.conservativeResize(positions.rows(), 5);
  for (size_t idx{0}; idx < natoms; idx++) {
    positions(idx, 3) = mat.getFixed(idx);
    positions(idx, 4) = idx;
  }
  for (size_t idx{0}; idx < positions.size(); idx++) {
    conf_prim.getInternalVector()[idx] =
        positions.reshaped<Eigen::RowMajor>()[idx];
  }
  return conf_prim;
}

gpr::Observation eon_matter_to_init_obs(shared_ptr<Matter> a_matter) {
  gpr::Observation o;
  o.clear();
  o.R.resize(a_matter->getPositions().rows(), a_matter->getPositions().cols());
  o.G.resize(a_matter->getForces().rows(), a_matter->getForces().cols());
  o.E.resize(1);
  o.E.set(a_matter->getPotentialEnergy());
  o.R.assignFromEigenMatrix(a_matter->getPositions());
  o.G.assignFromEigenMatrix(a_matter->getForces());
  return o;
}

std::pair<gpr::AtomsConfiguration, gpr::Coord>
eon_matter_to_frozen_conf_info(std::shared_ptr<Matter> a_matter,
                               double a_activeRadius) {
  gpr::Coord R_init;
  aux::ProblemSetUp problem_setup;
  auto retconf = eon_matter_to_atmconf(a_matter);
  gpr::EigenMatrix freePos = a_matter->getPositionsFree();
  R_init.resize(1, freePos.size());
  for (size_t idx{0}; idx < freePos.size(); ++idx) {
    R_init(0, idx) = freePos.reshaped<Eigen::RowMajor>()[idx];
  }
  problem_setup.activateFrozenAtoms(R_init, a_activeRadius, retconf);

  return make_pair(retconf, R_init);
}
} // namespace helpers::gproptim::input
