#include "GPRHelpers.h"

#include <set>
#include <unordered_map>

gpr::InputParameters
helper_functions::eon_parameters_to_gpr(Parameters *parameters) {
  gpr::InputParameters p;
  // Problem parameters
  p.actdist_fro.value = parameters->gprActiveRadius;
  p.dimer_sep.value = parameters->gprDimerSep;
  p.method_rot.value = parameters->gprDimerRotOptMethod;
  p.method_trans.value = parameters->gprDimerTransOptMethod;
  p.param_trans.value[0] = parameters->gprd_trans_options.step_length;
  p.param_trans.value[1] = parameters->gprd_trans_options.max_step_length;
  p.rotation_removal_projection_threshold.value =
      parameters->gprd_trans_options.rotrem_thresh;
  // Saddle point convergence parameter
  p.T_dimer.value = parameters->saddleConvergedForce;
  p.T_anglerot_init.value = parameters->gprDimerConvergedAngle;
  // ---
  p.initrot_nogp.value = parameters->gprDimerInitRotGP;
  p.num_iter_initrot.value = parameters->gprDimerInitRotationsMax;
  p.T_anglerot_gp.value = parameters->gprDimerRelaxConvAngle;
  p.num_iter_rot_gp.value = parameters->gprDimerRelaxRotationsMax;
  p.divisor_T_dimer_gp.value = parameters->gprDimerDivisorTdimerGP;
  p.disp_max.value = parameters->gprDimerMidpointMaxDisp;
  p.ratio_at_limit.value = parameters->gprDimerRatioAtLimit;
  p.num_bigiter.value = parameters->gprDimerMaxOuterIterations;
  p.num_iter.value = parameters->gprDimerMaxInnerIterations;
  p.islarge_num_iter.value = parameters->gprDimerManyIterations;
  // GPR Parameters
  p.gp_sigma2.value = parameters->gprDimerSigma2;
  p.jitter_sigma2.value = parameters->gprDimerJitterSigma2;
  p.sigma2.value = parameters->gprDimerNoiseSigma2;
  p.prior_mu.value = parameters->gprDimerPriorMu;
  p.prior_nu.value = parameters->gprDimerPriorNu;
  p.prior_s2.value = parameters->gprDimerPriorSigma2;
  // Hyperparameter optimization
  // Common
  p.optimization_alg.value = parameters->gpr_hypopt_options.hopt_method;
  p.check_derivative.value = parameters->gpr_hypopt_options.check_derivative;
  p.max_iter.value = parameters->gpr_hypopt_options.max_iter;
  p.tolerance_func.value = parameters->gpr_hypopt_options.tol_func;
  p.tolerance_sol.value = parameters->gpr_hypopt_options.tol_sol;
  // SCG
  p.lambda_limit.value = parameters->gpr_hypopt_options.scg.lambda_limit;
  p.lambda.value = parameters->gpr_hypopt_options.scg.lambda;
  // ADAM parameters
  p.learning_rate.value = parameters->gpr_hypopt_options.adam.lr;
  p.learning_rate_decay.value = parameters->gpr_hypopt_options.adam.lrd;
  p.beta1.value = parameters->gpr_hypopt_options.adam.b1;
  p.beta2.value = parameters->gpr_hypopt_options.adam.b2;
  p.epsilon.value = parameters->gpr_hypopt_options.adam.eps;
  p.weight_decay.value = parameters->gpr_hypopt_options.adam.weight_decay;
  p.amsgrad.value = parameters->gpr_hypopt_options.adam.amsgrad;
  // Early stopping
  if (parameters->early_stopping_options.dist_metrics == "emd") {
    p.es_dist_metric.value = DistanceMetricType::EMD;
  } else if (parameters->early_stopping_options.dist_metrics == "rmsd") {
    p.es_dist_metric.value = DistanceMetricType::RMSD;
  } else if (parameters->early_stopping_options.dist_metrics == "max1DLog") {
    p.es_dist_metric.value = DistanceMetricType::MAX_1D_LOG;
  }
  p.es_threshold.value = parameters->early_stopping_options.threshold;
  if (parameters->fps_options.metric == "emd") {
    p.fps_metric.value = DistanceMetricType::EMD;
  } else if (parameters->fps_options.metric == "rmsd") {
    p.fps_metric.value = DistanceMetricType::RMSD;
  } else if (parameters->fps_options.metric == "max1DLog") {
    p.fps_metric.value = DistanceMetricType::MAX_1D_LOG;
  }
  p.fps_history.value = parameters->fps_options.history;
  // Prune
  p.use_prune.value = parameters->gprUsePrune;
  p.start_prune_at.value = parameters->gprPruneBegin;
  p.nprune_vals.value = parameters->gprPruneNVals;
  p.prune_threshold.value = parameters->gprPruneThreshold;
  // Debugging
  p.report_level.value = parameters->gprReportLevel;
  p.debug_level.value = parameters->gprDebugLevel;
  p.debug_output_dir.value = parameters->gprDebugOutDir;
  p.debug_output_file_R.value = parameters->gprDebugPosFile;
  p.debug_output_file_E.value = parameters->gprDebugEnergyFile;
  p.debug_output_file_G.value = parameters->gprDebugGradFile;
  p.debug_output_file_extension.value = parameters->gprDebugOutExt;
  p.debug_offset_from_mid_point.value = parameters->gprDebugOffsetMidPoint;
  p.debug_dy.value = parameters->gprDebugDy;
  p.debug_dz.value = parameters->gprDebugDz;
  return p;
}

// FIXME: Take in the active / inactive pairs / atomtypes
gpr::AtomsConfiguration
helper_functions::eon_matter_to_atmconf(Matter *matter) {
  //   AtomsConfiguration a;
  //   aux::ProblemSetUp problem_setup;
  //   std::vector<int> atomnrs;
  //   a.positions.resize(matter->getPositions().rows(),
  //                      matter->getPositions().cols());
  //   a.is_frozen.resize(matter->numberOfAtoms());
  //   a.id.resize(matter->numberOfAtoms());
  //   a.positions.assignFromEigenMatrix(matter->getPositions());
  //   for (auto i = 0; i < matter->numberOfAtoms(); i++) {
  //     atomnrs.push_back(matter->getAtomicNr(i));
  //     a.is_frozen[i] = matter->getFixed(i);
  //     a.id[i] = i + 1;
  //   }
  //   a.atoms_froz_active.clear();
  //   a.atoms_mov.resize(matter->numberOfFreeAtoms());
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
  //   matter->numberOfFixedAtoms()); for (auto i = 0; i <
  //   matter->numberOfFixedAtoms(); ++i) {
  //     a.atoms_froz_inactive.positions.set(0, i,
  //                               {matter->getPosition(i, 0),
  //                                matter->getPosition(i, 1),
  //                                matter->getPosition(i, 2)});
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
  std::vector<int> atomnrs;
  std::unordered_map<int, int>
      atype_to_gprd_atype; //!> Remember that the atom type in EON is the real
                           //! atomic number, while in GPR Dimer it is a set of
                           //! values from 0 to n-1 so this is EON
  int fake_atype;          //!> False "atomtype" for GPR Dimer

  atoms_config.clear();
  atoms_config.positions.resize(matter->getPositions().rows(),
                                matter->getPositions().cols());
  atoms_config.is_frozen.resize(matter->numberOfAtoms());
  atoms_config.id.resize(matter->numberOfAtoms());
  atoms_config.positions.assignFromEigenMatrix(matter->getPositions());
  atoms_config.atomicNrs.resize(matter->numberOfAtoms());
  for (auto i = 0; i < matter->numberOfAtoms(); i++) {
    atomnrs.push_back(matter->getAtomicNr(i));
    atoms_config.atomicNrs[i] = matter->getAtomicNr(i);
    atoms_config.is_frozen[i] = matter->getFixed(i);
    atoms_config.id[i] = i + 1;
  }

  unique_atomtypes = std::set<int>(atomnrs.begin(), atomnrs.end());
  n_at = unique_atomtypes.size();
  fake_atype = 0;
  for (auto uatom : unique_atomtypes) {
    atype_to_gprd_atype.insert(
        std::pair<int, int>(static_cast<int>(uatom), fake_atype));
    fake_atype++;
  }

  number_of_mov_atoms = atoms_config.countMovingAtoms();
  number_of_fro_atoms = atoms_config.is_frozen.getSize() - number_of_mov_atoms;

  if (number_of_fro_atoms > 0 && number_of_mov_atoms > 0) {
    // Resize structures for moving and frozen atoms
    atoms_config.atoms_mov.resize(number_of_mov_atoms);
    atoms_config.atoms_froz_inactive.resize(number_of_fro_atoms);

    //!> Does a horrible to ensure that this is filled correctly. Essentially we
    //! use the Map of <EON atomtype, GPR faketype> to generate the fully filled
    //! vectors for moving and frozen_inactive
    //! FIXME: We should really just use the EON atomtype everywhere
    if (atype_to_gprd_atype.size() > 1) {
      int mov_counter = 0;
      int froz_inactive_counter = 0;
      for (auto i = 0; i < matter->numberOfAtoms(); i++) {
        if (matter->getFixed(i)) {
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
      std::exit(1);
    }
    //!> Now we will consider the case when everything is moving
    //! Everything is almost exactly the same, only we don't have frozen atoms
    atoms_config.atoms_mov.resize(number_of_mov_atoms);

    //!> FIXME: Same caveats as documented above
    if (atype_to_gprd_atype.size() > 1) {
      int mov_counter = 0;
      for (auto i = 0; i < matter->numberOfAtoms(); i++) {
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
  atoms_config.pairtype.resize(n_at, n_at);
  atoms_config.pairtype.set(EMPTY);
  atoms_config.n_pt = 0;

  // Set pairtype indices for moving+moving atom pairs (and update number of
  // active pairtypes)
  problem_setup.setPairtypeForMovingAtoms(
      atoms_config.atoms_mov.type, atoms_config.n_pt, atoms_config.pairtype);

  // // Activate frozen atoms within activation distance
  // problem_setup.activateFrozenAtoms(R_init, parameters.actdist_fro.value,
  //                                 atoms_config);

  return atoms_config;
}

gpr::Observation helper_functions::eon_matter_to_init_obs(Matter *matter) {
  gpr::Observation o;
  o.clear();

  const long num_atoms = matter->numberOfAtoms();
  const long num_dims = 3;

  // The GPR library expects a single "image" or configuration per row.
  // Positions and Forces for an N-atom system should be flattened
  // into a 1-row, (N*3)-column matrix.
  o.R.resize(1, num_atoms * num_dims);
  o.G.resize(1, num_atoms * num_dims);
  // The old behavior was
  //  o.R.resize(matter->getPositions().rows(), matter->getPositions().cols());
  //  o.G.resize(matter->getForces().rows(), matter->getForces().cols());
  // The Energy is a single value (1x1) for the whole configuration.
  o.E.resize(1, 1);

  o.R.assignFromEigenMatrix(matter->getPositionsV());
  o.G.assignFromEigenMatrix(matter->getForcesV());
  o.E.set(matter->getPotentialEnergy());
  return o;
}
