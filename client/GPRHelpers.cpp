#include "GPRHelpers.h"

#include <map>
#include <unordered_map>
#include <set>
#include <cassert>

gpr::InputParameters
helper_functions::eon_parameters_to_gprd(Parameters *parameters) {
  gpr::InputParameters p;
  // Problem parameters
  p.actdist_fro.value = parameters->gprActiveRadius;
  p.dimer_sep.value = parameters->gprDimerSep;
  p.method_rot.value = parameters->gprDimerRotOptMethod;
  p.method_trans.value = parameters->gprDimerTransOptMethod;
  p.param_trans.value[0] = parameters->gprDimerConvStep;
  p.param_trans.value[1] = parameters->gprDimerMaxStep;
  p.T_dimer.value = parameters->gprForceThreshold;
  p.initrot_nogp.value = parameters->gprDimerInitRotGP;
  p.T_anglerot_init.value = parameters->gprDimerConvergedAngle;
  p.num_iter_initrot.value = parameters->gprDimerInitRotationsMax;
  p.inittrans_nogp.value = parameters->gprDimerInitTransGP;
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
  p.check_derivative.value = parameters->gprOptCheckDerivatives;
  p.max_iter.value = parameters->gprOptMaxIterations;
  p.tolerance_func.value = parameters->gprOptTolFunc;
  p.tolerance_sol.value = parameters->gprOptTolSol;
  p.lambda_limit.value = parameters->gprOptLambdaLimit;
  p.lambda.value = parameters->gprOptLambdaInit;
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

gpr::InputParameters
helper_functions::eon_parameters_to_gprpot(Parameters *parameters) {
  gpr::InputParameters p;
  // Problem parameters
  p.actdist_fro.value = parameters->gprPotActiveRadius;
  // GPR Parameters
  p.gp_sigma2.value = parameters->gprPotSigma2;
  p.jitter_sigma2.value = parameters->gprPotJitterSigma2;
  p.sigma2.value = parameters->gprPotNoiseSigma2;
  p.prior_mu.value = parameters->gprPotPriorMu;
  p.prior_s2.value = parameters->gprPotPriorSigma2;
  p.prior_nu.value = parameters->gprPotPriorNu;
  p.check_derivative.value = parameters->gprPotOptCheckDerivatives;
  p.max_iter.value = parameters->gprPotOptMaxIterations;
  p.tolerance_func.value = parameters->gprPotOptTolFunc;
  p.tolerance_sol.value = parameters->gprPotOptTolSol;
  p.lambda_limit.value = parameters->gprPotOptLambdaLimit;
  p.lambda.value = parameters->gprPotOptLambdaInit;
  p.magnSigma2.value = parameters->gprPotmagnSigma2;
  p.constSigma2.value = parameters->gprPotconstSigma2;
  return p;
}
// FIXME: Take in the active / inactive pairs / atomtypes
gpr::AtomsConfiguration helper_functions::eon_matter_to_atmconf(Matter *matter) {
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
  int fake_atype; //!> False "atomtype" for GPR Dimer

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

    //!> Does a horrible hack to ensure that this is filled correctly. Essentially we
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
      log(" You need to have atoms move!!!\nIn stillness there is only "
          "death\n");
      exit(1);
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

  // Set pairtype indices for moving+moving atom pairs (and update number of
  // active pairtypes)
  problem_setup.setPairtypeForMovingAtoms(
      atoms_config.atoms_mov.type, atoms_config.n_pt, atoms_config.pairtype);

  // // Activate frozen atoms within activation distance
  // problem_setup.activateFrozenAtoms(R_init, parameters.actdist_fro.value,
  //                                 atoms_config);

  return atoms_config;
}

gpr::Observation helper_functions::eon_matter_to_init_obs(Matter& matter) {
  gpr::Observation obs;
  // TODO: Moving is Free, but also needs to have Frozen Active
  // Resizing can be only on free positions at the moment,
  // They are the same
  size_t nfree(matter.numberOfFreeAtoms()), dimens{3};
  obs.R.resize(1, nfree * dimens);
  obs.G.resize(1, nfree * dimens);
  obs.E.resize(1);
  auto pe_forces{matter.maybe_cached_energy_forces_free()};
  obs.E.set(std::get<double>(pe_forces));
  AtomMatrix freePos = matter.getPositionsFree();
  AtomMatrix freeForces = std::get<AtomMatrix>(pe_forces);
  for (size_t idx{0}; idx < nfree * dimens; ++idx){
      obs.R[idx] = freePos.data()[idx];
      obs.G[idx] = -1 * freeForces.data()[idx];
  }
  return obs;
}

std::pair<double, AtomMatrix> helper_functions::energy_and_forces(Matter *matter, Potential *pot){
  // TODO: Sanity checks, if the dynamic_cast fails it is a nullptr and so is
  // not evaluated
  // if (!pot->getName().compare("gpr_pot"s)){
  //   return helper_functions::gpr_energy_and_forces(matter, static_cast<GPRPotential*>(pot));
  // }
  // std::cout<<"Calling "<<pot->getName()<<std::endl;
  // if (auto* gppot = dynamic_cast<GPRPotential*>(pot)){
  //   return helper_functions::gpr_energy_and_forces(matter, gppot);
  // }
  int nAtoms = matter->numberOfAtoms();
  auto posdata = matter->getPositions();
  auto celldat = matter->getCell();
  AtomMatrix forces = AtomMatrix::Constant(nAtoms, 3, 0);
  double *pos = posdata.data();
  double *frcs = forces.data();
  double *bx = celldat.data();
  double energy{0};
  pot->force(nAtoms, pos, nullptr, frcs, &energy, bx, 1);
  AtomMatrix finForces{forces};
  for (int i = 0; i <nAtoms; i++){
    if(matter->getFixed(i)){
      finForces.row(i).setZero();
    }
  }
  return std::make_pair(energy, finForces);
}

std::pair<double, AtomMatrix> helper_functions::energy_and_forces_free(Matter *matter, Potential *pot){
  int nAtoms = matter->numberOfFreeAtoms();
  auto posdata = matter->getPositionsFree();
  auto celldat = matter->getCell();
  AtomMatrix forces = AtomMatrix::Constant(nAtoms, 3, 0);
  double *pos = posdata.data();
  double *frcs = forces.data();
  double *bx = celldat.data();
  double energy{0};
  pot->force(nAtoms, pos, nullptr, frcs, &energy, bx, 1);
  return std::make_pair(energy, forces);
}

std::pair<double, AtomMatrix> helper_functions::gpr_energy_and_forces(Matter *matter, GPRPotential *gprpot){
  int nAtoms = matter->numberOfAtoms();
  auto posdat = matter->getPositions();
  auto celldat = matter->getCell();
  AtomMatrix forces = AtomMatrix::Constant(nAtoms, 3, 0);
  auto calcEF = gprpot->force(posdat, matter->getAtomicNrs(), celldat, 1);
  double potentialEnergy = std::get<double>(calcEF);
  auto finForces = std::get<AtomMatrix>(calcEF);
  for (int i = 0; i < nAtoms; i++){
    if(matter->getFixed(i)){
      finForces.row(i).setZero();
    }
  }
  return std::make_pair(potentialEnergy, finForces);
}

std::pair<gpr::AtomsConfiguration, gpr::Coord> helper_functions::eon_matter_to_frozen_conf_info(Matter *matter, double activeRadius){
  gpr::Coord R_init;
  aux::ProblemSetUp problem_setup;
  auto retconf = helper_functions::eon_matter_to_atmconf(matter);
  R_init.resize(1, 3 * matter->numberOfFreeAtoms());
  int counter = 0;
  for(int i = 0; i < matter->getPositionsFree().rows(); ++i) {
    for(int j = 0; j < matter->getPositionsFree().cols(); ++j) {
      R_init[counter++] = matter->getPositionsFree()(i, j);
    }
  }
  problem_setup.activateFrozenAtoms(R_init, activeRadius,
                                    retconf);

  return std::make_pair(retconf, R_init);
}

void helper_functions::MatterHolder::getEnergyGradient(const Eigen::VectorXd& w,
                                   const gpr::EigenMatrix& x,
                                   const Eigen::VectorXd& x_ind,
                                   const Eigen::VectorXd& y,
                                   gpr::EnergyAndGradient& energy_and_gradient){
  if (this->mrr == nullptr){
    throw("Can't handle incorrectly initialized MatterHolder");
    // return false;
  }
  double erg{mrr->getPotentialEnergy()};
  Eigen::VectorXd grd = mrr->getForcesV();
  energy_and_gradient.energy = &erg;
  energy_and_gradient.gradient = &grd;
  // return true;
  }

std::vector<Matter> helper_functions::prepInitialPath(
           Parameters *params,
           std::string fname_reactant,
           std::string fname_product){
  // Prep final, initial images
  auto reactantFilename = helper_functions::getRelevantFile(fname_reactant);
  auto productFilename = helper_functions::getRelevantFile(fname_product);
  Matter initmatter(params), finalmatter(params);
  initmatter.con2matter(reactantFilename);
  finalmatter.con2matter(productFilename);
  initmatter.getPotentialEnergy();
  finalmatter.getPotentialEnergy();
  // Setup path
  const int natoms = initmatter.numberOfAtoms();
  const int nimages = params->nebImages;
  const int totImages = nimages + 2; // Final and end
  std::vector<Matter> imageArray;
  for (size_t idx{0}; idx < totImages; idx++){
    imageArray.emplace_back(initmatter);
  }
    imageArray.back() = finalmatter;
    auto posInit = imageArray.front().getPositionsFree();
    auto posFinal = imageArray.back().getPositionsFree();
    const AtomMatrix imageSep = imageArray.front().pbc(posFinal-posInit)/(imageArray.size());
    for (double idx{0}; auto &image: imageArray){
      if (idx == imageArray.size()-1 or (idx == 0)) { // Don't change the final and first image
        ++idx;
        continue;
      }
      image.setPositionsFree(posInit + imageSep * idx);
      image.setPotential(initmatter.getPotential());
      image.getPotentialEnergy();
      image.getForcesFree();
      image.useCache = true;
      ++idx;
    }
    return imageArray;
}

gpr::Observation helper_functions::prepInitialObs(std::vector<Matter> &vecmat) {
  gpr::Observation resObs;
  for (auto& mat : vecmat){
    resObs.append(helper_functions::eon_matter_to_init_obs(mat));
  }
  return resObs;
}

bool helper_functions::maybeUpdateObs(NudgedElasticBand& neb, gpr::Observation& prevObs, Parameters& params){
  // To prevent double calculations, we track the indices being compared
  // If all forces are within the convergence --> no update
  // If some images are within the convergence --> mark them as cachable (1)
  //   Then append the rest (2)
  bool updated{false};
  double fmax{0};
  // TODO: Handle climbingImage
  // TODO: Handle different convergence measures
  double nebConvergedForce {params.nebConvergedForce/10}, trupotdiff{0.0};
  auto potential = Potential::getPotential(&params);
  for (long idx {0}; idx <= neb.images+1; idx++){// excludes final, initial
    auto true_energy_forces = helper_functions::energy_and_forces(neb.image[idx], potential);
    trupotdiff = (neb.image[idx]->getForces() - std::get<AtomMatrix>(true_energy_forces)).norm();
    // std::cout<<trupotdiff<<std::endl;
    if (trupotdiff > nebConvergedForce){
      updated = true;
      neb.image[idx]->getPotential(); // Make sure energy is present
      neb.image[idx]->useCache = true; // For later
    };
  }
  if (updated){
    for (long idx {0}; idx <= neb.images+1; idx++){// prepare observation
      // std::cout<<"Appending to image "<<idx<<"\n";
      prevObs.append(helper_functions::eon_matter_to_init_obs(*neb.image[idx]));
      // prevObs.printSizes();
    }
  }
  return updated;
}

std::unique_ptr<NudgedElasticBand> helper_functions::prepGPRNEBround(GPRPotential& trainedGPR, Matter& reactant, Matter& product, Parameters& params){
  reactant.setPotential(&trainedGPR);
  product.setPotential(&trainedGPR);
  auto neb = std::make_unique<NudgedElasticBand>(dynamic_cast<Matter*>(&reactant), dynamic_cast<Matter*>(&product), dynamic_cast<Parameters*>(&params));
  for (long idx {0}; idx <= neb->images+1; idx++){// includes final, initial
    neb->image[idx]->setPotential(&trainedGPR);
  }
  return neb;
}

gpr::GaussianProcessRegression& helper_functions::initializeGPR(gpr::GaussianProcessRegression& gprfunc,
                                                                gpr::AtomsConfiguration& atoms_config,
                                                                gpr::Observation& obsPath,
                                                                std::pair<Parameters, Matter>& eon_matter_params){
    auto initmatter = std::get<Matter>(eon_matter_params);
    auto eonp = std::get<Parameters>(eon_matter_params);
    gpr::GPRSetup gpr_parameters;
    gpr_parameters.jitter_sigma2 = eonp.gprPotJitterSigma2;
    gpr_parameters.sigma2 = eonp.gprPotSigma2;
    gprfunc.setParameters(gpr_parameters);
    auto  potparams = helper_functions::eon_parameters_to_gprpot(&eonp);
    for (int i = 0; i < 9; i++) {
        potparams.cell_dimensions.value[i] = initmatter.getCell()(i);
    }
    gprfunc.initialize(potparams, atoms_config);
    return gprfunc;
}
