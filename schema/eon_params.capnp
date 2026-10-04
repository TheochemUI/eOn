# eOn simulation parameters — Cap'n Proto L0 field graph (SSoT).
#
# Authoring rule: add user-facing parameter names/types/defaults HERE first.
# Then run: python tools/params_ssot/codegen.py
# Do not invent peer field lists in Parameters.h / eon/schema.py / config.yaml
# for groups covered by this file without regenerating consumers.
#
# A struct tagged `project` is read and written by the generated INI and
# JSON adapters. Structs without that tag stay in the handwritten adapters.
# Users do not author Cap'n Proto binaries for ordinary runs.
#
# Wire ordinals are API: never renumber existing fields.

@0xbc8e0d2f4a71c9e3;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("eonc::params_ssot");

# ---------------------------------------------------------------------
# [Main]
# Defaults match client Parameters::main_options_t (NSDMI) unless noted.
# Server default job=akmc is a server-only dispatcher default (not client).
# ---------------------------------------------------------------------
struct MainOptions {
  # Client default process_search; server config.yaml uses akmc for dispatcher.
  job @0 :Text = "process_search";
  randomSeed @1 :Int64 = -1;
  temperature @2 :Float64 = 300.0;
  quiet @3 :Bool = false;
  writeLog @4 :Bool = true;
  checkpoint @5 :Bool = false;
  iniFilename @6 :Text = "config.ini";
  conFilename @7 :Text = "pos.con";
  finiteDifference @8 :Float64 = 0.01;
  maxForceCalls @9 :Int64 = 0;
  removeNetForce @10 :Bool = true;
  writeConForces @11 :Bool = false;
  # Client-only parallel force evaluation (not in server yaml)
  parallel @12 :Bool = true;
}

# ---------------------------------------------------------------------
# [Potential]
# ---------------------------------------------------------------------
struct PotentialOptions {
  potential @0 :Text = "lj";
  mpiPollPeriod @1 :Float64 = 0.25;
  lammpsLogging @2 :Bool = false;
  lammpsThreads @3 :Int32 = 0;
  emtRasmussen @4 :Bool = false;
  logPotential @5 :Bool = false;
  extPotPath @6 :Text = "./ext_pot";
  potentialsPath @7 :Text = "";
}

# ---------------------------------------------------------------------
# [Structure Comparison]
# ---------------------------------------------------------------------
struct StructureComparisonOptions {
  distanceDifference @0 :Float64 = 0.1;
  neighborCutoff @1 :Float64 = 3.3;
  checkRotation @2 :Bool = false;
  indistinguishableAtoms @3 :Bool = true;
  energyDifference @4 :Float64 = 0.01;
  removeTranslation @5 :Bool = true;
  # Server/python also expose (client match helpers may read server config):
  useCovalent @6 :Bool = false;
  covalentScale @7 :Float64 = 1.3;
  bruteNeighbors @8 :Bool = false;
}

# ---------------------------------------------------------------------
# [Process Search]
# ---------------------------------------------------------------------
struct ProcessSearchOptions {
  minimizeFirst @0 :Bool = true;
  minimizationOffset @1 :Float64 = 0.2;
}

# ---------------------------------------------------------------------
# [Optimizer] (+ nested method knobs)
# ---------------------------------------------------------------------
struct OptimizerLbfgsOptions {
  memory @0 :Int64 = 20;
  inverseCurvature @1 :Float64 = 0.01;
  maxInverseCurvature @2 :Float64 = 0.0;
  autoScale @3 :Bool = true;
  angleReset @4 :Bool = true;
  distanceReset @5 :Bool = true;
  curvature @6 :Text = "reset";
  projectRigid @7 :Bool = false;
  secant @8 :Text = "standard";
  precon @9 :Text = "none";
  h0 @10 :Text = "sy_yy";
  accept @11 :Text = "none";
  extraUpdates @12 :Int64 = 0;
  cautiousEps @13 :Float64 = 0.000001;
  cautiousAlpha @14 :Float64 = 0.01;
  preconA @15 :Float64 = 3.0;
  preconMu @16 :Float64 = 1.0;
  preconRcut @17 :Float64 = 0.0;
  step @18 :Text = "lbfgs";
}

struct OptimizerCgOptions {
  noOvershooting @0 :Bool = false;
  knockOutMaxMove @1 :Bool = false;
  lineSearch @2 :Bool = false;
  lineConverged @3 :Float64 = 0.1;
  lineSearchMaxIter @4 :Int64 = 10;
  maxIterBeforeReset @5 :Int64 = 0;
}

struct OptimizerQuickminOptions {
  steepestDescent @0 :Bool = false;
}

struct OptimizerSdOptions {
  alpha @0 :Float64 = 0.1;
  twoPoint @1 :Bool = false;
}

struct OptimizerXtsciOptions {
  # Engine-local solver. Tokens match xtsci-optimize xts_method_t:
  # lbfgs, bfgs, sr1, sr2, newton, rfo, steepest, adam, pso,
  # polak_ribiere, fletcher_reeves, hestenes_stiefel, dai_yuan,
  # conjugate_descent, hager_zhang, liu_storey, fr_pr,
  # fire, bb, dogleg, fire2
  method @0 :Text = "lbfgs";
  # How an L-BFGS session uses a host Hessian. Tokens match
  # xts_qn_step_t: lbfgs (two-loop, P is H0), newton, rfo.
  qnStep @1 :Text = "lbfgs";
  # Host pair / model Hessian. none, pair, pair_abs, pair_full,
  # exp, c1, lindh, lindh_full, fischer, schlegel, swart.
  # Built in eOn; xtsci only applies it.
  precon @2 :Text = "none";
  # How the session takes a proposed step. xts_accept_t: none, energy,
  # nonmonotone. none is one oracle at the new point.
  accept @3 :Text = "none";
  # HiGHS feasible-set QP on the host Hessian or two-loop direction.
  highs @4 :Bool = false;
  # Embedded manifold. Tokens match xts_manifold_t:
  # euclidean, rigid_quotient, mw_rigid, sphere, so3, stiefel, se3.
  # Isolated molecules: rigid_quotient (Sella Cartesian T+R) or
  # mw_rigid (Page-McIver / Sella IRC Eckart). so3 is length 9;
  # se3 is length 12.
  manifold @5 :Text = "euclidean";
}

struct OptimizerOptions {
  optMethod @0 :Text = "cg";
  convergenceMetric @1 :Text = "norm";
  maxIterations @2 :UInt64 = 1000;
  maxMove @3 :Float64 = 0.2;
  convergedForce @4 :Float64 = 0.01;
  timeStep @5 :Float64 = 1.0;       # input units (fs); runtime may scale by timeUnit
  maxTimeStep @6 :Float64 = 2.5;
  lbfgs @7 :OptimizerLbfgsOptions;
  cg @8 :OptimizerCgOptions;
  quickmin @9 :OptimizerQuickminOptions;
  sd @10 :OptimizerSdOptions;
  xtsci @11 :OptimizerXtsciOptions;
}

# ---------------------------------------------------------------------
# [RgpotPot]
# Names are the C++ members. ini-order is the on-disk chain; the first
# key present wins. overlay-order is the chain on [cpmd], which does not
# read the nwchem_* spellings. The XTBPot and cpmd overlays run after the
# [RgpotPot] read, and only when backend matches.
# ---------------------------------------------------------------------

# section: RgpotPot
# accessor: rgpot_options
# json: RgpotPot
# project: ini,json
# ini-when: potential=RGPOT
# overlay: section=XTBPot when=backend:xtb,xtbpot,gfn,gfnxtb map=xtb_paramset:paramset,xtb_accuracy:accuracy,xtb_electronic_temperature:electronic_temperature,xtb_max_iterations:max_iterations,xtb_uhf:uhf,xtb_charge:charge
# overlay: section=cpmd when=backend:cpmd,cpmdc,cpmdpot fields=functional,cutoff_ry,charge,multiplicity,title,memory_mb,input_block
struct RgpotPotOptions {
  backend @0 :Text = "nwchemc";
  # ini-order: basis, nwchem_basis
  basis @1 :Text = "sto-3g";
  # ini-order: theory, nwchem_theory
  theory @2 :Text = "scf";
  # ini-order: scf_type, nwchem_scf_type
  scf_type @3 :Text = "rhf";
  # ini-order: functional, cpmd_functional
  functional @4 :Text = "BLYP";
  # ini-order: cutOffRy, cutoff_ry, cpmd_cut_off_ry
  cutoff_ry @5 :Float64 = 70.0;
  # ini-order: charge, nwchem_charge
  # overlay-order: charge
  charge @6 :Int32 = 0;
  # ini-order: multiplicity, nwchem_multiplicity
  # overlay-order: multiplicity
  multiplicity @7 :Int32 = 1;
  engine_path @8 :Text = "";
  engine_library @9 :Text = "";
  engine_root @10 :Text = "";
  title @11 :Text = "";
  memory_mb @12 :Int32 = 0;
  scratch_dir @13 :Text = "";
  input_block @14 :Text = "";
  permanent_dir @15 :Text = "";
  params_path @16 :Text = "";
  ranks_per_image @17 :Int32 = 0;
  model_path @18 :Text = "";
  device @19 :Text = "cpu";
  length_unit @20 :Text = "angstrom";
  extensions_directory @21 :Text = "";
  check_consistency @22 :Bool = false;
  uncertainty_threshold @23 :Float64 = -1.0;
  torch_determinism_strict @24 :Bool = false;
  # ini-order: paramset, xtb_paramset
  xtb_paramset @25 :Text = "GFN2xTB";
  # ini-order: accuracy, xtb_accuracy
  xtb_accuracy @26 :Float64 = 1.0;
  # ini-order: electronic_temperature, xtb_electronic_temperature
  xtb_electronic_temperature @27 :Float64 = 300.0;
  # ini-order: max_iterations, xtb_max_iterations
  # cxx-cast: static_cast<int>
  xtb_max_iterations @28 :Int32 = 250;
  # ini-fallback-member: charge
  xtb_charge @29 :Float64 = 0.0;
  # ini-order: uhf, xtb_uhf
  # cxx-cast: static_cast<int>
  xtb_uhf @30 :Int32 = 0;
  # UMA (backend=uma): model task head of the AOTI package.
  task_name @31 :Text = "omol";
}

# ---------------------------------------------------------------------
# Root message: schema-backed parameters document
# ---------------------------------------------------------------------
struct EonParameters {
  main @0 :MainOptions;
  potential @1 :PotentialOptions;
  structureComparison @2 :StructureComparisonOptions;
  processSearch @3 :ProcessSearchOptions;
  optimizer @4 :OptimizerOptions;
  # Schema format version for forward-compat loaders
  schemaVersion @5 :UInt32 = 1;
  rgpot @6 :RgpotPotOptions;
}
