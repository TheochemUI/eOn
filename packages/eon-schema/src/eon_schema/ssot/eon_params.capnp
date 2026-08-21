# eOn simulation parameters — Cap'n Proto L0 field graph (SSoT).
#
# Authoring rule: add user-facing parameter names/types/defaults HERE first.
# Then run: python tools/params_ssot/codegen.py
# Do not invent peer field lists in Parameters.h / eon/schema.py / config.yaml
# for groups covered by this file without regenerating consumers.
#
# INI and JSON remain on-disk adapters into the runtime Parameters store;
# users do not need to author Cap'n Proto binaries for ordinary runs.
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
  accept @11 :Text = "energy";
  extraUpdates @12 :Int64 = 0;
  cautiousEps @13 :Float64 = 0.000001;
  cautiousAlpha @14 :Float64 = 0.01;
  preconA @15 :Float64 = 3.0;
  preconMu @16 :Float64 = 1.0;
  preconRcut @17 :Float64 = 0.0;
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
}
