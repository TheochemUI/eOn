# eOn job request/result envelope — Cap'n Proto L0 for kill-file-IPC.
#
# Geometry is flat C-order positions (len=3N) and box (len=9), not .con text.
# ConFrame remains the on-disk codec; this is the runtime / wire contract.
# Wire ordinals are API: never renumber existing fields.
#
# statusCode preserves historical results.dat integer termination_reason values.
# statusText is human-readable (describeStatus).

@0x9e4f1c2a8b7d6301;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("eonc::job_ssot");

# ---------------------------------------------------------------------
# Geometry: one configuration (reactant / saddle / product / image)
# ---------------------------------------------------------------------
struct Geometry {
  # Cartesian positions, C-order, length 3 * nAtoms
  positions @0 :List(Float64);
  # Cell matrix, row-major, length 9 (same convention as Matter / Structure)
  box @1 :List(Float64);
  atomicNumbers @2 :List(UInt16);
  # Empty list means all free; otherwise length nAtoms (true = frozen)
  frozen @3 :List(Bool);
  # Optional chemical symbols (when Z alone is insufficient for I/O)
  symbols @4 :List(Text);
  masses @5 :List(Float64);
  # Optional energy associated with this geometry (eV)
  energy @6 :Float64 = 0.0;
  hasEnergy @7 :Bool = false;
}

# ---------------------------------------------------------------------
# Job request (server → mechanics engine)
# ---------------------------------------------------------------------
struct JobRequest {
  jobId @0 :Text;
  # Matches Results.dat / Parameters job names (process_search, minimization, …)
  jobType @1 :Text;
  # Optional packed Parameters (JSON UTF-8) when not sharing an in-process object
  paramsJson @2 :Text;
  reactant @3 :Geometry;
  # Optional saddle seed / mode (3N)
  mode @4 :List(Float64);
  displacement @5 :Geometry;
  # NEB endpoints / path seeds
  product @6 :Geometry;
  images @7 :List(Geometry);
}

# ---------------------------------------------------------------------
# Force-call accounting (union of job writers)
# ---------------------------------------------------------------------
struct ForceCalls {
  total @0 :UInt64 = 0;
  minimization @1 :UInt64 = 0;
  saddle @2 :UInt64 = 0;
  prefactors @3 :UInt64 = 0;
  neb @4 :UInt64 = 0;
  dephase @5 :UInt64 = 0;
  dynamics @6 :UInt64 = 0;
  refine @7 :UInt64 = 0;
  sampling @8 :UInt64 = 0;
}

# ---------------------------------------------------------------------
# Job result (mechanics → server)
# ---------------------------------------------------------------------
struct ScalarExtra {
  key @0 :Text;
  value @1 :Float64;
}

struct JobResult {
  jobId @0 :Text;
  jobType @1 :Text;
  # Historical results.dat integer (e.g. MinMode status / RunStatus)
  statusCode @2 :Int32 = 0;
  statusText @3 :Text = "good";
  potentialType @4 :Text;
  randomSeed @5 :Int64 = -1;

  forceCalls @6 :ForceCalls;
  potentialEnergy @7 :Float64 = 0.0;
  potentialEnergySaddle @8 :Float64 = 0.0;
  potentialEnergyReactant @9 :Float64 = 0.0;
  potentialEnergyProduct @10 :Float64 = 0.0;
  barrierReactantToProduct @11 :Float64 = 0.0;
  barrierProductToReactant @12 :Float64 = 0.0;
  prefactorReactantToProduct @13 :Float64 = 0.0;
  prefactorProductToReactant @14 :Float64 = 0.0;
  displacementSaddleDistance @15 :Float64 = 0.0;

  # Optional dynamics saddle fields
  simulationTime @16 :Float64 = 0.0;
  mdTemperature @17 :Float64 = 0.0;
  hasDynamics @18 :Bool = false;

  reactant @19 :Geometry;
  saddle @20 :Geometry;
  product @21 :Geometry;
  minimized @22 :Geometry;
  # Mode at saddle (3N) or empty
  mode @23 :List(Float64);
  # NEB intermediate images (excluding or including ends — job-defined)
  path @24 :List(Geometry);

  wallTimeSeconds @25 :Float64 = 0.0;
  userTimeSeconds @26 :Float64 = 0.0;
  systemTimeSeconds @27 :Float64 = 0.0;
  extras @28 :List(ScalarExtra);
  clientVersion @29 :Text;
}
