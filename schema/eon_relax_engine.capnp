@0xc7a53f0e91b24d68;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("eonc::params_ssot");

# Wire format of the relax-engine library (libeon_relax_engine.so).
# Ordinals are API. eon_relax_create takes a flat-array of
# RelaxEngineParams. This file stays separate from eon_params.capnp,
# whose snake_case field names the capnp compiler rejects.

struct NebParams {
  imageCount @0 :Int64 = 5;
  maxIterations @1 :Int64 = 1000;
  forceTolerance @2 :Float64 = 0.01;
  minimizeEndpoints @3 :Bool = false;
  climbingImage @4 :Bool = true;
}

struct SaddleParams {
  maxIterations @0 :Int64 = 1000;
  convergedForce @1 :Float64 = 0.01;
}

struct RelaxEngineParams {
  # "neb" or "saddle". Unknown tokens are fail-closed at create.
  kind @0 :Text = "neb";
  neb @1 :NebParams;
  saddle @2 :SaddleParams;
  surfaceEpoch @3 :UInt64 = 0;
  randomSeed @4 :Int64 = -1;
  uncertainty @5 :Float64 = 0.05;
}

