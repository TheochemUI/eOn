#include "Parameters.h"
#include "Job.h"
#include "ProcessSearchJob.h"
#include "SaddleSearchJob.h"
#include "MinimizationJob.h"
#include "PointJob.h"
#include "HessianJob.h"
#include "ParallelReplicaJob.h"
#include "SafeHyperJob.h"
#include "TADJob.h"
#include "ReplicaExchangeJob.h"
#include "BasinHoppingJob.h"
#include "FiniteDifferenceJob.h"
#include "NudgedElasticBandJob.h"
#include "DynamicsJob.h"
#include "PrefactorJob.h"
#include "TestJob.h"
#include "GlobalOptimizationJob.h"
#include "StructureComparisonJob.h"
#include "MonteCarloJob.h"

const char Job::PROCESS_SEARCH[] =           "process_search";
const char Job::SADDLE_SEARCH[] =            "saddle_search";
const char Job::MINIMIZATION[] =             "minimization";
const char Job::POINT[] =                    "point";
const char Job::PARALLEL_REPLICA[] =         "parallel_replica";
const char Job::REPLICA_EXCHANGE[] =         "replica_exchange";
const char Job::BASIN_HOPPING[] =            "basin_hopping";
const char Job::HESSIAN[] =                  "hessian";
const char Job::FINITE_DIFFERENCE[] =        "finite_difference";
const char Job::NUDGED_ELASTIC_BAND[] =      "nudged_elastic_band";
const char Job::DYNAMICS[] =                 "molecular_dynamics";
const char Job::SAFE_HYPER[] =               "safe_hyper";
const char Job::TAD[] =                      "tad";
const char Job::PREFACTOR[] =                "prefactor";
const char Job::GLOBAL_OPTIMIZATION[] =      "global_optimization";
const char Job::STRUCTURE_COMPARISON[] =     "structure_comparison";
const char Job::MONTE_CARLO[] =              "monte_carlo";
const char Job::TEST[] =                     "test";

Job *Job::getJob(Parameters *parameters) {
    Job *job=NULL;
    if (parameters->job == Job::PROCESS_SEARCH) {
        job = new ProcessSearchJob(parameters);
    }else if (parameters->job == Job::SADDLE_SEARCH) {
        job = new SaddleSearchJob(parameters);
    }else if (parameters->job == Job::MINIMIZATION) {
        job = new MinimizationJob(parameters);
    }else if (parameters->job == Job::POINT) {
        job = new PointJob(parameters);
    }else if (parameters->job == Job::HESSIAN) {
        job = new HessianJob(parameters);
    }else if (parameters->job == Job::PARALLEL_REPLICA) {
        job =  new ParallelReplicaJob(parameters);
    }else if (parameters->job == Job::REPLICA_EXCHANGE) {
        job =  new ReplicaExchangeJob(parameters);
    }else if (parameters->job == Job::TAD){
        job = new TADJob(parameters);
    }else if (parameters->job == Job::SAFE_HYPER){
        job = new SafeHyperJob(parameters);
    }else if (parameters->job == Job::BASIN_HOPPING) {
        job =  new BasinHoppingJob(parameters);
    }else if (parameters->job == Job::FINITE_DIFFERENCE) {
        job =  new FiniteDifferenceJob(parameters);
    }else if (parameters->job == Job::NUDGED_ELASTIC_BAND) {
        job =  new NudgedElasticBandJob(parameters);
    }else if (parameters->job == Job::DYNAMICS) {
        job =  new DynamicsJob(parameters);
    }else if (parameters->job == Job::PREFACTOR) {
        job =  new PrefactorJob(parameters);
    }else if (parameters->job == Job::GLOBAL_OPTIMIZATION) {
        job =  new GlobalOptimizationJob(parameters);
    }else if (parameters->job == Job::STRUCTURE_COMPARISON) {
        job =  new StructureComparisonJob(parameters);
    }else if (parameters->job == Job::TEST) {
        job =  new TestJob(parameters);
    }else if (parameters->job == Job::MONTE_CARLO) {
        job = new MonteCarloJob(parameters);
    }

    return job;
}
