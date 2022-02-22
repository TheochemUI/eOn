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

const std::string Job::PROCESS_SEARCH =           "process_search"s;
const std::string Job::SADDLE_SEARCH =            "saddle_search"s;
const std::string Job::MINIMIZATION =             "minimization"s;
const std::string Job::POINT =                    "point"s;
const std::string Job::PARALLEL_REPLICA =         "parallel_replica"s;
const std::string Job::REPLICA_EXCHANGE =         "replica_exchange"s;
const std::string Job::BASIN_HOPPING =            "basin_hopping"s;
const std::string Job::HESSIAN =                  "hessian"s;
const std::string Job::FINITE_DIFFERENCE =        "finite_difference"s;
const std::string Job::NUDGED_ELASTIC_BAND =      "nudged_elastic_band"s;
const std::string Job::DYNAMICS =                 "molecular_dynamics"s;
const std::string Job::SAFE_HYPER =               "safe_hyper"s;
const std::string Job::TAD =                      "tad"s;
const std::string Job::PREFACTOR =                "prefactor"s;
const std::string Job::GLOBAL_OPTIMIZATION =      "global_optimization"s;
const std::string Job::STRUCTURE_COMPARISON =     "structure_comparison"s;
const std::string Job::MONTE_CARLO =              "monte_carlo"s;
// const std::string Job::TEST =                     "test"s;

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
    // }else if (parameters->job == Job::TEST) {
    //     job =  new TestJob(parameters);
    }else if (parameters->job == Job::MONTE_CARLO) {
        job = new MonteCarloJob(parameters);
    }

    return job;
}
