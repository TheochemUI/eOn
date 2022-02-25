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

// const std::string Job::TEST =                     "test"s;

Job *Job::getJob(Parameters *parameters) {
    Job *job=NULL;
    if (parameters->job == JobStrings::PROCESS_SEARCH) {
        job = new ProcessSearchJob(parameters);
    }else if (parameters->job == JobStrings::SADDLE_SEARCH) {
        job = new SaddleSearchJob(parameters);
    }else if (parameters->job == JobStrings::MINIMIZATION) {
        job = new MinimizationJob(parameters);
    }else if (parameters->job == JobStrings::POINT) {
        job = new PointJob(parameters);
    }else if (parameters->job == JobStrings::HESSIAN) {
        job = new HessianJob(parameters);
    }else if (parameters->job == JobStrings::PARALLEL_REPLICA) {
        job =  new ParallelReplicaJob(parameters);
    }else if (parameters->job == JobStrings::REPLICA_EXCHANGE) {
        job =  new ReplicaExchangeJob(parameters);
    }else if (parameters->job == JobStrings::TAD){
        job = new TADJob(parameters);
    }else if (parameters->job == JobStrings::SAFE_HYPER){
        job = new SafeHyperJob(parameters);
    }else if (parameters->job == JobStrings::BASIN_HOPPING) {
        job =  new BasinHoppingJob(parameters);
    }else if (parameters->job == JobStrings::FINITE_DIFFERENCE) {
        job =  new FiniteDifferenceJob(parameters);
    }else if (parameters->job == JobStrings::NUDGED_ELASTIC_BAND) {
        job =  new NudgedElasticBandJob(parameters);
    }else if (parameters->job == JobStrings::DYNAMICS) {
        job =  new DynamicsJob(parameters);
    }else if (parameters->job == JobStrings::PREFACTOR) {
        job =  new PrefactorJob(parameters);
    }else if (parameters->job == JobStrings::GLOBAL_OPTIMIZATION) {
        job =  new GlobalOptimizationJob(parameters);
    }else if (parameters->job == JobStrings::STRUCTURE_COMPARISON) {
        job =  new StructureComparisonJob(parameters);
    // }else if (parameters->job == JobStrings::TEST) {
    //     job =  new TestJob(parameters);
    }else if (parameters->job == JobStrings::MONTE_CARLO) {
        job = new MonteCarloJob(parameters);
    }

    return job;
}
