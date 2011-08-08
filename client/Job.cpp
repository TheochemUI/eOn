#include "Parameters.h"
#include "Job.h"
#include "ProcessSearchJob.h"
#include "SaddleSearchJob.h"
#include "MinimizationJob.h"
#include "PointJob.h"
#include "HessianJob.h"
#include "ParallelReplicaJob.h"
#include "DistributedReplicaJob.h"
#include "BasinHoppingJob.h"
#include "FiniteDifferenceJob.h"
#include "DimerRotationJob.h"
#include "NudgedElasticBandJob.h"
#include "DynamicsJob.h"
#include "TestJob.h"

const char Job::PROCESS_SEARCH[] =           "process_search";
const char Job::SADDLE_SEARCH[] =            "saddle_search";
const char Job::MINIMIZATION[] =             "minimization";
const char Job::POINT[] =                    "point";
const char Job::PARALLEL_REPLICA[] =         "parallel_replica";
const char Job::DISTRIBUTED_REPLICA[] =      "distributed_replica";
const char Job::BASIN_HOPPING[] =            "basin_hopping";
const char Job::HESSIAN[] =                  "hessian";
const char Job::FINITE_DIFFERENCE[] =        "finite_difference";
const char Job::DIMER_ROTATION[] =           "dimer_rotation";
const char Job::NUDGED_ELASTIC_BAND[] =      "nudged_elastic_band";
const char Job::DYNAMICS[] =                 "md";
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
    }else if (parameters->job == Job::DISTRIBUTED_REPLICA) {
        job =  new DistributedReplicaJob(parameters);
    }else if (parameters->job == Job::BASIN_HOPPING) {
        job =  new BasinHoppingJob(parameters);
    }else if (parameters->job == Job::FINITE_DIFFERENCE) {
        job =  new FiniteDifferenceJob(parameters);
    }else if (parameters->job == Job::DIMER_ROTATION) {
        job =  new DimerRotationJob(parameters);
    }else if (parameters->job == Job::NUDGED_ELASTIC_BAND) {
        job =  new NudgedElasticBandJob(parameters);
    }else if (parameters->job == Job::DYNAMICS) {
        job =  new DynamicsJob(parameters);
    }else if (parameters->job == Job::TEST) {
        job =  new TestJob(parameters);
    }

    return job;
}
