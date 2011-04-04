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
#include "DisplacementSamplingJob.h"
#include "TestJob.h"

const string Job::PROCESS_SEARCH =           "process_search";
const string Job::SADDLE_SEARCH =            "saddle_search";
const string Job::MINIMIZATION =             "minimization";
const string Job::POINT =                    "point";
const string Job::PARALLEL_REPLICA =         "parallel_replica";
const string Job::DISTRIBUTED_REPLICA =      "distributed_replica";
const string Job::BASIN_HOPPING =            "basin_hopping";
const string Job::HESSIAN =                  "hessian";
const string Job::FINITE_DIFFERENCE =        "finite_difference";
const string Job::DIMER_ROTATION =           "dimer_rotation";
const string Job::DISPLACEMENT_SAMPLING =    "displacement_sampling";
const string Job::TEST =                     "test";

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
    }else if (parameters->job == Job::DISPLACEMENT_SAMPLING) {
        job =  new DisplacementSamplingJob(parameters);
    }else if (parameters->job == Job::TEST) {
        job =  new TestJob(parameters);
    }

    return job;
}
