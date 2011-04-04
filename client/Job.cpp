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

Job *Job::getJob(Parameters *parameters) {
    Job *job=NULL;
    if (parameters->job == Parameters::PROCESS_SEARCH) {
        job = new ProcessSearchJob(parameters);
    }else if (parameters->job == Parameters::SADDLE_SEARCH) {
        job = new SaddleSearchJob(parameters);
    }else if (parameters->job == Parameters::MINIMIZATION) {
        job = new MinimizationJob(parameters);
    }else if (parameters->job == Parameters::POINT) {
        job = new PointJob(parameters);
    }else if (parameters->job == Parameters::HESSIAN) {
        job = new HessianJob(parameters);
    }else if (parameters->job == Parameters::PARALLEL_REPLICA) {
        job =  new ParallelReplicaJob(parameters);
    }else if (parameters->job == Parameters::DISTRIBUTED_REPLICA) {
        job =  new DistributedReplicaJob(parameters);
    }else if (parameters->job == Parameters::BASIN_HOPPING) {
        job =  new BasinHoppingJob(parameters);
    }else if (parameters->job == Parameters::FINITE_DIFFERENCE) {
        job =  new FiniteDifferenceJob(parameters);
    }else if (parameters->job == Parameters::DIMER_ROTATION) {
        job =  new DimerRotationJob(parameters);
    }else if (parameters->job == Parameters::DISPLACEMENT_SAMPLING) {
        job =  new DisplacementSamplingJob(parameters);
    }else if (parameters->job == Parameters::TEST) {
        job =  new TestJob(parameters);
    }

    return job;
}
