#include "py_wrapper.hpp"

PYBIND11_MODULE(pyeonclient, m) {
    py_parameters(m);
    py_log(m);
    py_matter(m);
    // Jobs
    py_job(m);
    py_saddlesearchjob(m);
    // Objective Functions
    py_objectivefunction(m);
    py_matterobjfunc(m);
    py_nebobjfunc(m);
    // Nudged Elastic Band
    py_nudgedelasticband(m);
    // Optimizer
    py_optimizer(m);
    py_lbfgs(m);
}
