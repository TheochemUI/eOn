#include "py_wrapper.hpp"

PYBIND11_MODULE(pyeonclient, m) {
    py_parameters(m);
    py_log(m);
    py_matter(m);
    // Jobs
    py_job(m);
    py_saddlesearchjob(m);
}
