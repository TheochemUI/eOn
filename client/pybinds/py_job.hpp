#ifndef PY_JOB_H
#define PY_JOB_H

#include "../Job.h"
#include "py_wrapper.hpp"

template <class JobBase = Job>
class PyJob : public JobBase {
public:
    using JobBase::JobBase;
    std::vector<std::string> run() override {
        PYBIND11_OVERRIDE_PURE(std::vector<std::string>, JobBase, run, );
    };
};

#endif /* PY_JOB_H */
