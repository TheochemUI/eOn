#ifndef PY_OPTIMIZER_H
#define PY_OPTIMIZER_H

#include "../Optimizer.h"
#include "py_wrapper.hpp"

template <class OptimizerBase = Optimizer>
class PyOptimizer : public OptimizerBase {
public:
    using OptimizerBase::OptimizerBase;
    int step(double maxMove) override {
        PYBIND11_OVERRIDE_PURE(int, OptimizerBase, step, maxMove);
    };
    int run(int maxIterations, double maxMove) override {
        PYBIND11_OVERRIDE_PURE(int, OptimizerBase, run, maxIterations, maxMove);
    };
};

#endif /* PY_OPTIMIZER_H */
