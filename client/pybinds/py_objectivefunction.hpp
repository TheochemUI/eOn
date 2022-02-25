#ifndef PY_OBJECTIVEFUNCTION_H
#define PY_OBJECTIVEFUNCTION_H

#include "../ObjectiveFunction.h"
#include "py_wrapper.hpp"

template <class ObjectiveFunctionBase = ObjectiveFunction>
class PyObjectiveFunction : public ObjectiveFunctionBase {
public:
    using ObjectiveFunctionBase::ObjectiveFunctionBase;
        double getEnergy() override {
        PYBIND11_OVERRIDE_PURE(double, ObjectiveFunctionBase, getEnergy, );
        };
        double getConvergence() override {
        PYBIND11_OVERRIDE_PURE(double, ObjectiveFunctionBase, getConvergence, );
        };
        bool isConverged() override {
        PYBIND11_OVERRIDE_PURE(bool, ObjectiveFunctionBase, isConverged,);
        };
        int degreesOfFreedom() override {
        PYBIND11_OVERRIDE_PURE(int, ObjectiveFunctionBase, degreesOfFreedom,);
        };
        VectorXd getPositions() override {
        PYBIND11_OVERRIDE_PURE(VectorXd, ObjectiveFunctionBase, getPositions,);
        };
        VectorXd getGradient(bool fdstep=false) override {
        PYBIND11_OVERRIDE_PURE(VectorXd, ObjectiveFunctionBase, getGradient, fdstep=false);
        };
        VectorXd difference(VectorXd a, VectorXd b) override {
        PYBIND11_OVERRIDE_PURE(VectorXd, ObjectiveFunctionBase, difference, a, b);
        };
        void setPositions(VectorXd x) override {
            PYBIND11_OVERRIDE_PURE(void, ObjectiveFunctionBase, x);
        };
};



#endif /* PY_OBJECTIVEFUNCTION_H */
