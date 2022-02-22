#include "py_wrapper.hpp"

PYBIND11_MODULE(eonclient, m) {
    py_parameters(m);
    py_potential(m);
}
