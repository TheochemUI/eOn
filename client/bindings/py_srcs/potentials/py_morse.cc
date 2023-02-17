// clang-format off
#include "../py_wrapper.hpp"
#include "../py_potential.hpp"
// Binding
#include "../../../potentials/Morse/Morse.h"
#include <utility>
// clang-format on

template <class MorseBase = Morse>
class PyMorse : public PyPotential<MorseBase> {
public:
    using PyPotential<MorseBase>::PyPotential; // Inherit constructor
    // Override pure virtual with non-pure
    std::pair<double, AtomMatrix>
    get_ef(const AtomMatrix positions, const VectorXi atomicNrs, const Matrix3d box) override {
        // This is needed to prevent substitution errors
        using Return = std::pair<double, AtomMatrix>;
        PYBIND11_OVERRIDE(Return, /* Return type */
                          MorseBase, /* Parent class */
                          get_ef, /* Function name (C++) */
                          positions, /* Arg(s) */
                          atomicNrs, box);
    }
};

void py_morse(py::module_ &m) {
    py::class_<Morse, Potential, PyMorse<>,  std::shared_ptr<Morse>>(m, "Morse")
        /*
        ** Constructors
        */
        .def(py::init<Parameters*>())

        /*
        ** Operators
        */

        /*
        ** Methods
        */
        .def("cleanMemory", &Morse::cleanMemory)
        // /** TODO: These force calls should be part of Potential, no need to make separate bindings **/
        /* Technically incorrect, but might be useful, does not account for the fixed / free atoms */
        .def("get_ef", &Morse::get_ef)

        /*
        ** Python helpers
        */

        .def("__repr__", [](const Morse &a) { return "<Morse object>"; });
}
