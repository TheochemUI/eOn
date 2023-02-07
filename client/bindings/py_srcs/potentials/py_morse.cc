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
    void initialize() override { PYBIND11_OVERRIDE(void, MorseBase, initialize, ); }
    void force(long nAtoms,
               const double *positions,
               const int *atomicNrs,
               double *forces,
               double *energy,
               const double *box,
               int nImages) override {
        PYBIND11_OVERRIDE(void,      /* Return type */
                          MorseBase, /* Parent class */
                          force,     /* Name of function in C++ (must match Python name) */
                          nAtoms,    /* Argument(s) */
                          positions,
                          atomicNrs,
                          forces,
                          energy,
                          box,
                          nImages);
    };
};

void py_morse(py::module_ &m) {
    py::class_<Morse, Potential, PyMorse<>>(m, "Morse")
        /*
        ** Constructors
        */
        .def(py::init<>())
        .def(py::init<double /*re*/, double /*De*/, double /*a*/, double /*cutoff*/>())

        /*
        ** Operators
        */

        /*
        ** Methods
        */
        .def("cleanMemory", &Morse::cleanMemory)
        .def("initialize", &Morse::initialize)
        // /** TODO: These force calls should be part of Potential, no need to make separate bindings **/
        .def(
            "force",
            [](Morse &pot, Matter &mat) {
                // Store and restore
                Potential* oldPot = mat.getPotential();
                mat.setPotential(&pot);
                double e_pot{mat.getPotentialEnergy()};
                AtomMatrix f_pot {mat.getForces()};
                mat.setPotential(oldPot);
                return std::make_pair(e_pot, f_pot);
            },
            py::arg("matter"))
        /* Technically incorrect, but might be useful, does not account for the fixed / free atoms */
        .def(
            "force",
            [](Morse &pot, size_t natoms, AtomMatrix pos, VectorXi atmnrs, Matrix3d box) {
                double e_pot{0};
                AtomMatrix f_pot = Eigen::MatrixXd::Ones(natoms, 3);
                pot.force(natoms, pos.data(), atmnrs.data(), f_pot.data(), &e_pot, box.data(), 1);
                return std::make_pair(e_pot, f_pot);
            },
            py::arg("natoms"),
            "pos"_a,
            "atmnrs"_a,
            "box"_a)
        /*
        ** Parameters
        */

        /*
        ** Python helpers
        */

        .def("__repr__", [](const Morse &a) { return "<Morse object>"; });
}
