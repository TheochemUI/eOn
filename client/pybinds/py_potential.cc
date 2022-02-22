#include "py_wrapper.hpp"
#include "../Potential.h"

// See -->

class PyPotential: public Potential {
public:
    /* Inherit the constructors */
    using Potential::Potential;

    /* Trampoline (need one for each virtual function) */
    void initialize() override {
        PYBIND11_OVERRIDE_PURE(
            void, /* Return type */
            Potential,      /* Parent class */
            initialize          /* Name of function in C++ (must match Python name) */
        );
    }
    void force(long nAtoms, const double *positions,
                           const int *atomicNrs, double *forces, double *energy,
                           const double *box, int nImages) override {
        PYBIND11_OVERRIDE_PURE(
            void, /* Return type */
            Potential,      /* Parent class */
            force,          /* Name of function in C++ (must match Python name) */
            nAtoms,     /* Argument(s) */
            positions,
            atomicNrs,
            forces,
            energy,
            box,
            nImages
        );
    }
};


void py_potential(py::module_ &m){
  py::class_<Potential, PyPotential>(m, "Potential")
      .def(py::init())
      /*
      ** Functions
       */
      .def_static("getPotential", &Potential::getPotential)
      .def("initialize", &Potential::initialize)
      // .def("force", py::overload_cast<AtomMatrix>(&Potential::force))
      .def("force", py::overload_cast<long /*nAtoms*/,
           const double* /*positions*/,
           const int* /*atomicNrs*/,
           double* /*forces*/,
           double* /*energy*/,
           const double* /*box*/,
           int /*nImages*/>(&Potential::force))

      /*
      ** Parameters
       */
      // Potential Names

      .def_property_readonly_static("POT_EMT",
                           [](py::object /*self*/){return Potential::POT_EMT;}
                           )

    // .def_readwrite_static("POT_EMT", [](py::object /*self*/){return Potential::POT_EMT;})
      // .def_readwrite_static("POT_EXT", &Potential::POT_EXT)
      // .def_readwrite_static("POT_LJ", &Potential::POT_LJ)
      // .def_readwrite_static("POT_LJCLUSTER", &Potential::POT_LJCLUSTER)
      // .def_readwrite_static("POT_MORSE_PT", &Potential::POT_MORSE_PT)
      // .def_readwrite_static("POT_NEW", &Potential::POT_NEW)
      // Book-keeping
      .def_readwrite_static("fcalls", &Potential::fcalls)
      .def_readwrite_static("fcallsTotal", &Potential::fcallsTotal)
      .def_readwrite_static("wu_fcallsTotal", &Potential::wu_fcallsTotal)
      .def_readwrite_static("totalUserTime", &Potential::totalUserTime)

      // Objects
      .def_readwrite("params", &Potential::params)
      .def_readwrite_static("pot", &Potential::pot)

      /*
      ** Python helpers
       */

      .def("__repr__", [](const Potential &a) {
        return "<Potential object>";
      });

}
