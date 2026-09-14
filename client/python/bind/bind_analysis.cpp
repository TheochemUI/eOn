/*
** Hessian + Prefactor — first-class vibrational analysis surface.
*/
#include "eigen_numpy.hpp"
#include "eon/Hessian.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Prefactor.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_analysis(nb::module_ &m) {
  using eonc::Hessian;
  using eonc::Matter;
  using eonc::Parameters;

  nb::class_<Hessian>(m, "Hessian",
                      "Finite-difference Hessian / frequencies on a PHVA "
                      "mobile atom set (pass indices; free/fixed is separate)")
      .def(
          "__init__",
          [](Hessian *self, const Parameters &params, Matter &matter) {
            new (self) Hessian(params, &matter);
          },
          nb::arg("parameters"), nb::arg("matter"), nb::keep_alive<1, 2>(),
          nb::keep_alive<1, 3>(),
          "Parameters + Matter (matter pointer kept alive via keep_alive)")
      .def(
          "get_hessian",
          [](Hessian &self, Matter &matter, const NpI64 &atoms) {
            VectorXi at = vectori_from_numpy_i64(atoms);
            MatrixXd H;
            {
              nb::gil_scoped_release release;
              H = self.getHessian(&matter, at);
            }
            return matrix_to_numpy(H);
          },
          nb::arg("matter"), nb::arg("atoms"),
          "Cartesian Hessian for atom indices (int64 1-d). Returns (3n,3n).")
      .def(
          "get_freqs",
          [](Hessian &self, Matter &matter, const NpI64 &atoms) {
            VectorXi at = vectori_from_numpy_i64(atoms);
            VectorXd f;
            {
              nb::gil_scoped_release release;
              f = self.getFreqs(&matter, at);
            }
            return vector_to_numpy(f);
          },
          nb::arg("matter"), nb::arg("atoms"),
          "Mass-weighted frequencies for atom indices (int64 1-d)")
      .def(
          "remove_zero_freqs",
          [](Hessian &self, const NpF64 &freqs) {
            return vector_to_numpy(
                self.removeZeroFreqs(vector_from_numpy(freqs)));
          },
          nb::arg("freqs"), nb::rv_policy::move,
          "Drop near-zero frequencies (Parameters.hessian zero threshold)");

  m.def(
      "get_prefactors",
      [](const Parameters &params, Matter &min1, Matter &saddle, Matter &min2) {
        double pref1 = 0.0, pref2 = 0.0;
        int rc = 0;
        {
          nb::gil_scoped_release release;
          rc = eonc::Prefactor::getPrefactors(params, &min1, &saddle, &min2,
                                              pref1, pref2);
        }
        if (rc != 0)
          throw std::runtime_error(
              "get_prefactors failed (bad Hessian or filter window); rc=" +
              std::to_string(rc));
        return nb::make_tuple(pref1, pref2);
      },
      nb::arg("parameters"), nb::arg("min1"), nb::arg("saddle"),
      nb::arg("min2"),
      "HTST/QQHTST prefactors for reactant min1 → saddle → product min2. "
      "Returns (pref1, pref2).");

  m.def(
      "moved_atoms",
      [](const Parameters &params, Matter &min1, Matter &saddle, Matter &min2) {
        VectorXi atoms;
        if (params.prefactor_options.filter_scheme ==
            eonc::Prefactor::FILTER_FRACTION)
          atoms = eonc::Prefactor::movedAtomsPct(params, &min1, &saddle, &min2);
        else
          atoms = eonc::Prefactor::movedAtoms(params, &min1, &saddle, &min2);
        return vectori_to_numpy(atoms);
      },
      nb::arg("parameters"), nb::arg("min1"), nb::arg("saddle"),
      nb::arg("min2"), "Atom indices selected by prefactor filter (int64 1-d)");

  m.def(
      "all_free_atoms",
      [](Matter &matter) {
        return vectori_to_numpy(eonc::Prefactor::allFreeAtoms(&matter));
      },
      nb::arg("matter"), "Indices of free (unfixed) atoms");
}

} // namespace eonc::pybind
