/*
** MinModeSaddleSearch — first-class single-ended saddle search.
** Default non-mutating: run() works on a Matter copy unless inplace=True.
*/
#include "bind_helpers.hpp"
#include "eigen_numpy.hpp"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include <memory>
#include <string>
#include <utility>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_saddle(nb::module_ &m) {
  using eonc::Matter;
  using eonc::MinModeSaddleSearch;
  using eonc::Parameters;
  using eonc::Potential;

  nb::enum_<MinModeSaddleSearch::Status>(m, "SaddleStatus")
      .value("GOOD", MinModeSaddleSearch::STATUS_GOOD)
      .value("INIT", MinModeSaddleSearch::STATUS_INIT)
      .value("BAD_NO_CONVEX", MinModeSaddleSearch::STATUS_BAD_NO_CONVEX)
      .value("BAD_HIGH_ENERGY", MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY)
      .value("BAD_MAX_CONCAVE_ITERATIONS",
             MinModeSaddleSearch::STATUS_BAD_MAX_CONCAVE_ITERATIONS)
      .value("BAD_MAX_ITERATIONS",
             MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS)
      .value("BAD_NOT_CONNECTED", MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED)
      .value("BAD_PREFACTOR", MinModeSaddleSearch::STATUS_BAD_PREFACTOR)
      .value("BAD_HIGH_BARRIER", MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER)
      .value("BAD_MINIMA", MinModeSaddleSearch::STATUS_BAD_MINIMA)
      .value("FAILED_PREFACTOR", MinModeSaddleSearch::STATUS_FAILED_PREFACTOR)
      .value("POTENTIAL_FAILED", MinModeSaddleSearch::STATUS_POTENTIAL_FAILED)
      .value("NONNEGATIVE_ABORT", MinModeSaddleSearch::STATUS_NONNEGATIVE_ABORT)
      .value("NONLOCAL_ABORT", MinModeSaddleSearch::STATUS_NONLOCAL_ABORT)
      .value("NEGATIVE_BARRIER", MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER)
      .value("BAD_MD_TRAJECTORY_TOO_SHORT",
             MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT)
      .value("BAD_NO_NEGATIVE_MODE_AT_SADDLE",
             MinModeSaddleSearch::STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE)
      .value("BAD_NO_BARRIER", MinModeSaddleSearch::STATUS_BAD_NO_BARRIER)
      .value("ZEROMODE_ABORT", MinModeSaddleSearch::STATUS_ZEROMODE_ABORT)
      .value("OPTIMIZER_ERROR", MinModeSaddleSearch::STATUS_OPTIMIZER_ERROR)
      .value("DIMER_LOST_MODE", MinModeSaddleSearch::STATUS_DIMER_LOST_MODE)
      .value("DIMER_RESTORED_BEST",
             MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST);

  // Overloaded so both the typed status properties and the raw int codes
  // returned by run() resolve.
  m.def(
      "saddle_status_message",
      [](MinModeSaddleSearch::Status status) {
        return std::string(MinModeSaddleSearch::statusMessage(status));
      },
      nb::arg("status"), "Human-readable MinModeSaddleSearch status string");

  m.def(
      "saddle_status_message",
      [](int status) {
        return std::string(MinModeSaddleSearch::statusMessage(status));
      },
      nb::arg("status"), "Human-readable MinModeSaddleSearch status string");

  // Free function: Matter-first saddle search with explicit inplace policy.
  m.def(
      "min_mode_saddle_search",
      [](std::shared_ptr<Matter> matter, const NpF64 &mode,
         double reactant_energy, const Parameters &params,
         std::shared_ptr<Potential> pot, bool inplace) {
        auto work = matter_work(matter, inplace);
        AtomMatrix mode_am = atom_matrix_from_numpy(mode);
        MinModeSaddleSearch ss(work, mode_am, reactant_energy, params, pot);
        int status = 0;
        {
          nb::gil_scoped_release release;
          status = ss.run();
        }
        return nb::make_tuple(work,
                              static_cast<MinModeSaddleSearch::Status>(status));
      },
      nb::arg("matter"), nb::arg("mode"), nb::arg("reactant_energy"),
      nb::arg("parameters"), nb::arg("potential"), nb::arg("inplace") = false,
      "Single-ended min-mode saddle search. Returns (Matter, status). "
      "Default copies matter (caller unchanged); inplace=True mutates matter.");

  nb::class_<MinModeSaddleSearch>(
      m, "MinModeSaddleSearch",
      "Single-ended min-mode following saddle search. Prefer "
      "min_mode_saddle_search(..., inplace=) free function. Class form still "
      "binds C++ object; run(inplace=) defaults to non-mutating.")
      .def(
          "__init__",
          [](MinModeSaddleSearch *self, std::shared_ptr<Matter> matter,
             const NpF64 &mode, double reactant_energy,
             const Parameters &params, std::shared_ptr<Potential> pot) {
            AtomMatrix mode_am = atom_matrix_from_numpy(mode);
            // Always construct on a private working copy; run(inplace) decides.
            auto work = std::make_shared<Matter>(*matter);
            new (self)
                MinModeSaddleSearch(std::move(work), mode_am, reactant_energy,
                                    params, std::move(pot));
            // stash original shared_ptr for inplace path via custom holder?
            // Keep simple: store only working copy; run returns that.
          },
          nb::arg("matter"), nb::arg("mode"), nb::arg("reactant_energy"),
          nb::arg("parameters"), nb::arg("potential"), nb::keep_alive<1, 2>(),
          nb::keep_alive<1, 5>(), nb::keep_alive<1, 6>(),
          "Constructs on a **copy** of matter (never mutates the argument). "
          "run() returns the working Matter + status.")
      .def(
          "run",
          [](MinModeSaddleSearch &self) {
            int status = 0;
            {
              nb::gil_scoped_release release;
              status = self.run();
            }
            return status;
          },
          "Run saddle search on the internal working Matter (copy of input). "
          "Returns status int. Working geometry is the bound matter used at "
          "construct (a copy).")
      .def(
          "run_retain_frames",
          [](MinModeSaddleSearch &self, long max_iterations) {
            int status = 0;
            {
              nb::gil_scoped_release release;
              status = self.runRetainFrames(max_iterations);
            }
            return status;
          },
          nb::arg("max_iterations") = -1L,
          "Like run(), but retain climb ConFrames in memory (same stamps as "
          "write_movies climb CON; no disk required). max_iterations=-1 uses "
          "Parameters.saddle_search max.")
      .def(
          "climb_frames",
          [](MinModeSaddleSearch &self) {
            return con_frames_to_python(self.climbFrames());
          },
          "list[readcon.ConFrame] retained by run_retain_frames().")
      .def(
          "take_climb_frames",
          [](MinModeSaddleSearch &self) {
            return con_frames_to_python(self.takeClimbFrames());
          },
          "Take ownership of retained climb frames (clears storage).")
      .def(
          "clear_climb_frames",
          [](MinModeSaddleSearch &self) { self.clearClimbFrames(); },
          "Drop retained climb frames.")
      // getStatus() is int on the SaddleSearchMethod interface; typed here so
      // `search.status == SaddleStatus.GOOD` compares against the same type.
      .def_prop_ro("status",
                   [](const MinModeSaddleSearch &self) {
                     return static_cast<MinModeSaddleSearch::Status>(
                         self.getStatus());
                   })
      .def_prop_ro("status_message",
                   [](const MinModeSaddleSearch &self) {
                     return std::string(
                         MinModeSaddleSearch::statusMessage(self.getStatus()));
                   })
      .def_prop_ro("iteration",
                   [](const MinModeSaddleSearch &self) {
                     return self.getIterationCount();
                   })
      .def_prop_ro(
          "force_calls",
          [](const MinModeSaddleSearch &self) { return self.getForceCalls(); })
      .def_prop_ro(
          "eigenvalue",
          [](MinModeSaddleSearch &self) { return self.getEigenvalue(); })
      .def_prop_ro(
          "eigenvector",
          [](MinModeSaddleSearch &self) {
            return matrix_to_numpy(self.getEigenvector());
          },
          nb::rv_policy::move)
      .def("__repr__", [](const MinModeSaddleSearch &self) {
        return "<MinModeSaddleSearch status=" +
               std::to_string(self.getStatus()) +
               " iter=" + std::to_string(self.getIterationCount()) + ">";
      });
}

} // namespace eonc::pybind
