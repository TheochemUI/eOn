#include "bind_helpers.hpp"
#include "eigen_numpy.hpp"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

#include <memory>
#include <stdexcept>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_matter(nb::module_ &m) {
  using eonc::Matter;
  using eonc::Parameters;
  using eonc::Potential;

  nb::class_<Matter>(m, "Matter",
                     "Atomic structure + potential (eOn C++ client core)")
      .def(nb::init<std::shared_ptr<Potential>, const Parameters &>(),
           nb::arg("potential"), nb::arg("parameters"), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>(),
           "Create empty Matter bound to a potential and parameters "
           "(Parameters must outlive Matter)")
      .def(nb::init<const Matter &>(), nb::arg("other"), "Copy construct")
      .def("resize", &Matter::resize, nb::arg("n_atoms"))
      .def_prop_ro("n_atoms", &Matter::numberOfAtoms)
      .def_prop_ro("n_free", &Matter::numberOfFreeAtoms)
      .def_prop_ro("n_fixed", &Matter::numberOfFixedAtoms)

      // positions: zero-copy view via public getPositions() (binding-only cast)
      .def_prop_rw(
          "positions",
          [](Matter &self) {
            return view_n3(matter_positions_ptr(self), self.numberOfAtoms());
          },
          [](Matter &self, const NpF64 &arr) {
            matter_set_positions_buf(self,
                                     require_n3(arr, self.numberOfAtoms()),
                                     self.numberOfAtoms());
          },
          nb::rv_policy::reference_internal,
          "Cartesian positions (n,3) float64 — zero-copy view of C++ storage")
      .def_prop_ro(
          "positions_free",
          [](const Matter &self) {
            return matrix_to_numpy(self.getPositionsFree());
          },
          nb::rv_policy::move)
      .def(
          "set_positions_free",
          [](Matter &self, const NpF64 &arr) {
            self.setPositionsFree(atom_matrix_from_numpy(arr));
          },
          nb::arg("positions"))

      // cell: getCell() is by value — small owned copy; set via Map
      .def_prop_rw(
          "cell",
          [](const Matter &self) { return matrix_to_numpy(self.getCell()); },
          [](Matter &self, const NpF64 &arr) {
            if (arr.ndim() != 2 || arr.shape(0) != 3 || arr.shape(1) != 3) {
              throw std::invalid_argument("cell must be float64 (3,3)");
            }
            matter_set_cell_buf(self, arr.data());
          },
          nb::rv_policy::move, "Cell matrix (3,3) row-major")

      .def_prop_rw(
          "velocities",
          [](const Matter &self) {
            return matrix_to_numpy(self.getVelocities());
          },
          [](Matter &self, const NpF64 &arr) {
            matter_set_velocities_buf(self,
                                      require_n3(arr, self.numberOfAtoms()),
                                      self.numberOfAtoms());
          },
          nb::rv_policy::move)

      // forces: zero-copy via public getForces() / getForcesRaw()
      .def_prop_ro(
          "forces",
          [](Matter &self) {
            {
              nb::gil_scoped_release release;
              (void)self.getForces(); // ensure force cache
            }
            return view_n3(matter_forces_ptr(self), self.numberOfAtoms());
          },
          nb::rv_policy::reference_internal,
          "Forces with fixed atoms zeroed (zero-copy view of cache)")
      .def_prop_ro(
          "forces_raw",
          [](Matter &self) {
            {
              nb::gil_scoped_release release;
              (void)self.getForcesRaw();
            }
            return view_n3(matter_forces_raw_ptr(self), self.numberOfAtoms());
          },
          nb::rv_policy::reference_internal)
      .def_prop_ro(
          "forces_free",
          [](Matter &self) {
            // Force eval under released GIL; numpy conversion must hold GIL.
            AtomMatrix free_forces;
            {
              nb::gil_scoped_release release;
              free_forces = self.getForcesFree();
            }
            return matrix_to_numpy(free_forces);
          },
          nb::rv_policy::move)
      .def(
          "set_forces",
          [](Matter &self, const NpF64 &arr) {
            self.setForces(atom_matrix_from_numpy(arr));
          },
          nb::arg("forces"))

      .def_prop_rw(
          "masses",
          [](const Matter &self) { return vector_to_numpy(self.getMasses()); },
          [](Matter &self, const NpF64 &arr) {
            if (arr.ndim() != 1 ||
                static_cast<long>(arr.shape(0)) != self.numberOfAtoms()) {
              throw std::invalid_argument(
                  "masses must be float64 length n_atoms");
            }
            matter_set_masses_buf(self, arr.data(), self.numberOfAtoms());
          },
          nb::rv_policy::move)
      .def_prop_rw(
          "atomic_numbers",
          [](const Matter &self) {
            return vectori_to_numpy(self.getAtomicNrs());
          },
          [](Matter &self, nb::object obj) {
            auto np = nb::module_::import_("numpy");
            nb::object cont =
                np.attr("ascontiguousarray")(obj, nb::arg("dtype") = "int32");
            auto arr = nb::cast<NpI32>(cont);
            if (arr.ndim() != 1 ||
                static_cast<long>(arr.shape(0)) != self.numberOfAtoms()) {
              throw std::invalid_argument(
                  "atomic_numbers must be length n_atoms");
            }
            matter_set_z_buf(self, reinterpret_cast<const int *>(arr.data()),
                             self.numberOfAtoms());
          },
          nb::rv_policy::move, "Atomic numbers (Z)")
      .def(
          "get_fixed",
          [](const Matter &self, long atom) { return self.getFixed(atom); },
          nb::arg("atom"))
      .def(
          "set_fixed",
          [](Matter &self, long atom, int fixed) {
            self.setFixed(atom, fixed);
          },
          nb::arg("atom"), nb::arg("fixed"))
      .def_prop_rw(
          "fixed",
          [](const Matter &self) {
            const long n = self.numberOfAtoms();
            VectorXi v(n);
            for (long i = 0; i < n; ++i) {
              v(i) = self.getFixed(i);
            }
            return vectori_to_numpy(v);
          },
          [](Matter &self, nb::object obj) {
            auto np = nb::module_::import_("numpy");
            nb::object cont =
                np.attr("ascontiguousarray")(obj, nb::arg("dtype") = "int32");
            auto arr = nb::cast<NpI32>(cont);
            if (arr.ndim() != 1 ||
                static_cast<long>(arr.shape(0)) != self.numberOfAtoms()) {
              throw std::invalid_argument(
                  "fixed mask must be 1-d int array of length n_atoms");
            }
            matter_set_fixed_buf(self,
                                 reinterpret_cast<const int *>(arr.data()),
                                 self.numberOfAtoms());
          },
          nb::rv_policy::move, "Fixed flags (1=fixed, 0=free)")

      // --- free mask ---
      .def_prop_ro(
          "free_mask",
          [](const Matter &self) { return matrix_to_numpy(self.getFree()); },
          nb::rv_policy::move, "Nx3 free mask (1 free, 0 fixed)")

      // --- energy / mechanics ---
      .def_prop_ro(
          "potential_energy",
          [](Matter &self) {
            nb::gil_scoped_release release;
            return self.getPotentialEnergy();
          },
          "Potential energy (GIL released for metatomic/torch autograd)")
      .def_prop_ro("kinetic_energy", &Matter::getKineticEnergy)
      .def_prop_ro("mechanical_energy", &Matter::getMechanicalEnergy)
      .def_prop_ro("energy_variance", &Matter::getEnergyVariance)
      .def_prop_ro(
          "max_force",
          [](Matter &self) {
            nb::gil_scoped_release release;
            return self.maxForce();
          },
          "Max force norm (GIL released for metatomic/torch)")
      .def_prop_ro("force_calls", &Matter::getForceCalls)
      .def("reset_force_calls", &Matter::resetForceCalls)
      .def_prop_ro("needs_force_update", &Matter::needsForceUpdate)
      .def_prop_ro("potential_calls", &Matter::getPotentialCalls)

      // --- PBC ---
      .def_prop_rw("periodic", &Matter::getPeriodic, &Matter::setPeriodic)
      .def_prop_rw("pbc_convention", &Matter::getPbcConvention,
                   &Matter::setPbcConvention)
      .def(
          "pbc",
          [](const Matter &self, const NpF64 &diff) {
            return matrix_to_numpy(self.pbc(atom_matrix_from_numpy(diff)));
          },
          nb::arg("diff"), nb::rv_policy::move,
          "Minimum-image map for displacement array (n,3)")
      .def("apply_periodic_boundary_if_enabled",
           &Matter::applyPeriodicBoundaryIfEnabled)

      // --- potential handle ---
      .def("set_potential", &Matter::setPotential, nb::arg("potential"))
      .def("get_potential", &Matter::getPotential)

      // --- geometry helpers ---
      .def("distance",
           nb::overload_cast<long, long>(&Matter::distance, nb::const_),
           nb::arg("i"), nb::arg("j"))
      .def(
          "distance_to_matter",
          [](const Matter &self, const Matter &other, long index) {
            return self.distance(other, index);
          },
          nb::arg("other"), nb::arg("index"))
      .def("distance_to", &Matter::distanceTo, nb::arg("other"))
      .def("per_atom_norm", &Matter::perAtomNorm, nb::arg("other"))
      .def("compare", &Matter::compare, nb::arg("other"),
           nb::arg("indistinguishable") = false)

      // --- relax (default non-mutating; inplace=True mutates self) ---
      .def(
          "relax",
          [](Matter &self, bool inplace, bool quiet, bool write_movie,
             bool checkpoint, const std::string &prefix_movie,
             const std::string &prefix_checkpoint,
             bool retain_frames) -> nb::object {
            if (inplace) {
              bool ok = false;
              {
                nb::gil_scoped_release release;
                ok = self.relax(quiet, write_movie, checkpoint, prefix_movie,
                                prefix_checkpoint, retain_frames);
              }
              return nb::make_tuple(nb::cast(self, nb::rv_policy::reference),
                                    ok);
            }
            Matter out(self);
            bool ok = false;
            {
              nb::gil_scoped_release release;
              ok = out.relax(quiet, write_movie, checkpoint, prefix_movie,
                             prefix_checkpoint, retain_frames);
            }
            return nb::make_tuple(std::move(out), ok);
          },
          nb::arg("inplace") = false, nb::arg("quiet") = false,
          nb::arg("write_movie") = false, nb::arg("checkpoint") = false,
          nb::arg("prefix_movie") = "", nb::arg("prefix_checkpoint") = "",
          nb::arg("retain_frames") = false,
          "Minimize. Default: copy then relax (self unchanged). "
          "inplace=True mutates self. retain_frames=True keeps iteration "
          "ConFrames in memory (same stamps as write_movie, no disk required). "
          "Returns (Matter, converged: bool).")
      .def(
          "movie_frames",
          [](Matter &self) { return con_frames_to_python(self.movieFrames()); },
          "list[readcon.ConFrame] retained by relax(retain_frames=True).")
      .def(
          "take_movie_frames",
          [](Matter &self) {
            return con_frames_to_python(self.takeMovieFrames());
          },
          "Take ownership of retained movie frames (clears storage).")
      .def(
          "clear_movie_frames", [](Matter &self) { self.clearMovieFrames(); },
          "Drop retained minimization movie frames.")

      // --- I/O (same ConFileIO / readcon path as CLI client) ---
      .def(
          "con2matter",
          [](Matter &self, const std::string &path) {
            return self.con2matter(path);
          },
          nb::arg("path"), "Load .con file into this Matter")
      .def(
          "matter2con",
          [](Matter &self, const std::string &path, bool append) {
            return self.matter2con(path, append, nullptr);
          },
          nb::arg("path"), nb::arg("append") = false, "Write Matter to .con")
      .def(
          "convel2matter",
          [](Matter &self, const std::string &path) {
            return self.convel2matter(path);
          },
          nb::arg("path"))
      .def(
          "matter2convel",
          [](Matter &self, const std::string &path) {
            return self.matter2convel(path);
          },
          nb::arg("path"))
      .def(
          "matter2xyz",
          [](Matter &self, const std::string &path, bool append) {
            return self.matter2xyz(path, append);
          },
          nb::arg("path"), nb::arg("append") = false)
      .def(
          "write_tibble",
          [](Matter &self, const std::string &path) {
            return self.writeTibble(path);
          },
          nb::arg("path"))

      // --- atom index / header (bindings-friendly) ---
      .def("get_atom_index", &Matter::getAtomIndex, nb::arg("atom"))
      .def("set_atom_index", &Matter::setAtomIndex, nb::arg("atom"),
           nb::arg("index"))
      .def_prop_ro("header_con",
                   [](const Matter &self) { return self.getHeaderCon(); })
      .def(
          "set_header_con_line",
          [](Matter &self, size_t i, const std::string &line) {
            self.setHeaderConLine(i, line);
          },
          nb::arg("i"), nb::arg("line"))

      .def("__len__", &Matter::numberOfAtoms)
      .def("__repr__", [](const Matter &self) {
        // Never trigger a force call: repr runs on every REPL echo, and a
        // stale cache would evaluate the potential with the GIL held.
        // Read potential_energy when the number is wanted.
        const std::string energy =
            self.needsForceUpdate() ? std::string("stale")
                                    : std::to_string(self.getPotentialEnergy());
        return "<Matter n=" + std::to_string(self.numberOfAtoms()) +
               " E=" + energy + ">";
      });
}

} // namespace eonc::pybind
