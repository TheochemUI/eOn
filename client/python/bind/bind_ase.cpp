/*
** ASE ↔ Matter in nanobind: bulk buffer paths (no Python-list geometry builds).
*/
#include "bind_helpers.hpp"
#include "eigen_numpy.hpp"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;

namespace {

// Looser dtype than c_contig-required NpF64 — ASE arrays are usually C-contig
// but nanobind's c_contig cast is strict and throws std::bad_cast.
using ArrF64 = nb::ndarray<nb::numpy, double, nb::device::cpu>;
using ArrI32 = nb::ndarray<nb::numpy, int32_t, nb::device::cpu>;
using ArrI64 = nb::ndarray<nb::numpy, int64_t, nb::device::cpu>;

nb::object np_contig(nb::object obj, const char *dtype) {
  auto np = nb::module_::import_("numpy");
  return np.attr("ascontiguousarray")(obj, nb::arg("dtype") = dtype);
}

const double *f64_ptr(const ArrF64 &a) { return a.data(); }

} // namespace

void bind_ase(nb::module_ &m) {
  using eonc::Matter;
  using eonc::Parameters;
  using eonc::Potential;

  m.def(
      "matter_from_ase",
      [](nb::object atoms, std::shared_ptr<Potential> pot,
         Parameters &params) -> std::shared_ptr<Matter> {
        if (!pot) {
          throw std::invalid_argument(
              "matter_from_ase: potential is required (use "
              "make_potential_from_ase(atoms.calc) or pass pot)");
        }
        auto np = nb::module_::import_("numpy");

        nb::object pos_o = np_contig(atoms.attr("get_positions")(), "float64");
        nb::object cell_o =
            np_contig(np.attr("asarray")(atoms.attr("get_cell")()), "float64");
        nb::object mass_o = np_contig(atoms.attr("get_masses")(), "float64");
        nb::object z_o = np_contig(atoms.attr("get_atomic_numbers")(), "int32");

        ArrF64 pos = nb::cast<ArrF64>(pos_o);
        ArrF64 cell = nb::cast<ArrF64>(cell_o);
        ArrF64 mass = nb::cast<ArrF64>(mass_o);
        ArrI32 z = nb::cast<ArrI32>(z_o);

        if (pos.ndim() != 2 || pos.shape(1) != 3) {
          throw std::invalid_argument(
              "matter_from_ase: positions must be (n,3)");
        }
        const long n = static_cast<long>(pos.shape(0));
        if (mass.ndim() != 1 || static_cast<long>(mass.shape(0)) != n) {
          throw std::invalid_argument(
              "matter_from_ase: masses length must match n_atoms");
        }
        if (z.ndim() != 1 || static_cast<long>(z.shape(0)) != n) {
          throw std::invalid_argument(
              "matter_from_ase: atomic numbers length must match n_atoms");
        }
        if (cell.ndim() != 2 || cell.shape(0) != 3 || cell.shape(1) != 3) {
          throw std::invalid_argument("matter_from_ase: cell must be (3,3)");
        }

        auto matter = std::make_shared<Matter>(pot, params);
        matter->resize(n);
        matter_set_cell_buf(*matter, f64_ptr(cell));
        matter_set_positions_buf(*matter, f64_ptr(pos), n);
        matter_set_masses_buf(*matter, f64_ptr(mass), n);
        matter_set_z_buf(*matter, reinterpret_cast<const int *>(z.data()), n);

        // FixAtoms / FixCartesian → per-axis fixed mask
        try {
          nb::object constraints = atoms.attr("constraints");
          const size_t nc = nb::len(constraints);
          for (size_t ci = 0; ci < nc; ++ci) {
            nb::object c = constraints[nb::int_(ci)];
            if (nb::hasattr(c, "mask") && nb::hasattr(c, "a")) {
              const long i = nb::cast<long>(c.attr("a"));
              if (i < 0 || i >= n) {
                continue;
              }
              nb::object mask_o = np_contig(c.attr("mask"), "int64");
              ArrI64 mask = nb::cast<ArrI64>(mask_o);
              if (mask.ndim() != 1 || mask.shape(0) != 3) {
                continue;
              }
              auto cur = matter->getFixedMask(i);
              for (size_t ax = 0; ax < 3; ++ax) {
                if (mask.data()[ax] != 0) {
                  cur[ax] = true;
                }
              }
              matter->setFixedMask(i, cur);
              continue;
            }
            if (!nb::hasattr(c, "get_indices")) {
              continue;
            }
            nb::object idx_o = np_contig(c.attr("get_indices")(), "int64");
            ArrI64 idx = nb::cast<ArrI64>(idx_o);
            if (idx.ndim() != 1) {
              continue;
            }
            for (size_t k = 0; k < idx.shape(0); ++k) {
              long i = static_cast<long>(idx.data()[k]);
              if (i >= 0 && i < n) {
                matter->setFixed(i, 1);
              }
            }
          }
        } catch (const nb::python_error &) {
          // ignore constraint extraction failures
        } catch (const nb::cast_error &) {
        }

        bool periodic = false;
        try {
          nb::object any =
              np.attr("any")(np.attr("asarray")(atoms.attr("pbc")));
          periodic = nb::cast<bool>(any);
        } catch (...) {
          periodic = true;
        }
        matter->setPeriodic(periodic);
        return matter;
      },
      nb::arg("atoms"), nb::arg("potential"), nb::arg("parameters"),
      nb::keep_alive<0, 2>(), nb::keep_alive<0, 3>(),
      "Build Matter from ASE Atoms via bulk ndarray buffers (C++).");

  m.def(
      "matter_to_ase",
      [](Matter &matter, nb::object pbc_opt) -> nb::object {
        auto ase = nb::module_::import_("ase");
        auto Atoms = ase.attr("Atoms");
        auto np = nb::module_::import_("numpy");

        const long n = matter.numberOfAtoms();
        nb::object positions = np.attr("array")(
            nb::cast(view_n3(matter_positions_ptr(matter), n),
                     nb::rv_policy::reference),
            nb::arg("dtype") = "float64", nb::arg("copy") = true);
        Matrix3d cell_m = matter.getCell();
        nb::object cell = np.attr("array")(
            nb::cast(view_33(cell_m.data()), nb::rv_policy::reference),
            nb::arg("dtype") = "float64", nb::arg("copy") = true);
        nb::object masses = np.attr("asarray")(
            nb::cast(vector_to_numpy(matter.getMasses()), nb::rv_policy::move),
            nb::arg("dtype") = "float64");
        nb::object numbers = np.attr("ascontiguousarray")(
            nb::cast(vectori_to_numpy(matter.getAtomicNrs()),
                     nb::rv_policy::move),
            nb::arg("dtype") = "int64");

        bool use_pbc = matter.getPeriodic();
        if (!pbc_opt.is_none()) {
          use_pbc = nb::cast<bool>(pbc_opt);
        }

        nb::object atoms =
            Atoms(nb::arg("numbers") = numbers,
                  nb::arg("positions") = positions, nb::arg("cell") = cell,
                  nb::arg("pbc") = use_pbc, nb::arg("masses") = masses);

        auto constraints_mod = ase.attr("constraints");
        auto FixAtoms = constraints_mod.attr("FixAtoms");
        auto FixCartesian = constraints_mod.attr("FixCartesian");
        nb::list cons;
        nb::list fully_fixed;
        for (long i = 0; i < n; ++i) {
          const auto mask = matter.getFixedMask(i);
          const bool fx = mask[0];
          const bool fy = mask[1];
          const bool fz = mask[2];
          if (fx && fy && fz) {
            fully_fixed.append(i);
          } else if (fx || fy || fz) {
            cons.append(FixCartesian(
                i, nb::make_tuple(fx ? 1 : 0, fy ? 1 : 0, fz ? 1 : 0)));
          }
        }
        if (nb::len(fully_fixed) > 0) {
          cons.append(FixAtoms(nb::arg("indices") = fully_fixed));
        }
        if (nb::len(cons) > 0) {
          atoms.attr("set_constraint")(cons);
        }
        return atoms;
      },
      nb::arg("matter"), nb::arg("pbc") = nb::none(),
      "Matter → ase.Atoms (geometry from Matter buffers).");
}

} // namespace eonc::pybind
