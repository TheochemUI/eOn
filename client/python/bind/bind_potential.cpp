#include "bind_helpers.hpp"
#include "eon/BaseStructures.h"
#include "eon/Eigen.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include <algorithm>
#include <cstring>
#include <nanobind/ndarray.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

#include <magic_enum/magic_enum.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;

namespace {

/// ASE Calculator as eOn Potential — designed to drive Matter, not Atoms.
///
/// Model:
///   Matter owns geometry (Eigen). ASE calc is the force engine.
///   MD/relax/NEB call Matter → Potential::force(R*, …) with R* =
///   Matter.positionsData(). This class keeps one peer ase.Atoms and:
///
///   1. **bind_matter(m)**: wire ASE ``arrays['positions']`` / cell to
///      Matter buffers when possible (shared storage; no per-step copy).
///   2. **system_changes** tracked only in this class (pointer identity +
///      cell memcmp) — no Matter API for ASE cache state.
///   3. Fallback path when force() is called with foreign buffers (e.g.
///      Potential.get_ef from Python arrays): bulk view + set_positions.
///
/// Chemist API (Python)::
///
///   pot = pyec.ase_potential(EMT())
///   m = pyec.Matter(pot, params);  m.resize(...); fill geometry
///   pot.bind_matter(m)            # optional; from_ase does this
///   m.potential_energy / m.relax()  # ASE under the hood
class AseCalcPotential final : public eonc::Potential {
  nb::object calc_;
  nb::object atoms_;
  nb::object np_;
  nb::object props_;
  nb::object all_changes_;
  nb::object changes_positions_; // ['positions']
  nb::object changes_cell_;      // ['cell']
  nb::object changes_pos_cell_;  // ['positions', 'cell']
  nb::object changes_empty_;     // []
  std::vector<int> z_cache_;
  long n_cached_{-1};
  bool modules_ready_{false};

  // Bound Matter (not owned). Cache state is binding-local only.
  eonc::Matter *linked_{nullptr};
  bool shared_positions_{false};
  bool cell_synced_{false};
  bool ever_calculated_{false};
  const double *last_R_{nullptr};
  const double *last_box_{nullptr};
  double last_cell_[9]{};

  void ensure_modules() {
    if (modules_ready_)
      return;
    np_ = nb::module_::import_("numpy");
    auto calc_mod = nb::module_::import_("ase.calculators.calculator");
    all_changes_ = calc_mod.attr("all_changes");
    props_ = nb::list();
    props_.attr("append")("energy");
    props_.attr("append")("forces");
    changes_positions_ = nb::list();
    changes_positions_.attr("append")("positions");
    changes_cell_ = nb::list();
    changes_cell_.attr("append")("cell");
    changes_pos_cell_ = nb::list();
    changes_pos_cell_.attr("append")("positions");
    changes_pos_cell_.attr("append")("cell");
    changes_empty_ = nb::list();
    modules_ready_ = true;
  }

  static nb::object view_f64(const double *data, size_t d0, size_t d1) {
    nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu> arr(
        const_cast<double *>(data), {d0, d1});
    return nb::cast(arr, nb::rv_policy::reference);
  }

  void ensure_atoms(long nAtoms, const int *atomicNrs, bool pbc) {
    const size_t n = static_cast<size_t>(nAtoms);
    const bool size_changed = (n_cached_ != nAtoms);
    bool z_changed = size_changed;
    if (!z_changed) {
      z_changed = !std::equal(z_cache_.begin(), z_cache_.end(), atomicNrs);
    }
    if (!size_changed && !z_changed && atoms_.is_valid() && !atoms_.is_none())
      return;

    z_cache_.assign(atomicNrs, atomicNrs + n);
    n_cached_ = nAtoms;
    shared_positions_ = false;
    cell_synced_ = false;
    ever_calculated_ = false;

    nb::ndarray<nb::numpy, int32_t, nb::c_contig, nb::device::cpu> z_arr(
        z_cache_.data(), {n});
    nb::object numbers = nb::cast(z_arr, nb::rv_policy::reference);
    nb::object zeros = np_.attr("zeros");
    nb::object positions =
        zeros(nb::make_tuple(nAtoms, 3), nb::arg("dtype") = "float64");
    nb::object cell = np_.attr("eye")(3, nb::arg("dtype") = "float64");

    nb::object Atoms = nb::module_::import_("ase").attr("Atoms");
    atoms_ =
        Atoms(nb::arg("numbers") = numbers, nb::arg("positions") = positions,
              nb::arg("cell") = cell, nb::arg("pbc") = pbc);
    atoms_.attr("calc") = calc_;
  }

  /// Point ASE positions array at Matter storage (zero-copy when ASE allows).
  bool try_share_positions(eonc::Matter &m) {
    try {
      nb::object pos_view =
          view_f64(matter_positions_ptr(m),
                   static_cast<size_t>(m.numberOfAtoms()), size_t{3});
      // Direct arrays dict install — np.asarray on a view does not copy.
      atoms_.attr("arrays").attr("__setitem__")("positions", pos_view);
      // Sanity: must still look like (n,3)
      nb::object p = atoms_.attr("positions");
      auto sh = nb::cast<nb::tuple>(p.attr("shape"));
      if (nb::len(sh) != 2 || nb::cast<long>(sh[0]) != m.numberOfAtoms()) {
        shared_positions_ = false;
        return false;
      }
      shared_positions_ = true;
      return true;
    } catch (...) {
      shared_positions_ = false;
      return false;
    }
  }

  static void copy_forces_out(nb::object forces_obj, long nAtoms, double *F) {
    const size_t n3 = static_cast<size_t>(nAtoms) * 3u;
    try {
      auto farr =
          nb::cast<nb::ndarray<nb::numpy, double, nb::device::cpu>>(forces_obj);
      if (farr.ndim() == 2 && static_cast<long>(farr.shape(0)) == nAtoms &&
          farr.shape(1) == 3) {
        if (farr.stride(0) == 3 && farr.stride(1) == 1) {
          std::memcpy(F, farr.data(), n3 * sizeof(double));
          return;
        }
        const double *base = farr.data();
        const auto s0 = static_cast<size_t>(farr.stride(0));
        const auto s1 = static_cast<size_t>(farr.stride(1));
        for (long i = 0; i < nAtoms; ++i)
          for (long j = 0; j < 3; ++j)
            F[i * 3 + j] =
                base[static_cast<size_t>(i) * s0 + static_cast<size_t>(j) * s1];
        return;
      }
    } catch (const nb::cast_error &) {
    }
    auto np = nb::module_::import_("numpy");
    nb::object cont =
        np.attr("ascontiguousarray")(forces_obj, nb::arg("dtype") = "float64");
    auto farr =
        nb::cast<nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>>(
            cont);
    if (farr.ndim() != 2 || static_cast<long>(farr.shape(0)) != nAtoms ||
        farr.shape(1) != 3) {
      throw std::runtime_error("ASE forces must be shape (n_atoms, 3) float64");
    }
    std::memcpy(F, farr.data(), n3 * sizeof(double));
  }

  nb::object pick_system_changes(bool pos_changed, bool cell_changed,
                                 bool first) {
    if (first)
      return all_changes_;
    if (pos_changed && cell_changed)
      return changes_pos_cell_;
    if (pos_changed)
      return changes_positions_;
    if (cell_changed)
      return changes_cell_;
    return changes_empty_;
  }

public:
  explicit AseCalcPotential(nb::object calc)
      : eonc::Potential(eonc::PotType::ASE_POT),
        calc_(std::move(calc)) {
    if (!calc_.is_valid() || calc_.is_none()) {
      throw std::invalid_argument(
          "make_potential_from_ase: calculator is None");
    }
  }

  /// Attach this calculator to a Matter for shared geometry + gen cache.
  /// Call after Matter is resized and Z/cell filled (from_ase does this).
  void bind_matter(eonc::Matter &m) {
    nb::gil_scoped_acquire gil;
    ensure_modules();
    linked_ = &m;
    const long n = m.numberOfAtoms();
    if (n <= 0) {
      return;
    }
    // Z from public API (by-value VectorXi)
    VectorXi z = m.getAtomicNrs();
    ensure_atoms(n, z.data(), m.getPeriodic());
    Matrix3d cell = m.getCell();
    nb::object cell_view = view_f64(cell.data(), size_t{3}, size_t{3});
    atoms_.attr("set_cell")(cell_view, nb::arg("scale_atoms") = false);
    std::memcpy(last_cell_, cell.data(), 9 * sizeof(double));
    cell_synced_ = true;
    try_share_positions(m);
    if (!shared_positions_) {
      nb::object pos_view =
          view_f64(matter_positions_ptr(m), static_cast<size_t>(n), size_t{3});
      atoms_.attr("set_positions")(pos_view);
    }
    ever_calculated_ = false;
    last_R_ = matter_positions_ptr(m);
    // cell is by-value on Matter; force path uses box pointer from force()
    last_box_ = nullptr;
  }

  void unbind_matter() {
    linked_ = nullptr;
    shared_positions_ = false;
    ever_calculated_ = false;
  }

  [[nodiscard]] bool is_bound() const noexcept { return linked_ != nullptr; }

  void force(long nAtoms, const double *R, const int *atomicNrs, double *F,
             double *U, double * /*variance*/, const double *box) override {
    if (nAtoms <= 0) {
      *U = 0.0;
      return;
    }
    if (!R || !atomicNrs || !F || !U || !box) {
      throw std::invalid_argument("ase_potential force: null buffer");
    }
    nb::gil_scoped_acquire gil;
    try {
      ensure_modules();
      const bool pbc = linked_ ? linked_->getPeriodic() : true;
      ensure_atoms(nAtoms, atomicNrs, pbc);

      // Linked when force is driven by Matter::computePotential (R is
      // Matter positions storage). Cache is binding-local only.
      const bool from_linked = linked_ && nAtoms == linked_->numberOfAtoms() &&
                               R == matter_positions_ptr(*linked_);

      bool first = !ever_calculated_;
      bool pos_changed = true;
      bool cell_changed =
          (std::memcmp(box, last_cell_, 9 * sizeof(double)) != 0);

      if (from_linked) {
        if (!shared_positions_) {
          try_share_positions(*linked_);
        }
        if (shared_positions_) {
          // Recompute only when Matter says geometry is dirty (public API).
          pos_changed = linked_->needsForceUpdate();
        } else {
          nb::object pos_view =
              view_f64(R, static_cast<size_t>(nAtoms), size_t{3});
          atoms_.attr("set_positions")(pos_view);
          pos_changed = true;
        }
        if (cell_changed || !cell_synced_) {
          nb::object cell_view = view_f64(box, size_t{3}, size_t{3});
          atoms_.attr("set_cell")(cell_view, nb::arg("scale_atoms") = false);
          std::memcpy(last_cell_, box, 9 * sizeof(double));
          cell_synced_ = true;
          cell_changed = true;
        }
      } else {
        nb::object pos_view =
            view_f64(R, static_cast<size_t>(nAtoms), size_t{3});
        nb::object cell_view = view_f64(box, size_t{3}, size_t{3});
        atoms_.attr("set_positions")(pos_view);
        atoms_.attr("set_cell")(cell_view, nb::arg("scale_atoms") = false);
        std::memcpy(last_cell_, box, 9 * sizeof(double));
        shared_positions_ = false;
        first = true;
        pos_changed = true;
        cell_changed = true;
      }

      last_R_ = R;
      last_box_ = box;

      nb::object sys_changes =
          pick_system_changes(pos_changed, cell_changed, first);
      ever_calculated_ = true;

      nb::object forces_obj;
      try {
        calc_.attr("calculate")(atoms_, props_, sys_changes);
        nb::object results = calc_.attr("results");
        *U = nb::cast<double>(results["energy"]);
        forces_obj = results["forces"];
      } catch (const nb::python_error &) {
        *U = nb::cast<double>(atoms_.attr("get_potential_energy")());
        forces_obj = atoms_.attr("get_forces")();
      }
      copy_forces_out(forces_obj, nAtoms, F);
    } catch (const nb::python_error &e) {
      throw std::runtime_error(std::string("ASE calculator Python error: ") +
                               e.what());
    } catch (const std::exception &e) {
      throw std::runtime_error(std::string("ASE calculator error: ") +
                               e.what());
    }
  }

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return false;
  }
};

// Type-erased holder so Python can call bind_matter on ASE pots only.
struct AsePotentialHandle {
  std::shared_ptr<AseCalcPotential> pot;
  void bind_matter(eonc::Matter &m) { pot->bind_matter(m); }
  void unbind_matter() { pot->unbind_matter(); }
  [[nodiscard]] bool is_bound() const { return pot && pot->is_bound(); }
  std::shared_ptr<eonc::Potential> as_potential() const { return pot; }
};

} // namespace

void bind_potential(nb::module_ &m) {
  nb::class_<eonc::Potential>(m, "Potential",
                              "Abstract potential (force + energy)")
      .def_prop_ro("type", &eonc::Potential::getType)
      .def_prop_ro("force_call_counter",
                   [](const eonc::Potential &self) {
                     return self.forceCallCounter.load();
                   })
      .def("is_surrogate", &eonc::Potential::isSurrogate)
      .def("requires_isolated_molecule_layout",
           &eonc::Potential::requiresIsolatedMoleculeLayout)
      .def("is_thread_safe", &eonc::Potential::isThreadSafe)
      .def(
          "get_ef",
          [](eonc::Potential &self,
             nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu> pos,
             nb::ndarray<nb::numpy, int64_t, nb::c_contig, nb::device::cpu> z,
             nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>
                 box) {
            if (pos.ndim() != 2 || pos.shape(1) != 3)
              throw std::invalid_argument("pos must be (n,3)");
            if (z.ndim() != 1 || z.shape(0) != pos.shape(0))
              throw std::invalid_argument("z must be (n,)");
            if (box.ndim() != 2 || box.shape(0) != 3 || box.shape(1) != 3)
              throw std::invalid_argument("box must be (3,3)");
            const auto n = static_cast<Eigen::Index>(pos.shape(0));
            AtomMatrix R =
                AtomMatrix::Map(const_cast<double *>(pos.data()), n, 3);
            Matrix3d cell = Matrix3d::Map(const_cast<double *>(box.data()));
            VectorXi atmnrs(n);
            for (Eigen::Index i = 0; i < n; ++i)
              atmnrs(i) = static_cast<int>(z.data()[i]);
            // Torch autograd in RGPOT metatomic engines refuses the GIL.
            // Release around the evaluation only; argument validation above
            // and result construction below touch Python objects.
            // An ASE-backed potential re-enters Python from force(), where
            // nb::gil_scoped_acquire re-takes the GIL on this thread.
            double energy = 0.0;
            AtomMatrix forces;
            {
              nb::gil_scoped_release release;
              auto ef = self.get_ef(R, atmnrs, cell);
              energy = std::get<0>(ef);
              forces = std::move(std::get<1>(ef));
            }
            const size_t rows = static_cast<size_t>(forces.rows());
            const size_t cols = 3;
            double *buf = new double[rows * cols];
            std::memcpy(buf, forces.data(), rows * cols * sizeof(double));
            nb::capsule owner(buf, [](void *p) noexcept {
              delete[] static_cast<double *>(p);
            });
            auto f_np = nb::ndarray<nb::numpy, double, nb::c_contig>(
                buf, {rows, cols}, owner);
            return nb::make_tuple(energy, f_np);
          },
          nb::arg("positions"), nb::arg("atomic_numbers"), nb::arg("cell"),
          "Evaluate energy and forces: returns (energy, forces (n,3))")
      .def("__repr__", [](const eonc::Potential &self) {
        return "<Potential " +
               std::string(magic_enum::enum_name(self.getType())) + ">";
      });

  nb::class_<AsePotentialHandle>(m, "AsePotential",
                                 "ASE Calculator bound as eOn Potential. "
                                 "Use bind_matter(m) so MD/relax hit shared "
                                 "geometry + system_changes caching.")
      .def(
          "bind_matter", &AsePotentialHandle::bind_matter, nb::arg("matter"),
          nb::keep_alive<1, 2>(),
          "Wire ASE peer Atoms to this Matter (shared positions when possible)")
      .def("unbind_matter", &AsePotentialHandle::unbind_matter)
      .def_prop_ro("is_bound", &AsePotentialHandle::is_bound)
      .def_prop_ro("potential", &AsePotentialHandle::as_potential,
                   "Underlying Potential for Matter(pot, params)");

  m.def(
      "make_potential",
      [](eonc::PotType ptype, eonc::Parameters &params) {
        auto pot = eonc::helpers::makePotential(ptype, params);
        if (!pot) {
          throw std::runtime_error("make_potential returned null");
        }
        return pot;
      },
      nb::arg("pot_type"), nb::arg("parameters"),
      "Construct a Potential from PotType + Parameters");

  m.def(
      "make_potential",
      [](const std::string &name, eonc::Parameters &params) {
        auto v = magic_enum::enum_cast<eonc::PotType>(
            name, magic_enum::case_insensitive);
        if (!v) {
          throw std::invalid_argument("unknown PotType: " + name);
        }
        params.potential_options.potential = *v;
        auto pot = eonc::helpers::makePotential(*v, params);
        if (!pot) {
          throw std::runtime_error("make_potential returned null for " + name);
        }
        return pot;
      },
      nb::arg("name"), nb::arg("parameters"),
      "Construct a Potential from config-style name + Parameters");

  m.def(
      "make_potential",
      [](eonc::Parameters &params) {
        auto pot = eonc::helpers::makePotential(params);
        if (!pot) {
          throw std::runtime_error("make_potential returned null");
        }
        return pot;
      },
      nb::arg("parameters"), "Construct a Potential from Parameters.potential");

  // ASE Calculator → Potential (runtime ASE; no -Dwith_ase compile flag).
  m.def(
      "make_potential_from_ase",
      [](nb::object calculator) -> std::shared_ptr<eonc::Potential> {
        return std::make_shared<AseCalcPotential>(std::move(calculator));
      },
      nb::arg("calculator"),
      "Wrap ASE Calculator as Potential. Prefer ase_potential() + "
      "bind_matter.");

  m.def(
      "ase_potential",
      [](nb::object calculator) -> AsePotentialHandle {
        AsePotentialHandle h;
        h.pot = std::make_shared<AseCalcPotential>(std::move(calculator));
        return h;
      },
      nb::arg("calculator"),
      "ASE Calculator as first-class pyeonclient potential. "
      "m = Matter(handle.potential, params); handle.bind_matter(m).");

  m.def(
      "attach_ase_calculator",
      [](eonc::Matter &matter, nb::object calculator) {
        auto pot = std::make_shared<AseCalcPotential>(std::move(calculator));
        pot->bind_matter(matter);
        matter.setPotential(pot);
      },
      nb::arg("matter"), nb::arg("calculator"),
      "Set Matter's potential to ASE calc and bind shared geometry.");

  m.def(
      "bind_ase_matter",
      [](std::shared_ptr<eonc::Potential> pot, eonc::Matter &matter) -> bool {
        auto *ase = dynamic_cast<AseCalcPotential *>(pot.get());
        if (!ase) {
          return false;
        }
        ase->bind_matter(matter);
        return true;
      },
      nb::arg("potential"), nb::arg("matter"),
      "If potential is an ASE wrap, bind it to Matter for shared geometry. "
      "Returns True if bound.");
}

} // namespace eonc::pybind
