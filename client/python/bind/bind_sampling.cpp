/*
** Class-first Matter sampling / long-timescale surface (import as pyec).
** Free run_* removed from public API; construct then run(inplace=False).
*/
#include "bind_helpers.hpp"
#include "eigen_numpy.hpp"
#include "eon/Dynamics.h"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/MonteCarlo.h"
#include "eon/ParallelReplicaJob.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/ReplicaExchangeJob.h"
#include "eon/SafeHyperJob.h"
#include "eon/TADJob.h"

#ifdef WITH_GP_SURROGATE
#include "eon/GPSurrogateJob.h"
#include "eon/NudgedElasticBand.h"
#endif

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;

namespace {

void ensure_dynamics_steps(eonc::Parameters &params, long default_steps) {
  if (params.dynamics_options.steps <= 0)
    params.dynamics_options.steps = default_steps;
  if (params.dynamics_options.time_step <= 0.0) {
    params.dynamics_options.time_step =
        params.dynamics_options.time_step_input / params.constants.timeUnit;
    if (params.dynamics_options.time_step <= 0.0)
      params.dynamics_options.time_step = 1.0 / params.constants.timeUnit;
  }
}

void ensure_long_timescale_params(eonc::Parameters &params,
                                  long default_steps) {
  ensure_dynamics_steps(params, default_steps);
  auto &dyn = params.dynamics_options;
  const double dt = dyn.time_step;
  const double horizon = dt * static_cast<double>(dyn.steps);

  auto &pr = params.parallel_replica_options;
  if (pr.state_check_interval <= 0.0)
    pr.state_check_interval = std::max(dt, horizon);
  if (pr.record_interval <= 0.0)
    pr.record_interval = std::max(dt, pr.state_check_interval / 2.0);
  if (pr.record_interval > pr.state_check_interval)
    pr.record_interval = pr.state_check_interval;
  if (pr.state_check_interval < dt)
    pr.state_check_interval = dt;
  if (pr.record_interval < dt)
    pr.record_interval = dt;
  if (pr.dephase_time <= 0.0)
    pr.dephase_time = dt;
  if (pr.corr_time <= 0.0)
    pr.corr_time = dt;
  if (dyn.steps > 0 && dyn.steps < 50) {
    pr.refine_transition = false;
    pr.dephase_loop_stop = true;
    if (pr.dephase_loop_max <= 0 || pr.dephase_loop_max > 5)
      pr.dephase_loop_max = 2;
  }

  auto &rex = params.replica_exchange_options;
  if (rex.replicas < 2)
    rex.replicas = 2;
  if (rex.sampling_time <= 0.0 ||
      (dyn.steps > 0 && dyn.steps < 50 && rex.sampling_time > horizon * 2.0))
    rex.sampling_time = horizon > 0.0 ? horizon : dt;
  if (rex.exchange_period <= 0.0 || rex.exchange_period > rex.sampling_time)
    rex.exchange_period = std::max(dt, rex.sampling_time / 2.0);
  if (rex.temperature_low <= 0.0)
    rex.temperature_low = params.main_options.temperature > 0.0
                              ? params.main_options.temperature
                              : 300.0;
  if (rex.temperature_high <= rex.temperature_low)
    rex.temperature_high = rex.temperature_low * 1.5;
  if (rex.exchange_trials <= 0)
    rex.exchange_trials = std::max(1L, rex.replicas - 1);
  if (dyn.steps > 0 && dyn.steps < 50 && rex.replicas > 3)
    rex.replicas = 2;
}

struct SeededAlgo {
  std::shared_ptr<eonc::Matter> seed;
  eonc::Parameters params;
  std::shared_ptr<eonc::Potential> pot;
  std::shared_ptr<eonc::Matter> result;

  SeededAlgo(std::shared_ptr<eonc::Matter> m, eonc::Parameters p,
             std::shared_ptr<eonc::Potential> pot_in)
      : seed(std::move(m)),
        params(std::move(p)),
        pot(std::move(pot_in)) {
    if (!seed)
      throw std::runtime_error("matter is required");
    if (!pot)
      pot = seed->getPotential();
  }
};

/// Hand a Matter back to Python holding the seed's wrapper, which owns the
/// Parameters the two share (Matter keeps only a non-owning pointer).
///
/// The patient is the seed, not the algorithm object. The algorithm already
/// holds the seed through nb::keep_alive<1, 2> on its constructor, so an
/// algorithm-as-patient edge closes an uncollectable cycle whenever the result
/// *is* the seed: every inplace=True run, and every `matter` read before the
/// first run. Tying to the seed makes that case a self-tie, which
/// tie_lifetime skips, and leaves no edge back from the seed to the result.
template <typename Algo>
nb::object algo_matter_out(const Algo &self,
                           std::shared_ptr<eonc::Matter> out) {
  nb::object obj = nb::cast(out);
  tie_lifetime(obj, matter_object(*self.seed));
  return obj;
}

} // namespace

void bind_sampling(nb::module_ &m) {
  using eonc::Dynamics;
  using eonc::DynamicsConfig;
  using eonc::Matter;
  using eonc::MinModeSaddleSearch;
  using eonc::MonteCarlo;
  using eonc::Parameters;
  using eonc::Potential;

  m.def(
      "built_with_gp_surrogate",
      []() {
#ifdef WITH_GP_SURROGATE
        return true;
#else
        return false;
#endif
      },
      "True if compiled with WITH_GP_SURROGATE / meson with_gp_surrogate");

  // --- MolecularDynamics ---
  struct PyMD : SeededAlgo {
    using SeededAlgo::SeededAlgo;
    std::shared_ptr<Matter> run(bool inplace) {
      ensure_dynamics_steps(params, 10);
      auto work = matter_work(seed, inplace);
      if (pot)
        work->setPotential(pot);
      Dynamics dyn(work.get(), DynamicsConfig::fromParams(params));
      {
        nb::gil_scoped_release release;
        dyn.run();
      }
      result = matter_result(seed, work, inplace);
      return result;
    }
  };

  nb::class_<PyMD>(m, "MolecularDynamics",
                   "Velocity-Verlet MD on Matter (C++ Dynamics). "
                   "run(inplace=False) default non-mutating.")
      .def(nb::init<std::shared_ptr<Matter>, Parameters,
                    std::shared_ptr<Potential>>(),
           nb::arg("matter"), nb::arg("parameters"),
           nb::arg("potential") = nb::none(), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 4>())
      .def(
          "run",
          [](PyMD &self, bool inplace) {
            return algo_matter_out(self, self.run(inplace));
          },
          nb::arg("inplace") = false,
          "Run MD; returns Matter. inplace=True mutates seed.")
      .def_prop_ro("matter", [](const PyMD &s) {
        return algo_matter_out(s, s.result ? s.result : s.seed);
      });

  // --- MonteCarlo ---
  struct PyMC : SeededAlgo {
    using SeededAlgo::SeededAlgo;
    std::shared_ptr<Matter> run(bool inplace) {
      auto work = matter_work(seed, inplace);
      MonteCarlo mc(work, params);
      int steps = params.monte_carlo_options.steps;
      if (steps <= 0)
        steps = 10;
      {
        nb::gil_scoped_release release;
        mc.run(steps, params.main_options.temperature,
               params.monte_carlo_options.step_size);
      }
      result = matter_result(seed, work, inplace);
      return result;
    }
  };

  nb::class_<PyMC>(m, "MonteCarlo",
                   "Metropolis MC on Matter. run(inplace=False).")
      .def(nb::init<std::shared_ptr<Matter>, Parameters,
                    std::shared_ptr<Potential>>(),
           nb::arg("matter"), nb::arg("parameters"),
           nb::arg("potential") = nb::none(), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 4>())
      .def(
          "run",
          [](PyMC &self, bool inplace) {
            return algo_matter_out(self, self.run(inplace));
          },
          nb::arg("inplace") = false)
      .def_prop_ro("matter", [](const PyMC &s) {
        return algo_matter_out(s, s.result ? s.result : s.seed);
      });

  // --- BasinHopping ---
  struct PyBH : SeededAlgo {
    using SeededAlgo::SeededAlgo;
    std::shared_ptr<Matter> run(bool inplace) {
      if (!pot)
        throw std::runtime_error("BasinHopping: potential required");
      auto work = matter_work(seed, inplace);
      work->setPotential(pot);
      long nsteps = params.basin_hopping_options.steps;
      if (nsteps <= 0)
        nsteps = 5;
      double max_disp = params.basin_hopping_options.displacement;
      double kT = params.constants.kB * params.main_options.temperature;
      std::mt19937_64 rng(
          params.main_options.randomSeed >= 0
              ? static_cast<uint64_t>(params.main_options.randomSeed)
              : 0xC0FFEEULL);
      std::uniform_real_distribution<double> uni(0.0, 1.0);
      std::normal_distribution<double> gauss(0.0, 1.0);
      {
        nb::gil_scoped_release release;
        work->relax(true);
        double Ecur = work->getPotentialEnergy();
        Matter best(*work);
        double Ebest = Ecur;
        for (long step = 0; step < nsteps; ++step) {
          Matter trial(*work);
          AtomMatrix pos = trial.getPositions();
          for (long i = 0; i < trial.numberOfAtoms(); ++i) {
            if (trial.getFixed(i))
              continue;
            for (int k = 0; k < 3; ++k)
              pos(i, k) += max_disp * gauss(rng);
          }
          trial.setPositions(pos);
          trial.relax(true);
          double Etrial = trial.getPotentialEnergy();
          double dE = Etrial - Ecur;
          if (dE <= 0.0 || uni(rng) < std::exp(-dE / kT)) {
            *work = trial;
            Ecur = Etrial;
            if (Ecur < Ebest) {
              best = *work;
              Ebest = Ecur;
            }
          }
        }
        *work = best;
      }
      result = matter_result(seed, work, inplace);
      return result;
    }
  };

  nb::class_<PyBH>(m, "BasinHopping",
                   "Basin-hopping Metropolis + displace + relax. "
                   "run(inplace=False).")
      .def(nb::init<std::shared_ptr<Matter>, Parameters,
                    std::shared_ptr<Potential>>(),
           nb::arg("matter"), nb::arg("parameters"), nb::arg("potential"),
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 4>())
      .def(
          "run",
          [](PyBH &self, bool inplace) {
            return algo_matter_out(self, self.run(inplace));
          },
          nb::arg("inplace") = false)
      .def_prop_ro("matter", [](const PyBH &s) {
        return algo_matter_out(s, s.result ? s.result : s.seed);
      });

  // --- ProcessSearch ---
  struct PyProcessSearch {
    std::shared_ptr<Matter> seed;
    AtomMatrix mode;
    Parameters params;
    std::shared_ptr<Potential> pot;
    std::shared_ptr<Matter> reactant_out;
    std::shared_ptr<Matter> saddle_out;
    int status{0};

    PyProcessSearch(std::shared_ptr<Matter> m, const NpF64 &mode_np,
                    Parameters p, std::shared_ptr<Potential> pot_in)
        : seed(std::move(m)),
          mode(atom_matrix_from_numpy(mode_np)),
          params(std::move(p)),
          pot(std::move(pot_in)) {
      if (!seed)
        throw std::runtime_error("ProcessSearch: matter required");
      if (!pot)
        pot = seed->getPotential();
    }

    nb::object run(bool inplace) {
      auto reactant = matter_work(seed, inplace);
      reactant->setPotential(pot);
      if (params.process_search_options.minimize_first) {
        nb::gil_scoped_release release;
        reactant->relax(true);
      }
      double E0 = 0.0;
      {
        nb::gil_scoped_release release;
        E0 = reactant->getPotentialEnergy();
      }
      auto saddle = std::make_shared<Matter>(*reactant);
      double mag = params.saddle_search_options.displace_magnitude;
      if (mag > 0.0) {
        AtomMatrix pos = saddle->getPositions();
        pos += mag * mode;
        saddle->setPositions(pos);
      }
      MinModeSaddleSearch ss(saddle, mode, E0, params, pot);
      {
        nb::gil_scoped_release release;
        status = ss.run();
      }
      reactant_out = reactant;
      saddle_out = saddle;
      // Both share the seed's Parameters pointer, and keep_alive cannot reach
      // into the tuple.
      nb::object seed_obj = matter_object(*seed);
      nb::object reactant_obj = nb::cast(reactant);
      nb::object saddle_obj = nb::cast(saddle);
      tie_lifetime(reactant_obj, seed_obj);
      tie_lifetime(saddle_obj, seed_obj);
      return nb::make_tuple(reactant_obj, saddle_obj, status);
    }
  };

  nb::class_<PyProcessSearch>(m, "ProcessSearch",
                              "Min-mode process search. run(inplace=False) -> "
                              "(reactant, saddle, status).")
      .def(nb::init<std::shared_ptr<Matter>, const NpF64 &, Parameters,
                    std::shared_ptr<Potential>>(),
           nb::arg("matter"), nb::arg("mode"), nb::arg("parameters"),
           nb::arg("potential"), nb::keep_alive<1, 2>(), nb::keep_alive<1, 5>())
      .def("run", &PyProcessSearch::run, nb::arg("inplace") = false)
      // Typed so `search.status == SaddleStatus.GOOD` compares against the
      // same type; run() keeps returning the raw int code.
      .def_prop_ro("status", [](const PyProcessSearch &s) {
        return static_cast<MinModeSaddleSearch::Status>(s.status);
      });

  // --- Structure helpers (predicates; neither argument is modified) ---
  m.def(
      "structures_equal",
      [](const Matter &a, const Matter &b, bool indistinguishable) {
        // Matter::compare is non-const: with remove_translation set it runs
        // translationRemove, which calls setPositions on its left operand,
        // translating and PBC-wrapping it and invalidating the force cache.
        // A predicate must not do that to the caller, so it runs on a copy.
        Matter probe(a);
        return probe.compare(b, indistinguishable);
      },
      nb::arg("a"), nb::arg("b"), nb::arg("indistinguishable") = false,
      "True when the two structures match under Parameters.structure_"
      "comparison. Neither argument is modified.");

  m.def(
      "structure_distance",
      [](Matter &a, Matter &b) { return a.distanceTo(b); }, nb::arg("a"),
      nb::arg("b"));

  // --- TAD ---
  struct PyTAD : SeededAlgo {
    using SeededAlgo::SeededAlgo;
    std::shared_ptr<Matter> run(bool inplace) {
      ensure_long_timescale_params(params, 5);
      auto work = matter_work(seed, inplace);
      if (pot)
        work->setPotential(pot);
      auto job_pot = pot ? pot : work->getPotential();
      eonc::TADJob job(job_pot, params);
      std::shared_ptr<Matter> out;
      {
        nb::gil_scoped_release release;
        out = job.runFromMatter(work);
      }
      result = matter_result(seed, out, inplace);
      return result;
    }
  };

  nb::class_<PyTAD>(m, "TAD",
                    "Temperature Accelerated Dynamics (TADJob). "
                    "run(inplace=False).")
      .def(nb::init<std::shared_ptr<Matter>, Parameters,
                    std::shared_ptr<Potential>>(),
           nb::arg("matter"), nb::arg("parameters"),
           nb::arg("potential") = nb::none(), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 4>())
      .def(
          "run",
          [](PyTAD &self, bool inplace) {
            return algo_matter_out(self, self.run(inplace));
          },
          nb::arg("inplace") = false)
      .def_prop_ro("matter", [](const PyTAD &s) {
        return algo_matter_out(s, s.result ? s.result : s.seed);
      });

  // Mechanical shells for PR family (not product focus)
  auto bind_pr_like = [&](const char *name, auto make_and_run) {
    // use explicit classes below instead
    (void)name;
    (void)make_and_run;
  };
  (void)bind_pr_like;

  struct PyPR : SeededAlgo {
    using SeededAlgo::SeededAlgo;
    std::shared_ptr<Matter> run(bool inplace) {
      ensure_long_timescale_params(params, 5);
      auto work = matter_work(seed, inplace);
      if (pot)
        work->setPotential(pot);
      auto job_pot = pot ? pot : work->getPotential();
      eonc::ParallelReplicaJob job(job_pot, params);
      std::shared_ptr<Matter> out;
      {
        nb::gil_scoped_release release;
        out = job.runFromMatter(work);
      }
      result = matter_result(seed, out, inplace);
      return result;
    }
  };
  nb::class_<PyPR>(m, "ParallelReplica",
                   "Parallel Replica (mechanical class; run(inplace=False)).")
      .def(nb::init<std::shared_ptr<Matter>, Parameters,
                    std::shared_ptr<Potential>>(),
           nb::arg("matter"), nb::arg("parameters"),
           nb::arg("potential") = nb::none(), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 4>())
      .def(
          "run",
          [](PyPR &self, bool inplace) {
            return algo_matter_out(self, self.run(inplace));
          },
          nb::arg("inplace") = false);

  struct PySH : SeededAlgo {
    using SeededAlgo::SeededAlgo;
    std::shared_ptr<Matter> run(bool inplace) {
      ensure_long_timescale_params(params, 5);
      auto work = matter_work(seed, inplace);
      if (pot)
        work->setPotential(pot);
      auto job_pot = pot ? pot : work->getPotential();
      eonc::SafeHyperJob job(job_pot, params);
      std::shared_ptr<Matter> out;
      {
        nb::gil_scoped_release release;
        out = job.runFromMatter(work);
      }
      result = matter_result(seed, out, inplace);
      return result;
    }
  };
  nb::class_<PySH>(m, "SafeHyperdynamics",
                   "Safe Hyperdynamics. run(inplace=False).")
      .def(nb::init<std::shared_ptr<Matter>, Parameters,
                    std::shared_ptr<Potential>>(),
           nb::arg("matter"), nb::arg("parameters"),
           nb::arg("potential") = nb::none(), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 4>())
      .def(
          "run",
          [](PySH &self, bool inplace) {
            return algo_matter_out(self, self.run(inplace));
          },
          nb::arg("inplace") = false);

  struct PyREX : SeededAlgo {
    using SeededAlgo::SeededAlgo;
    std::shared_ptr<Matter> run(bool inplace) {
      ensure_long_timescale_params(params, 5);
      auto work = matter_work(seed, inplace);
      if (pot)
        work->setPotential(pot);
      auto job_pot = pot ? pot : work->getPotential();
      eonc::ReplicaExchangeJob job(job_pot, params);
      std::shared_ptr<Matter> out;
      {
        nb::gil_scoped_release release;
        out = job.runFromMatter(work);
      }
      result = matter_result(seed, out, inplace);
      return result;
    }
  };
  nb::class_<PyREX>(m, "ReplicaExchange",
                    "Replica Exchange. run(inplace=False).")
      .def(nb::init<std::shared_ptr<Matter>, Parameters,
                    std::shared_ptr<Potential>>(),
           nb::arg("matter"), nb::arg("parameters"),
           nb::arg("potential") = nb::none(), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 4>())
      .def(
          "run",
          [](PyREX &self, bool inplace) {
            return algo_matter_out(self, self.run(inplace));
          },
          nb::arg("inplace") = false);
}

} // namespace eonc::pybind
