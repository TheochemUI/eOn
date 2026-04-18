/*
 * This file is part of eOn.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Copyright (c) 2010--present, eOn Development Team
 * All rights reserved.
 *
 * Repo:
 * https://github.com/TheochemUI/eOn
 */
#include "ARTnSaddleSearch.h"
#include "Eigen.h"

#include <cstdlib>
#include <filesystem>
#include <limits>
#include <sstream>
namespace eonc {

ARTnSaddleSearch::ARTnSaddleSearch(std::shared_ptr<Matter> matterPassed,
                                   std::shared_ptr<Potential> potPassed,
                                   AtomMatrix modeInitial,
                                   const Parameters &paramsPassed)
    : SaddleSearchMethod(potPassed, paramsPassed),
      matter{matterPassed},
      mode{modeInitial},
      eigenvector{AtomMatrix::Zero(matterPassed->numberOfAtoms(), 3)} {
  log = eonc::log::get();
  if (!log) {
    throw std::runtime_error("ARTnSaddleSearch: Logger not initialized");
  }
}

ARTnSaddleSearch::~ARTnSaddleSearch() {
#ifdef WITH_ARTN
  // Clean up is done within the search loop, not in destructor
#endif
}

int ARTnSaddleSearch::run() {
#ifdef WITH_ARTN
  auto &res = get_artn_resource();
  const int nat = matter->numberOfAtoms();

  if (mode.rows() != nat || mode.cols() != 3) {
    mode = AtomMatrix::Zero(nat, 3);
  }

  // Pre-declare force-loop storage outside the lock so it outlives the
  // setup critical section. These reads do not touch pARTn state.
  AtomMatrix positions = matter->getPositions();
  AtomMatrix forces = AtomMatrix::Zero(nat, 3);
  AtomMatrix displacement = AtomMatrix::Zero(nat, 3);

  // Eigen::Maps give Fortran a column-major [3, nat] view over the same
  // memory as the row-major [nat, 3] AtomMatrix (zero-copy).
  Eigen::Map<AtomMatrixF> pos_map(positions.data(), 3, nat);
  Eigen::Map<AtomMatrixF> force_map(forces.data(), 3, nat);
  Eigen::Map<AtomMatrixF> disp_map(displacement.data(), 3, nat);
  Eigen::Map<AtomMatrixF> mode_map(mode.data(), 3, nat);

  const double push_step = params.artn_options.push_step_size;
  const double mode_norm = mode.norm();
  AtomMatrixF mode_fort;
  if (mode_norm > 1e-10) {
    mode_fort = mode_map;
    Eigen::Map<VectorXd> mode_vec_map(mode_fort.data(), mode_fort.size());
    mode_vec_map *= (push_step / mode_norm);
  }
  int dim_mode[2] = {3, nat};

  // 1. Library Initialization, Configuration, and Initial Push (Locked)
  //    Held as a single critical section so concurrent ARTn searches in
  //    the same process cannot interleave create(), set_param(), setup,
  //    or push_init calls on pARTn's non-thread-safe global state.
  {
    std::lock_guard<std::mutex> lock(res.library_mutex);

    try {
      res.require_loaded();
    } catch (const std::exception &e) {
      QUILL_LOG_ERROR(log, "ARTn library not available: {}", e.what());
      status = STATUS_BAD_ARTN_ERROR;
      return status;
    }

    res.get_create_fn()();

    // Set engine units first (required before any other params)
    const char *units = "lammps/metal";
    int size0 = 0;
    int result_units = res.get_set_param_fn()("engine_units", 0, &size0, units);
    if (result_units != 0) {
      QUILL_LOG_ERROR(log, "set_param(engine_units) failed with code {}",
                      result_units);
    }

    // Set parameters from eOn config
    double push_step = params.artn_options.push_step_size;
    int result_push =
        res.get_set_param_fn()("push_step_size", 0, &size0, &push_step);
    if (result_push != 0) {
      QUILL_LOG_ERROR(log, "set_param(push_step_size) failed with code {}",
                      result_push);
    }

    double force_thr = params.artn_options.force_threshold;
    int result_force =
        res.get_set_param_fn()("forc_thr", 0, &size0, &force_thr);
    if (result_force != 0) {
      QUILL_LOG_ERROR(log, "set_param(forc_thr) failed with code {}",
                      result_force);
    }

    // Only set filin if the input file exists; pARTn fails setup if the
    // file is specified but missing.
    const char *filin = "artn_input.dat";
    if (std::filesystem::exists(filin)) {
      int result_filin = res.get_set_param_fn()("filin", 0, &size0, filin);
      if (result_filin != 0) {
        QUILL_LOG_ERROR(log, "set_param(filin) failed with code {}",
                        result_filin);
      }
    }

    int verbosity = 3;
    int result_verbose =
        res.get_set_param_fn()("verbose", 0, &size0, &verbosity);
    if (result_verbose != 0) {
      QUILL_LOG_ERROR(log, "set_param(verbose) failed with code {}",
                      result_verbose);
    }

    // ninit controls initial push steps before Lanczos eigenmode estimation.
    // 0 = skip push, go straight to Lanczos (appropriate when eOn provides
    // the displacement direction via push_init); >0 = push that many steps.
    // -1 sentinel means "leave pARTn's own default in place", so we only
    // call set_param when the user asked for a specific value.
    if (params.artn_options.ninit >= 0) {
      int ninit = params.artn_options.ninit;
      int result_ninit = res.get_set_param_fn()("ninit", 0, &size0, &ninit);
      if (result_ninit != 0) {
        QUILL_LOG_ERROR(log, "set_param(ninit) failed with code {}",
                        result_ninit);
      }
    }

    // nperp_limitation: controls perp-relax steps per Lanczos cycle.
    // pARTn defaults are tuned for exploration from minimum. For refinement
    // near a saddle, -1 (unlimited) or 20-30 (for ML potentials) is better.
    if (params.artn_options.nperp_limitation != "default") {
      // Parse comma-separated integers into a vector
      std::vector<int> nperp_vals;
      std::istringstream ss(params.artn_options.nperp_limitation);
      std::string token;
      while (std::getline(ss, token, ',')) {
        nperp_vals.push_back(std::stoi(token));
      }
      if (!nperp_vals.empty()) {
        int nperp_size = static_cast<int>(nperp_vals.size());
        int result_nperp = res.get_set_param_fn()(
            "nperp_limitation", 1, &nperp_size, nperp_vals.data());
        if (result_nperp != 0) {
          QUILL_LOG_WARNING(log, "set_param(nperp_limitation) failed: {}",
                            result_nperp);
        }
      }
    }

    // lanczos_min_size: minimum Lanczos iterations before convergence check.
    // Default 3 for exploration; 1 for refinement near a saddle.
    if (params.artn_options.lanczos_min_size >= 0) {
      int lms = params.artn_options.lanczos_min_size;
      res.get_set_param_fn()("lanczos_min_size", 0, &size0, &lms);
    }

    // nsmooth: number of smooth interpolation steps. 0 disables.
    if (params.artn_options.nsmooth >= 0) {
      int ns = params.artn_options.nsmooth;
      res.get_set_param_fn()("nsmooth", 0, &size0, &ns);
    }

    // nnewchance: retries permitted when Lanczos returns a positive lowest
    // eigenvalue (convex region, no unstable mode). pARTn defaults to 0,
    // i.e. immediate "EIGENVALUE LOST" failure -- too brittle for small
    // clusters where the eigenvalue can flip positive transiently before
    // the saddle direction stabilizes. eOn's default of 3 (Parameters.h)
    // gives Lanczos a few random-restart chances before giving up.
    if (params.artn_options.nnewchance >= 0) {
      int nnc = params.artn_options.nnewchance;
      res.get_set_param_fn()("nnewchance", 0, &size0, &nnc);
    }

    // Setup
    bool cerr = false;
    res.get_setup_fn()(nat, &cerr);
    if (cerr) {
      QUILL_LOG_ERROR(log, "ARTn setup failed (nat={})", nat);
      res.get_destroy_fn()();
      status = STATUS_BAD_ARTN_ERROR;
      return status;
    }

    // Initial push vector (if eOn supplied a non-trivial mode). Must be set
    // after setup_artn() and before the first artn_step(); keep inside the
    // same critical section so no concurrent ARTn search can re-init between
    // setup and push.
    if (mode_norm > 1e-10) {
      int result =
          res.get_set_param_fn()("push_init", 2, dim_mode, mode_fort.data());
      if (result != 0) {
        QUILL_LOG_WARNING(log, "set_param(push_init) failed with code {}",
                          result);
      }
    }
  }

  // Per-atom metadata for the Fortran step (no pARTn state, unlocked).
  std::vector<int> ityp(nat);
  std::vector<int> if_pos(3 * nat, 1); // all atoms free by default
  double box_f[9];
  bool lconv = false;

  for (int i = 0; i < nat; i++) {
    if (matter->getFixed(i)) {
      Eigen::Map<Eigen::Vector3i>(&if_pos[i * 3]).setZero();
    }
    ityp[i] = matter->getAtomicNr(i);
  }

  // Convert box to column-major 3x3 - copy to array
  Matrix3d cell = matter->getCell();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      box_f[j * 3 + i] = cell(i, j);

  int maxIter = params.artn_options.max_iterations;

  // while-loop (not for) so `iteration` reports the index of the converged
  // step rather than one past it: on convergence during step k, we break
  // before the increment and the post-loop value is k, matching log output.
  iteration = 0;
  while (iteration < maxIter && !lconv) {
    // 2. Parallel PES Evaluation (UNLOCKED)
    // Concurrent ARTn searches in the same process would share the same
    // PES call, so the potential evaluation itself is not serialized here.
    double energy = matter->getPotentialEnergy();
    forces = matter->getForces(); // Returns RowMajor Nx3
    this->forcecalls++;

    // 3. ARTn State Update (LOCKED)
    // Serialize only the interaction with the non-thread-safe backend.
    //
    // Perf note (artn-plugin >= 9dab2053): the inner Lanczos eigenvector
    // reconstruction now uses intrinsic matmul on a reshaped Vmat slice,
    // which allocates a temporary [3*nat, ilanc] array per Lanczos iteration.
    // Negligible at our sizes (small molecules, ilanc < O(20)); revisit if
    // we ever drive artn against large supercell DFT.
    {
      std::unique_lock<std::mutex> lock(res.library_mutex);
      res.get_artn_step_fn()(nat, energy, force_map.data(), ityp.data(),
                             pos_map.data(), box_f, if_pos.data(),
                             disp_map.data(), &lconv);
    }

    if (!lconv) {
      // Apply displacement and update Matter (PES-specific, typically
      // thread-safe) Add displacement to positions
      positions += displacement;
      matter->setPositions(positions);
      iteration++;
    }
  }

  // 4. Data Retrieval (LOCKED)
  if (lconv) {
    std::lock_guard<std::mutex> lock(res.library_mutex);

    // pARTn exposes a dedicated C get_error() that returns both the error
    // code and a c_malloc'd message pointer (see m_artn_error.f90 in
    // artn-plugin). Prefer it when available; fall back to the has_error
    // flag via get_data for older libartn builds that predate the wrapper.
    int artn_err = 0;
    std::string artn_err_msg;
    if (auto *get_error_fn_ = res.get_get_error_fn()) {
      void *cmsg = nullptr;
      artn_err = get_error_fn_(&cmsg);
      if (artn_err != 0 && cmsg != nullptr) {
        artn_err_msg.assign(static_cast<const char *>(cmsg));
        std::free(cmsg);
      }
    } else {
      bool *has_error_ptr = nullptr;
      int result_has_error = res.get_get_data_fn()(
          "has_error", reinterpret_cast<void **>(&has_error_ptr));
      if (result_has_error != 0) {
        QUILL_LOG_WARNING(log, "get_data(has_error) failed with code {}",
                          result_has_error);
      } else if (has_error_ptr) {
        artn_err = *has_error_ptr ? 1 : 0;
        std::free(has_error_ptr);
      }
    }
    bool has_error = artn_err != 0;

    bool has_sad = false;
    bool *has_sad_ptr = nullptr;
    int result_has_sad = res.get_get_data_fn()(
        "has_sad", reinterpret_cast<void **>(&has_sad_ptr));
    if (result_has_sad != 0) {
      QUILL_LOG_WARNING(log, "get_data(has_sad) failed with code {}",
                        result_has_sad);
    } else if (has_sad_ptr) {
      has_sad = *has_sad_ptr;
      std::free(has_sad_ptr);
    }

    if (has_error && !artn_err_msg.empty()) {
      QUILL_LOG_WARNING(log, "pARTn reported error {}: {}", artn_err,
                        artn_err_msg);
    }

    // If a saddle was found, accept it even if force didn't fully converge
    if (has_sad) {
      QUILL_LOG_INFO(
          log,
          "ARTn found saddle after {} iterations (has_error={}, has_sad={})",
          this->iteration, has_error, has_sad);

      // Retrieve the saddle coordinates tracked internally by pARTn so the
      // Matter object matches the reported eigenpair and subsequent endpoint
      // minimizations start from the actual saddle.
      double *tau_sad_ptr = nullptr;
      int result_tau_sad = res.get_get_data_fn()(
          "tau_sad", reinterpret_cast<void **>(&tau_sad_ptr));
      if (result_tau_sad == 0 && tau_sad_ptr) {
        matter->setPositions(eonc::from_fortran_layout_vector(
            std::vector<double>(tau_sad_ptr, tau_sad_ptr + 3 * nat), nat));
        std::free(tau_sad_ptr);
      } else {
        QUILL_LOG_WARNING(
            log, "Failed to retrieve tau_sad (result={}, ptr_valid={})",
            result_tau_sad, tau_sad_ptr != nullptr);
      }

      // Retrieve eigenvalue
      double *eigval_ptr = nullptr;
      int result_eigval = res.get_get_data_fn()(
          "eigval_sad", reinterpret_cast<void **>(&eigval_ptr));
      if (result_eigval == 0 && eigval_ptr) {
        eigenvalue = *eigval_ptr;
        std::free(eigval_ptr);
      } else {
        QUILL_LOG_WARNING(
            log, "Failed to retrieve eigenvalue (result={}, ptr_valid={})",
            result_eigval, eigval_ptr != nullptr);
        eigenvalue =
            std::numeric_limits<double>::quiet_NaN(); // Use NaN to indicate
                                                      // missing value
      }

      // Retrieve eigenvector (3*nat flat array, column-major from Fortran).
      // get_data allocates via c_malloc and writes the pointer to cval.
      // The C header says void* but the Fortran intent(out) semantics
      // require void** (see artn_c_wrappers.f90:324 and LAMMPS example).
      double *evec_ptr = nullptr;
      int result_evec = res.get_get_data_fn()(
          "eigen_sad", reinterpret_cast<void **>(&evec_ptr));
      if (result_evec == 0 && evec_ptr) {
        // Use direct Eigen::Map to convert from Fortran layout
        eigenvector = eonc::from_fortran_layout_vector(
            std::vector<double>(evec_ptr, evec_ptr + 3 * nat), nat);

        // get_data allocates via c_malloc (artn_c_wrappers.f90), safe to free
        std::free(evec_ptr);
      } else {
        QUILL_LOG_ERROR(
            log, "Failed to retrieve eigenvector (result={}, ptr_valid={})",
            result_evec, evec_ptr != nullptr);
        // Set eigenvector to zero matrix if retrieval failed
        eigenvector = AtomMatrix::Zero(nat, 3);
      }

      status = STATUS_GOOD;
      res.get_destroy_fn()();
      return status;
    }

    // No saddle found - this is a real error
    QUILL_LOG_WARNING(
        log, "ARTn stopped after {} iterations (has_error={}, has_sad={})",
        iteration, has_error, has_sad);
    status = STATUS_BAD_ARTN_ERROR;
    res.get_destroy_fn()();
    return status;
  }

  QUILL_LOG_WARNING(log, "ARTn did not converge after {} iterations",
                    iteration);
  status = STATUS_BAD_MAX_ITERATIONS;

  // Clean up in all cases
  {
    std::lock_guard<std::mutex> lock(res.library_mutex);
    res.get_destroy_fn()();
  }
  return status;

#else
  QUILL_LOG_ERROR(log, "ARTn support not compiled");
  status = STATUS_BAD_ARTN_ERROR;
  return status;
#endif
}

double ARTnSaddleSearch::getEigenvalue() {
  if (log && std::isnan(eigenvalue)) {
    QUILL_LOG_WARNING(log, "Requesting uninitialized/invalid eigenvalue");
  }
  return eigenvalue;
}

AtomMatrix ARTnSaddleSearch::getEigenvector() {
  if (log && std::isnan(eigenvalue)) {
    QUILL_LOG_WARNING(log, "Requesting uninitialized eigenvector");
  }
  return eigenvector;
}

std::string_view ARTnSaddleSearch::describeStatus(int status) const {
  switch (status) {
  case STATUS_GOOD:
    return "Success";
  case STATUS_BAD_MAX_ITERATIONS:
    return "Too many iterations";
  case STATUS_BAD_ARTN_ERROR:
    return "ARTn backend error";
  default:
    return "Unknown status";
  }
}

} // namespace eonc
