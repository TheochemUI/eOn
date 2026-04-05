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

#include <limits>
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

  // 1. Library Initialization & Configuration (Locked)
  {
    std::lock_guard<std::mutex> lock(res.library_mutex);

    try {
      res.require_loaded();
    } catch (const std::exception &e) {
      QUILL_LOG_ERROR(log, "ARTn library not available: {}", e.what());
      return -1;
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

    const char *filin = "artn_input.dat";
    int result_filin = res.get_set_param_fn()("filin", 0, &size0, filin);
    if (result_filin != 0) {
      QUILL_LOG_ERROR(log, "set_param(filin) failed with code {}",
                      result_filin);
    }

    int verbosity = 3;
    int result_verbose =
        res.get_set_param_fn()("verbose", 0, &size0, &verbosity);
    if (result_verbose != 0) {
      QUILL_LOG_ERROR(log, "set_param(verbose) failed with code {}",
                      result_verbose);
    }

    // Go directly to Lanczos - the mode will be used as the starting guess
    int ninit = 0;
    int result_ninit = res.get_set_param_fn()("ninit", 0, &size0, &ninit);
    if (result_ninit != 0) {
      QUILL_LOG_ERROR(log, "set_param(ninit) failed with code {}",
                      result_ninit);
    }

    // Setup
    bool cerr = false;
    res.get_setup_fn()(nat, &cerr);
    if (cerr) {
      QUILL_LOG_ERROR(log, "ARTn setup failed");
      return -1;
    }
  }

  // Prepare arrays as Eigen objects to avoid std::vector copies
  // Using the memory layout equivalence between RowMajor(N,3) and
  // ColumnMajor(3,N)
  AtomMatrix positions = matter->getPositions();
  AtomMatrix forces(nat, 3);
  forces.setZero(); // Initialize to zero to prevent undefined behavior
  AtomMatrix displacement(nat, 3);
  displacement.setZero(); // Initialize to zero to prevent undefined behavior

  // Prepare Eigen::Maps to provide the correct interface to Fortran
  Eigen::Map<AtomMatrixF> pos_map(positions.data(), 3, nat);
  Eigen::Map<AtomMatrixF> force_map(forces.data(), 3, nat);
  Eigen::Map<AtomMatrixF> disp_map(displacement.data(), 3, nat);
  Eigen::Map<AtomMatrixF> mode_map(mode.data(), 3, nat); // For initial mode

  // Prepare integer arrays (still needed for metadata)
  std::vector<int> ityp(nat);
  std::vector<int> if_pos(3 * nat, 1); // all atoms free by default
  double box_f[9];
  bool lconv = false;
  int dim_mode[2] = {3, nat};

  // Use Eigen for norm calculation
  const double mode_norm = mode.norm();
  const double push_step =
      params.artn_options.push_step_size; // Access push_step here

  if (mode_norm > 1e-10) {
    // Create a temporary mode vector for the initial push
    AtomMatrixF mode_fort = mode_map; // Copy to Fortran layout
    // Vectorized scaling via Eigen Map
    Eigen::Map<VectorXd> mode_vec_map(mode_fort.data(), mode_fort.size());
    mode_vec_map *= (push_step / mode_norm);

    // Set initial push vector for ARTn (locked)
    {
      std::lock_guard<std::mutex> lock(res.library_mutex);
      int result =
          res.get_set_param_fn()("push_init", 2, dim_mode, mode_fort.data());
      if (result != 0) {
        QUILL_LOG_WARNING(log, "set_param(push_init) failed with code {}",
                          result);
      }
    }
  }

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

  for (iteration = 0; iteration < maxIter && !lconv; iteration++) {
    // 2. Parallel PES Evaluation (UNLOCKED)
    // Multiple jthreads can execute this block concurrently.
    double energy = matter->getPotentialEnergy();
    forces = matter->getForces(); // Returns RowMajor Nx3
    this->forcecalls++;

    // 3. ARTn State Update (LOCKED)
    // Serialize only the interaction with the non-thread-safe backend.
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
    }
  }

  // 4. Data Retrieval (LOCKED)
  if (lconv) {
    std::lock_guard<std::mutex> lock(res.library_mutex);

    // Check has_error and has_sad to distinguish success from error
    bool has_error = false;
    bool has_sad = false;
    bool *has_error_ptr = nullptr;
    bool *has_sad_ptr = nullptr;
    int result_has_error = res.get_get_data_fn()(
        "has_error", reinterpret_cast<void **>(&has_error_ptr));
    if (result_has_error != 0) {
      QUILL_LOG_WARNING(log, "get_data(has_error) failed with code {}",
                        result_has_error);
      has_error_ptr = nullptr;
    } else {
      has_error = has_error_ptr ? *has_error_ptr : false;
    }

    int result_has_sad = res.get_get_data_fn()(
        "has_sad", reinterpret_cast<void **>(&has_sad_ptr));
    if (result_has_sad != 0) {
      QUILL_LOG_WARNING(log, "get_data(has_sad) failed with code {}",
                        result_has_sad);
      has_sad_ptr = nullptr;
    } else {
      has_sad = has_sad_ptr ? *has_sad_ptr : false;
    }

    // If a saddle was found, accept it even if force didn't fully converge
    if (has_sad) {
      QUILL_LOG_INFO(
          log,
          "ARTn found saddle after {} iterations (has_error={}, has_sad={})",
          this->iteration, has_error, has_sad);

      // Retrieve eigenvalue
      double *eigval_ptr = nullptr;
      int result_eigval = res.get_get_data_fn()(
          "eigval_sad", reinterpret_cast<void **>(&eigval_ptr));
      if (result_eigval == 0 && eigval_ptr) {
        eigenvalue = *eigval_ptr;
      } else {
        QUILL_LOG_WARNING(
            log, "Failed to retrieve eigenvalue (result={}, ptr_valid={})",
            result_eigval, eigval_ptr != nullptr);
        eigenvalue =
            std::numeric_limits<double>::quiet_NaN(); // Use NaN to indicate
                                                      // missing value
      }

      // Retrieve eigenvector
      double *evec_ptr = nullptr;
      int result_evec = res.get_get_data_fn()(
          "eigen_sad", reinterpret_cast<void **>(&evec_ptr));
      if (result_evec == 0 && evec_ptr) {
        // Use direct Eigen::Map to convert from Fortran layout
        eigenvector = eonc::from_fortran_layout_vector(
            std::vector<double>(evec_ptr, evec_ptr + 3 * nat), nat);

        // WARNING: Do NOT free evec_ptr - ARTn may have allocated it with
        // Fortran ALLOCATE Calling std::free() on Fortran-allocated memory
        // causes segfaults Memory must be freed by ARTn library via exposed
        // deallocation function For now, we must rely on ARTn to manage this
        // memory internally
      } else {
        QUILL_LOG_ERROR(
            log, "Failed to retrieve eigenvector (result={}, ptr_valid={})",
            result_evec, evec_ptr != nullptr);
        // Set eigenvector to zero matrix if retrieval failed
        eigenvector = AtomMatrix::Zero(nat, 3);
      }

      status = 0;           // GOOD
      res.get_clean_fn()(); // Clean up before releasing lock
      return 0;
    }

    // No saddle found - this is a real error
    QUILL_LOG_WARNING(
        log, "ARTn stopped after {} iterations (has_error={}, has_sad={})",
        iteration, has_error, has_sad);
    status = 2;           // BAD_ARTN_ERROR
    res.get_clean_fn()(); // Clean up before releasing lock
    return 2;
  }

  QUILL_LOG_WARNING(log, "ARTn did not converge after {} iterations",
                    iteration);
  status = 1; // BAD_MAX_ITERATIONS

  // Clean up in all cases
  {
    std::lock_guard<std::mutex> lock(res.library_mutex);
    res.get_clean_fn()();
  }
  return 1;

#else
  QUILL_LOG_ERROR(log, "ARTn support not compiled");
  return -1;
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

} // namespace eonc
