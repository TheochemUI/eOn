/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "MetatomicPotential.h"
#include "../../Parameters.h"
#include "../../fpe_handler.h"
#include "vesin.h"

#include <torch/csrc/jit/runtime/graph_executor.h>

#include <cstdint>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std::string_literals;

static torch::optional<std::string> normalize_variant(const std::string &s) {
  if (s.empty() || s == "off")
    return torch::nullopt;
  return s;
}

MetatomicPotential::MetatomicPotential(const Parameters &params)
    : Potential(PotType::METATOMIC),
      m_metatomic_opts{params.metatomic_options},
      model_(torch::jit::Module()),
      device_type_(c10::DeviceType::CPU),
      device_(torch::Device(device_type_)) {

  // Determinism knobs (see
  // https://rgoswami.me/snippets/pytorch-deterministic-regression/): JIT
  // profiling specializes graphs after the first few forwards, so a fresh model
  // instance can differ at ULP level from a warm one — bad for parallel NEB.
  // cuBLAS / index_add_ on CUDA are likewise nondeterministic unless forced.
  // [Metatomic] deterministic=true (default) applies the safe defaults;
  // deterministic_strict=true requires CUBLAS_WORKSPACE_CONFIG (e.g. :4096:8)
  // and fails on nondeterministic ops instead of warning.
  if (m_metatomic_opts.deterministic) {
    torch::jit::getProfilingMode() = false;
    const bool strict = m_metatomic_opts.deterministic_strict ||
                        (std::getenv("CUBLAS_WORKSPACE_CONFIG") != nullptr);
    at::globalContext().setDeterministicAlgorithms(true,
                                                   /*warn_only=*/!strict);
    at::globalContext().setBenchmarkCuDNN(false);
    QUILL_LOG_INFO(m_log,
                   "[MetatomicPotential] Deterministic algorithms enabled "
                   "(strict={})",
                   strict);
  } else {
    QUILL_LOG_INFO(m_log,
                   "[MetatomicPotential] Deterministic algorithms disabled "
                   "(faster, may diverge across runs / NEB images)");
  }

  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();

  QUILL_LOG_INFO(m_log, "[MetatomicPotential] Initializing...");

  // 1. Load the model from the path specified in parameters
  torch::optional<std::string> extensions_directory = torch::nullopt;
  if (!m_metatomic_opts.extensions_directory.empty()) {
    extensions_directory = m_metatomic_opts.extensions_directory;
  }

  try {
    QUILL_LOG_INFO(m_log, "[MetatomicPotential] Loading model from '{}'",
                   m_metatomic_opts.model_path);
    this->model_ = metatomic_torch::load_atomistic_model(
        m_metatomic_opts.model_path, extensions_directory);
  } catch (const std::exception &e) {
    QUILL_LOG_ERROR(m_log, "[MetatomicPotential] Failed to load model: {}",
                    e.what());
    throw;
  }

  // 2. Extract capabilities and neighbor list requests from the model
  this->capabilities_ =
      this->model_.run_method("capabilities")
          .toCustomClass<metatomic_torch::ModelCapabilitiesHolder>();
  auto requests_ivalue = this->model_.run_method("requested_neighbor_lists");
  for (const auto &request_ivalue : requests_ivalue.toList()) {
    auto request =
        request_ivalue.get()
            .toCustomClass<metatomic_torch::NeighborListOptionsHolder>();
    this->nl_requests_.push_back(request);
  }

  // 3. Determine and set up the device (CPU/CUDA/MPS)
  torch::optional<std::string> desired = torch::nullopt;
  if (!m_metatomic_opts.device.empty()) {
    desired = m_metatomic_opts.device;
  }
  device_type_ = metatomic_torch::pick_device(
      this->capabilities_->supported_devices, desired);

  device_ = torch::Device(device_type_);
  QUILL_LOG_INFO(m_log, "[MetatomicPotential] Using device: {}", device_.str());

  this->model_.to(this->device_);

  // 4. Set data type (float32/float64) based on model capabilities
  if (this->capabilities_->dtype() == "float64") {
    this->dtype_ = torch::kFloat64;
  } else if (this->capabilities_->dtype() == "float32") {
    this->dtype_ = torch::kFloat32;
  } else {
    throw std::runtime_error("Unsupported dtype: " +
                             this->capabilities_->dtype());
  }
  QUILL_LOG_INFO(m_log, "[MetatomicPotential] Using dtype: {}",
                 this->capabilities_->dtype().c_str());

  // 5. Resolve energy / force output keys: explicit keys (#215) or variants
  // (#296 for non_conservative_force)
  auto outputs = this->capabilities_->outputs();

  auto v_base = normalize_variant(m_metatomic_opts.variant.base);
  auto v_energy = m_metatomic_opts.variant.energy.empty()
                      ? v_base
                      : normalize_variant(m_metatomic_opts.variant.energy);
  auto v_energy_uq =
      m_metatomic_opts.variant.energy_uncertainty.empty()
          ? v_energy
          : normalize_variant(m_metatomic_opts.variant.energy_uncertainty);
  auto v_force = m_metatomic_opts.variant.force.empty()
                     ? v_energy
                     : normalize_variant(m_metatomic_opts.variant.force);

  if (!m_metatomic_opts.energy_output.empty()) {
    this->energy_key_ = m_metatomic_opts.energy_output;
  } else {
    this->energy_key_ =
        metatomic_torch::pick_output("energy", outputs, v_energy);
  }

  this->non_conservative_ = m_metatomic_opts.non_conservative;
  this->random_rotation_ = m_metatomic_opts.random_rotation;
  this->n_symmetry_rotations_ = m_metatomic_opts.n_symmetry_rotations;
  if (this->non_conservative_) {
    if (!m_metatomic_opts.force_output.empty()) {
      this->nc_forces_key_ = m_metatomic_opts.force_output;
      if (!outputs.contains(this->nc_forces_key_)) {
        throw std::runtime_error(
            "Missing explicit force_output in metatomic model: " +
            this->nc_forces_key_);
      }
    } else {
      this->nc_forces_key_ = metatomic_torch::pick_output(
          "non_conservative_force", outputs, v_force);
    }
    QUILL_LOG_INFO(m_log,
                   "[MetatomicPotential] Non-conservative forces from '{}'",
                   this->nc_forces_key_);
  }
  if (this->n_symmetry_rotations_ > 0) {
    QUILL_LOG_INFO(m_log,
                   "[MetatomicPotential] Symmetry averaging over {} rotations",
                   this->n_symmetry_rotations_);
  } else if (this->random_rotation_) {
    QUILL_LOG_INFO(
        m_log, "[MetatomicPotential] Per-call random SO(3) rotation enabled");
  }

  if (!outputs.contains(this->energy_key_)) {
    QUILL_LOG_ERROR(
        m_log,
        "[MetatomicPotential] The model does not provide an '{}' output.",
        this->energy_key_);
    throw std::runtime_error("Missing energy output in metatomic model");
  }

  // 6. Set up evaluation options to request total energy
  this->evaluations_options_ =
      torch::make_intrusive<metatomic_torch::ModelEvaluationOptionsHolder>();
  evaluations_options_->set_length_unit(m_metatomic_opts.length_unit);

  auto model_output = outputs.at(this->energy_key_);
  auto requested_output =
      torch::make_intrusive<metatomic_torch::ModelOutputHolder>();

  // Per-atom granularity is sample_kind == "atom" (get/set_per_atom removed).
  requested_output->set_sample_kind(model_output->sample_kind());
  requested_output->explicit_gradients = {};
  requested_output->set_unit("eV");
  evaluations_options_->outputs.insert(this->energy_key_, requested_output);

  // Request non-conservative forces when enabled (#296)
  if (this->non_conservative_ && !this->nc_forces_key_.empty()) {
    auto nc_info = outputs.at(this->nc_forces_key_);
    auto requested_nc =
        torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
    requested_nc->set_sample_kind(nc_info->sample_kind());
    requested_nc->explicit_gradients = {};
    requested_nc->set_unit("eV/Angstrom");
    evaluations_options_->outputs.insert(this->nc_forces_key_, requested_nc);
  }

  // 7. Optionally request energy uncertainty if threshold is positive
  if (m_metatomic_opts.uncertainty_threshold > 0) {
    this->uncertainty_threshold_ = m_metatomic_opts.uncertainty_threshold;
    const bool explicit_uq_key =
        !m_metatomic_opts.energy_uncertainty_output.empty();
    if (explicit_uq_key) {
      // User-specified key: hard fail if missing (not soft-disabled).
      this->energy_uncertainty_key_ =
          m_metatomic_opts.energy_uncertainty_output;
      if (!outputs.contains(this->energy_uncertainty_key_)) {
        QUILL_LOG_ERROR(m_log,
                        "[MetatomicPotential] energy_uncertainty_output '{}' "
                        "is not provided by the model.",
                        this->energy_uncertainty_key_);
        throw std::runtime_error(
            "Missing explicit energy_uncertainty_output in metatomic model: " +
            this->energy_uncertainty_key_);
      }
    } else {
      try {
        this->energy_uncertainty_key_ = metatomic_torch::pick_output(
            "energy_uncertainty", outputs, v_energy_uq);
      } catch (const std::exception &e) {
        QUILL_LOG_DEBUG(
            m_log, "[MetatomicPotential] No uncertainty output available: {}",
            e.what());
        this->uncertainty_threshold_ = -1.0;
      }
    }

    if (this->uncertainty_threshold_ > 0) {
      auto uncertainty_info = outputs.at(this->energy_uncertainty_key_);
      if (uncertainty_info->sample_kind() == "atom") {
        auto requested_uncertainty =
            torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
        requested_uncertainty->set_sample_kind("atom");
        requested_uncertainty->explicit_gradients = {};
        requested_uncertainty->set_unit("eV");
        evaluations_options_->outputs.insert(this->energy_uncertainty_key_,
                                             requested_uncertainty);
        QUILL_LOG_INFO(m_log,
                       "[MetatomicPotential] Requested per-atom "
                       "'{}' from model (threshold = {})",
                       this->energy_uncertainty_key_,
                       this->uncertainty_threshold_);
      } else {
        QUILL_LOG_DEBUG(m_log,
                        "[MetatomicPotential] Model provides '{}' "
                        "but sample_kind is not \"atom\"; skipping uncertainty "
                        "checks.",
                        this->energy_uncertainty_key_);
        this->uncertainty_threshold_ = -1.0;
      }
    }
  }

  this->check_consistency_ = m_metatomic_opts.check_consistency;
  QUILL_LOG_INFO(m_log, "[MetatomicPotential] Initialization complete.");

  fpeh.restore_fpe();
}

// --- helpers for random / symmetry rotations (#287, #292) ---

namespace {

// Uniform random rotation in SO(3) via QR of a Gaussian matrix with positive
// determinant (Arvo / Shoemake style, sufficient for stochastic averaging).
torch::Tensor random_so3(torch::Device device, torch::ScalarType dtype) {
  auto A =
      torch::randn({3, 3}, torch::TensorOptions().dtype(dtype).device(device));
  auto qr = torch::linalg_qr(A);
  auto Q = std::get<0>(qr);
  auto R = std::get<1>(qr);
  auto d = torch::sign(torch::diagonal(R));
  Q = Q * d.unsqueeze(0);
  if (torch::det(Q).item<double>() < 0) {
    Q.select(1, 0).mul_(-1);
  }
  return Q;
}

} // namespace

// --- MetatomicPotential::force ---

void MetatomicPotential::force(long nAtoms, const double *positions,
                               const int *atomicNrs, double *forces,
                               double *energy, double *variance,
                               const double *box) {
  // Serialize concurrent calls -- PyTorch model inference on the same
  // Module instance is not thread-safe
  std::lock_guard<std::mutex> lock(inference_mutex_);

  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();

  if (!atomicNrs) {
    throw std::runtime_error(
        "[MetatomicPotential] `atomicNrs` must be provided.");
  }

  const long n_avg = this->n_symmetry_rotations_ > 0
                         ? this->n_symmetry_rotations_
                         : (this->random_rotation_ ? 1 : 1);
  const bool use_rotation =
      this->random_rotation_ || this->n_symmetry_rotations_ > 0;
  // n_symmetry_rotations averages; random_rotation alone is one rotated eval
  const long n_passes =
      this->n_symmetry_rotations_ > 0 ? this->n_symmetry_rotations_ : 1;

  auto f64_options =
      torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
  std::vector<int32_t> types_vec(atomicNrs, atomicNrs + nAtoms);
  auto atomic_types_cpu =
      torch::tensor(types_vec, torch::TensorOptions().dtype(torch::kInt32));

  double energy_acc = 0.0;
  auto forces_acc = torch::zeros({nAtoms, 3}, f64_options);
  bool variance_set = false;

  for (long i_pass = 0; i_pass < n_passes; ++i_pass) {
    torch::Tensor R = torch::eye(
        3, torch::TensorOptions().dtype(this->dtype_).device(this->device_));
    if (use_rotation) {
      R = random_so3(this->device_, this->dtype_);
    }
    // R is applied to row vectors: pos' = pos @ R^T  (equiv. R @ pos for cols)
    auto R_cpu = R.to(torch::kCPU).to(torch::kFloat64);
    auto R_T = R.transpose(0, 1);

    auto pos_cpu = torch::from_blob(const_cast<double *>(positions),
                                    {nAtoms, 3}, f64_options)
                       .clone();
    auto cell_cpu =
        torch::from_blob(const_cast<double *>(box), {3, 3}, f64_options)
            .clone();
    if (use_rotation) {
      pos_cpu = pos_cpu.matmul(R_cpu.transpose(0, 1));
      // Rotate cell vectors (rows) the same way
      cell_cpu = cell_cpu.matmul(R_cpu.transpose(0, 1));
    }

    std::vector<double> pos_buf(static_cast<size_t>(nAtoms) * 3);
    std::vector<double> cell_buf(9);
    std::memcpy(pos_buf.data(), pos_cpu.contiguous().data_ptr<double>(),
                pos_buf.size() * sizeof(double));
    std::memcpy(cell_buf.data(), cell_cpu.contiguous().data_ptr<double>(),
                9 * sizeof(double));

    auto torch_positions =
        torch::from_blob(pos_buf.data(), {nAtoms, 3}, f64_options)
            .to(this->dtype_)
            .to(this->device_)
            .set_requires_grad(!this->non_conservative_);

    auto torch_cell = torch::from_blob(cell_buf.data(), {3, 3}, f64_options)
                          .to(this->dtype_)
                          .to(this->device_);

    auto cell_norms = torch::norm(torch_cell, 2, /*dim=*/1);
    auto torch_pbc = cell_norms.abs() > 1e-9;
    bool periodic[3] = {torch_pbc[0].item<bool>(), torch_pbc[1].item<bool>(),
                        torch_pbc[2].item<bool>()};

    auto atomic_types = atomic_types_cpu.to(this->device_);

    auto system = torch::make_intrusive<metatomic_torch::SystemHolder>(
        atomic_types, torch_positions, torch_cell, torch_pbc);

    for (const auto &request : this->nl_requests_) {
      auto neighbors = this->computeNeighbors(request, nAtoms, pos_buf.data(),
                                              cell_buf.data(), periodic);
      metatomic_torch::register_autograd_neighbors(system, neighbors,
                                                   this->check_consistency_);
      system->add_neighbor_list(request, neighbors);
    }

    torch::Tensor forces_tensor;
    try {
      auto ivalue_output = this->model_.forward({
          std::vector<metatomic_torch::System>{system},
          evaluations_options_,
          this->check_consistency_,
      });
      auto dict_output = ivalue_output.toGenericDict();
      auto output_map = dict_output.at(this->energy_key_)
                            .toCustomClass<metatensor_torch::TensorMapHolder>();

      if (this->uncertainty_threshold_ > 0 && i_pass == 0) {
        try {
          if (dict_output.contains(this->energy_uncertainty_key_)) {
            auto uncertainty_map =
                dict_output.at(this->energy_uncertainty_key_)
                    .toCustomClass<metatensor_torch::TensorMapHolder>();
            auto uncertainty_block =
                metatensor_torch::TensorMapHolder::block_by_id(uncertainty_map,
                                                               0);
            auto flat_uncertainty =
                uncertainty_block->values().reshape({-1}).to(torch::kCPU);
            if (variance != nullptr && flat_uncertainty.numel() > 0) {
              try {
                *variance =
                    flat_uncertainty.to(torch::kFloat64).mean().item<double>();
                variance_set = true;
              } catch (...) {
                QUILL_LOG_DEBUG(m_log,
                                "[MetatomicPotential] Failed to compute mean "
                                "uncertainty for variance.");
              }
            }
            auto atoms_above_threshold =
                flat_uncertainty > this->uncertainty_threshold_;
            if (torch::any(atoms_above_threshold).item<bool>()) {
              auto samples = uncertainty_block->samples();
              auto atom_indices_all = samples->column("atom").to(torch::kCPU);
              auto atom_indices_above =
                  atom_indices_all.index({atoms_above_threshold});
              std::ostringstream ss;
              ss << "atoms at index [";
              auto n_report = std::min<int64_t>(10, atom_indices_above.size(0));
              for (int64_t i = 0; i < n_report; ++i) {
                if (i > 0)
                  ss << ", ";
                ss << atom_indices_above[i].item<int32_t>();
              }
              ss << "]";
              if (atom_indices_above.size(0) > n_report) {
                ss << " and " << (atom_indices_above.size(0) - n_report)
                   << " more";
              }
              QUILL_LOG_WARNING(
                  m_log,
                  "[MetatomicPotential] The uncertainty on atomic energies for "
                  "{} are larger than the threshold of {}. (Key: {}) Be "
                  "careful "
                  "when analyzing the results, and consider retraining the "
                  "model to better describe these configurations.",
                  ss.str(), this->uncertainty_threshold_,
                  this->energy_uncertainty_key_);
            }
          }
        } catch (const std::exception &e) {
          QUILL_LOG_WARNING(m_log,
                            "[MetatomicPotential] Failed to check {}: {}",
                            this->energy_uncertainty_key_, e.what());
        }
      }

      auto energy_block =
          metatensor_torch::TensorMapHolder::block_by_id(output_map, 0);
      auto energy_tensor = energy_block->values();
      energy_acc += energy_tensor.sum().item<double>();

      if (this->non_conservative_ && !this->nc_forces_key_.empty()) {
        auto nc_map = dict_output.at(this->nc_forces_key_)
                          .toCustomClass<metatensor_torch::TensorMapHolder>();
        auto nc_block =
            metatensor_torch::TensorMapHolder::block_by_id(nc_map, 0);
        forces_tensor = nc_block->values()
                            .reshape({nAtoms, 3})
                            .to(torch::kCPU)
                            .to(torch::kFloat64);
      } else {
        energy_tensor.backward(torch::ones_like(energy_tensor));
        auto positions_grad = system->positions().grad();
        forces_tensor = (-positions_grad).to(torch::kCPU).to(torch::kFloat64);
      }
    } catch (const std::exception &e) {
      QUILL_LOG_ERROR(m_log, "[MetatomicPotential] Model evaluation failed: {}",
                      e.what());
      throw;
    }

    // Rotate forces back to original frame: F = F' @ R  (since pos' = pos @
    // R^T)
    if (use_rotation) {
      forces_tensor = forces_tensor.matmul(R_cpu);
    }
    forces_acc += forces_tensor;
  }

  const double inv_n = 1.0 / static_cast<double>(n_passes);
  *energy = energy_acc * inv_n;
  forces_acc = forces_acc * inv_n;
  (void)variance_set;
  (void)n_avg;

  std::memcpy(forces, forces_acc.contiguous().data_ptr<double>(),
              nAtoms * 3 * sizeof(double));

  fpeh.restore_fpe();
}

// --- MetatomicPotential::computeNeighbors (helper) ---

metatensor_torch::TensorBlock MetatomicPotential::computeNeighbors(
    metatomic_torch::NeighborListOptions request, long nAtoms,
    const double *positions, const double *box, const bool periodic[3]) {

  auto cutoff = request->engine_cutoff(m_metatomic_opts.length_unit);

  VesinOptions options{}; // zero-initialize all fields (incl. algorithm=0=Auto)
  options.cutoff = cutoff;
  options.full = request->full_list();
  options.return_shifts = true;
  options.return_distances = false; // we don't need distances
  options.return_vectors = true;    // metatomic uses vectors for autograd

  VesinNeighborList *vesin_neighbor_list = new VesinNeighborList();

  VesinDevice cpu{VesinCPU, 0};
  const char *error_message = nullptr;
  int status = vesin_neighbors(reinterpret_cast<const double (*)[3]>(positions),
                               static_cast<size_t>(nAtoms),
                               reinterpret_cast<const double (*)[3]>(box),
                               const_cast<bool *>(periodic), cpu, options,
                               vesin_neighbor_list, &error_message);

  if (status != EXIT_SUCCESS) {
    std::string err_str = "vesin_neighbors failed";
    if (error_message != nullptr) {
      err_str += ": " + std::string(error_message);
    }
    delete vesin_neighbor_list;
    throw std::runtime_error(err_str);
  }

  // Convert from vesin to metatomic format
  auto n_pairs = static_cast<int64_t>(vesin_neighbor_list->length);
  auto labels_options_cpu =
      torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU);

  auto pair_samples_values = torch::empty({n_pairs, 5}, labels_options_cpu);
  auto pair_samples_values_ptr = pair_samples_values.accessor<int32_t, 2>();
  for (int64_t i = 0; i < n_pairs; i++) {
    pair_samples_values_ptr[i][0] =
        static_cast<int32_t>(vesin_neighbor_list->pairs[i][0]);
    pair_samples_values_ptr[i][1] =
        static_cast<int32_t>(vesin_neighbor_list->pairs[i][1]);
    pair_samples_values_ptr[i][2] = vesin_neighbor_list->shifts[i][0];
    pair_samples_values_ptr[i][3] = vesin_neighbor_list->shifts[i][1];
    pair_samples_values_ptr[i][4] = vesin_neighbor_list->shifts[i][2];
  }

  // Custom deleter to free vesin's memory when the torch tensor is destroyed
  auto deleter = [=](void *) {
    vesin_free(vesin_neighbor_list);
    delete vesin_neighbor_list;
  };

  auto pair_vectors = torch::from_blob(
      vesin_neighbor_list->vectors, {n_pairs, 3, 1}, deleter,
      torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU));

  auto neighbor_samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
      std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a",
                               "cell_shift_b", "cell_shift_c"},
      pair_samples_values.to(this->device_));

  auto labels_options_dev =
      torch::TensorOptions().dtype(torch::kInt32).device(this->device_);
  auto neighbor_component =
      torch::make_intrusive<metatensor_torch::LabelsHolder>(
          "xyz", torch::tensor({0, 1, 2}, labels_options_dev).reshape({3, 1}));
  auto neighbor_properties =
      torch::make_intrusive<metatensor_torch::LabelsHolder>(
          "distance", torch::zeros({1, 1}, labels_options_dev));

  return torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
      pair_vectors.to(this->dtype_).to(this->device_), neighbor_samples,
      std::vector<metatensor_torch::Labels>{neighbor_component},
      neighbor_properties);
}

// --- MetatomicPotential::forceBatch ---
// Processes N systems sequentially through the same model instance.
// Numerically identical to N individual force() calls. The single-instance
// design avoids JIT profiling divergence and N model copies. True batched
// model.forward({sys0..sysN}) is a future optimization (see #if 0 block below).

void MetatomicPotential::forceBatch(long nSystems, long nAtoms,
                                    const double *const *positions,
                                    const int *const *atomicNrs,
                                    double *const *forces, double *energies,
                                    double *variances,
                                    const double *const *boxes) {
  // Sequential evaluation through force() -- numerically identical to
  // N individual computePotential() calls. The mutex inside force()
  // serializes, and all calls share the same model instance + JIT state.
  for (long s = 0; s < nSystems; s++) {
    double var = 0;
    force(nAtoms, positions[s], atomicNrs[s], forces[s], &energies[s], &var,
          boxes[s]);
    if (variances)
      variances[s] = var;
    forceCallCounter++;
    PotRegistry::get().on_force_call(ptype);
  }
}

// --- True batched forward (future optimization) ---
// model.forward({sys0..sysN}) verified identical in Python.
// C++ energy extraction needs work to match single-system path exactly.
#if 0
void MetatomicPotential::forceBatchNative(long nSystems, long nAtoms,
                                     const double *const *positions,
                                     const int *const *atomicNrs,
                                     double *const *forces, double *energies,
                                     double *variances,
                                     const double *const *boxes) {
  std::lock_guard<std::mutex> lock(inference_mutex_);

  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();

  auto f64_options =
      torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);

  std::vector<metatomic_torch::System> systems;
  std::vector<torch::Tensor> pos_tensors;
  systems.reserve(static_cast<size_t>(nSystems));
  pos_tensors.reserve(static_cast<size_t>(nSystems));

  for (long s = 0; s < nSystems; s++) {
    auto torch_positions =
        torch::from_blob(const_cast<double *>(positions[s]), {nAtoms, 3},
                         f64_options)
            .to(this->dtype_)
            .to(this->device_)
            .set_requires_grad(true);
    pos_tensors.push_back(torch_positions);

    auto torch_cell =
        torch::from_blob(const_cast<double *>(boxes[s]), {3, 3}, f64_options)
            .to(this->dtype_)
            .to(this->device_);

    auto cell_norms = torch::norm(torch_cell, 2, /*dim=*/1);
    auto torch_pbc = cell_norms.abs() > 1e-9;
    bool periodic[3] = {torch_pbc[0].item<bool>(), torch_pbc[1].item<bool>(),
                        torch_pbc[2].item<bool>()};

    if (!atomicNrs[s]) {
      throw std::runtime_error(
          "[MetatomicPotential] `atomicNrs` must be provided.");
    }
    std::vector<int32_t> types_vec(atomicNrs[s], atomicNrs[s] + nAtoms);
    auto atomic_types =
        torch::tensor(types_vec, torch::TensorOptions().dtype(torch::kInt32))
            .to(this->device_);

    auto system = torch::make_intrusive<metatomic_torch::SystemHolder>(
        atomic_types, torch_positions, torch_cell, torch_pbc);

    // Compute and register neighbor lists for this system
    for (const auto &request : this->nl_requests_) {
      auto neighbors = this->computeNeighbors(request, nAtoms, positions[s],
                                               boxes[s], periodic);
      metatomic_torch::register_autograd_neighbors(system, neighbors,
                                                    this->check_consistency_);
      system->add_neighbor_list(request, neighbors);
    }

    systems.push_back(system);
  }

  // Single batched forward pass
  metatensor_torch::TensorMap output_map;
  try {
    auto ivalue_output = this->model_.forward({
        systems,
        evaluations_options_,
        this->check_consistency_,
    });
    auto dict_output = ivalue_output.toGenericDict();
    output_map = dict_output.at(this->energy_key_)
                     .toCustomClass<metatensor_torch::TensorMapHolder>();
  } catch (const std::exception &e) {
    QUILL_LOG_ERROR(m_log,
                    "[MetatomicPotential] Batched model evaluation failed: {}",
                    e.what());
    throw;
  }

  // Extract per-system energies from the output TensorMap.
  // For per-atom output: samples have ["system", "atom"] dimensions.
  // For system-level output: samples have ["system"] dimension.
  // In both cases, sum over all non-system dimensions to get per-system energy.
  auto energy_block =
      metatensor_torch::TensorMapHolder::block_by_id(output_map, 0);
  auto energy_values = energy_block->values();
  auto samples = energy_block->samples();

  // Check if this is per-atom or per-system output
  bool per_atom = samples->size() > 0 && samples->names().size() > 1;

  if (per_atom) {
    // Per-atom output: sum energies by system index
    auto system_col = samples->column("system").to(torch::kCPU);
    auto flat_energies = energy_values.reshape({-1}).to(torch::kCPU);

    // Sum per-system energies for output
    for (long s = 0; s < nSystems; s++) {
      auto mask = (system_col == s);
      energies[s] =
          flat_energies.index({mask}).sum().to(torch::kFloat64).item<double>();
    }

    // Backward: sum ALL energies, single backward call.
    // Each system's positions.grad() gets only its own contribution
    // because energy_i depends only on positions_i.
    energy_values.sum().backward();
  } else {
    // System-level output: values shape is (nSystems, 1) or similar
    auto cpu_energies = energy_values.to(torch::kCPU).to(torch::kFloat64);
    for (long s = 0; s < nSystems; s++) {
      energies[s] = cpu_energies[s].sum().item<double>();
    }
    energy_values.backward(torch::ones_like(energy_values));
  }

  // Extract per-system forces from position gradients
  for (long s = 0; s < nSystems; s++) {
    auto positions_grad = pos_tensors[s].grad();
    auto forces_tensor =
        -positions_grad.to(torch::kCPU).to(torch::kFloat64);
    std::memcpy(forces[s], forces_tensor.contiguous().data_ptr<double>(),
                nAtoms * 3 * sizeof(double));
  }

  // Variances: not yet supported in batched path
  if (variances) {
    for (long s = 0; s < nSystems; s++) {
      variances[s] = 0.0;
    }
  }

  fpeh.restore_fpe();
}
#endif
