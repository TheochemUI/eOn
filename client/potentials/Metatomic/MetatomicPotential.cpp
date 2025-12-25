#include "MetatomicPotential.h"
#include "../../Parameters.h"
#include "vesin.h"

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>

// TODO(rg): handle variants
// https://docs.metatensor.org/metatomic/latest/outputs/variants.html#output-variants

using namespace std::string_literals;
MetatomicPotential::MetatomicPotential(std::shared_ptr<Parameters> params)
    : Potential(PotType::METATOMIC, params),
      model_(torch::jit::Module()),
      device_type_(c10::DeviceType::CPU),
      device_(torch::Device(device_type_)) {

  m_params = params;
  m_log->info("[MetatomicPotential] Initializing...");

  // 1. Load the model from the path specified in parameters
  torch::optional<std::string> extensions_directory = torch::nullopt;
  if (!m_params->metatomic_options.extensions_directory.empty()) {
    extensions_directory = m_params->metatomic_options.extensions_directory;
  }

  try {
    m_log->info("[MetatomicPotential] Loading model from '{}'",
                m_params->metatomic_options.model_path);
    this->model_ = metatomic_torch::load_atomistic_model(
        m_params->metatomic_options.model_path, extensions_directory);
  } catch (const std::exception &e) {
    m_log->error("[MetatomicPotential] Failed to load model: {}", e.what());
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
  // TODO(rg):: Eventually switch to upstream selection
  // https://github.com/metatensor/metatomic/issues/21
  torch::optional<std::string> desired = torch::nullopt;
  if (!m_params->metatomic_options.device.empty()) {
    desired = m_params->metatomic_options.device;
  }
  device_type_ = metatomic_torch::pick_device(
      this->capabilities_->supported_devices, desired);

  device_ = torch::Device(device_type_);
  m_log->info("[MetatomicPotential] Using device: {}", device_.str());

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
  m_log->info("[MetatomicPotential] Using dtype: {}",
              this->capabilities_->dtype().c_str());

  // 5. Set up evaluation options to request total energy
  this->evaluations_options_ =
      torch::make_intrusive<metatomic_torch::ModelEvaluationOptionsHolder>();
  evaluations_options_->set_length_unit(
      m_params->metatomic_options.length_unit);

  auto outputs = this->capabilities_->outputs();
  if (!outputs.contains("energy")) {
    throw std::runtime_error("Metatomic model must provide an 'energy' output "
                             "to be used as a potential.");
  }
  auto energy_output_info = outputs.at("energy");

  auto requested_output =
      torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
  requested_output->per_atom = false; // We want a single total energy value
  requested_output->explicit_gradients = {}; // Use autograd for forces

  evaluations_options_->outputs.insert("energy", requested_output);

  // 5b. Optionally request energy uncertainty (per-atom) if the model offers it
  // and the user didn't disable it. We read a threshold from parameters to know
  // whether to warn (non-positive threshold disables checking).
  // Default: if user provided a positive threshold, use it; otherwise keep it
  // disabled (-1).
  if (m_params->metatomic_options.uncertainty_threshold > 0) {
    this->uncertainty_threshold_ =
        m_params->metatomic_options.uncertainty_threshold;
  } else {
    this->uncertainty_threshold_ = -1.0; // disabled
  }

  auto uncertainty_it = outputs.find("energy_uncertainty");
  if (uncertainty_it != outputs.end() && this->uncertainty_threshold_ > 0) {
    // If the model has per-atom uncertainty output, request it for evaluation.
    // Prefer only per-atom uncertainties for per-atom checking.
    try {
      // The capabilities entry may provide whether the output is per-atom.
      // We attempt to detect that; if not per-atom, we do not request it here.
      auto uncertainty_info = uncertainty_it->value();
      bool per_atom = false;
      try {
        per_atom = uncertainty_info->per_atom;
      } catch (...) {
        per_atom = false;
      }

      if (per_atom) {
        auto requested_uncertainty =
            torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
        requested_uncertainty->per_atom = true;
        requested_uncertainty->explicit_gradients = {};
        requested_uncertainty->set_quantity("energy");
        // leave unit unset or let metatomic handle unit conversion if available
        evaluations_options_->outputs.insert("energy_uncertainty",
                                             requested_uncertainty);
        m_log->info("[MetatomicPotential] Requested per-atom "
                    "'energy_uncertainty' from model (threshold = {})",
                    this->uncertainty_threshold_);
      } else {
        m_log->debug("[MetatomicPotential] Model provides 'energy_uncertainty' "
                     "but not per-atom; skipping per-atom uncertainty checks.");
      }
    } catch (...) {
      // If introspection fails, we skip requesting uncertainty
      m_log->debug("[MetatomicPotential] Could not determine "
                   "'energy_uncertainty' per-atom capability; skipping.");
    }
  }

  this->check_consistency_ = m_params->metatomic_options.check_consistency;
  m_log->info("[MetatomicPotential] Initialization complete.");
}

// --- MetatomicPotential::force ---

void MetatomicPotential::force(long nAtoms, const double *positions,
                               const int *atomicNrs, double *forces,
                               double *energy, double *variance,
                               const double *box) {
  // 1. Convert input arrays to torch::Tensors
  auto f64_options =
      torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);

  auto torch_positions = torch::from_blob(const_cast<double *>(positions),
                                          {nAtoms, 3}, f64_options)
                             .to(this->dtype_)
                             .to(this->device_)
                             .set_requires_grad(true);

  auto torch_cell =
      torch::from_blob(const_cast<double *>(box), {3, 3}, f64_options)
          .to(this->dtype_)
          .to(this->device_);

  auto cell_norms = torch::norm(torch_cell, 2, /*dim=*/1);
  auto torch_pbc = cell_norms.abs() > 1e-9;
  bool periodic = torch::all(torch_pbc).item<bool>();

  bool types_changed = false;
  if (static_cast<size_t>(nAtoms) != last_atomic_nrs_.size()) {
    types_changed = true;
  } else if (nAtoms > 0 && std::memcmp(atomicNrs, last_atomic_nrs_.data(),
                                       nAtoms * sizeof(int)) != 0) {
    types_changed = true;
  }

  if (types_changed) {
    m_log->trace("[MetatomicPotential] Atomic numbers changed, re-creating "
                 "types tensor.");
    if (!atomicNrs) {
      throw std::runtime_error(
          "[MetatomicPotential] `atomicNrs` must be provided.");
    }
    std::vector<int32_t> types_vec(atomicNrs, atomicNrs + nAtoms);
    // XXX(rg): Reordering of labels might take place for non-conservative
    // forces / per atom energies
    this->atomic_types_ =
        torch::tensor(types_vec, torch::TensorOptions().dtype(torch::kInt32))
            .to(this->device_);

    // Update the cache
    last_atomic_nrs_.assign(atomicNrs, atomicNrs + nAtoms);
  }

  // 2. Create the metatomic::System object
  auto system = torch::make_intrusive<metatomic_torch::SystemHolder>(
      this->atomic_types_, torch_positions, torch_cell, torch_pbc);

  // 3. Compute and add neighbor lists to the system
  for (const auto &request : this->nl_requests_) {
    auto neighbors =
        this->computeNeighbors(request, nAtoms, positions, box, periodic);
    metatomic_torch::register_autograd_neighbors(system, neighbors,
                                                 this->check_consistency_);
    system->add_neighbor_list(request, neighbors);
  }

  // 4. Execute the model
  metatensor_torch::TensorMap output_map;
  try {
    auto ivalue_output = this->model_.forward({
        std::vector<metatomic_torch::System>{system},
        evaluations_options_,
        this->check_consistency_,
    });
    auto dict_output = ivalue_output.toGenericDict();
    output_map = dict_output.at("energy")
                     .toCustomClass<metatensor_torch::TensorMapHolder>();

    // If we requested per-atom uncertainty and provided a positive threshold,
    // check the returned uncertainty values and log/warn if any atoms exceed
    // the threshold.
    if (this->uncertainty_threshold_ > 0) {
      try {
        // Check if the dict contains an uncertainty entry
        if (dict_output.contains("energy_uncertainty")) {
          auto uncertainty_map =
              dict_output.at("energy_uncertainty")
                  .toCustomClass<metatensor_torch::TensorMapHolder>();
          auto uncertainty_block =
              metatensor_torch::TensorMapHolder::block_by_id(uncertainty_map,
                                                             0);

          // Expecting a per-atom tensor shaped [n_atoms, 1] (or similar)
          auto uncertainty_values = uncertainty_block->values();
          // Flatten to 1D of per-atom uncertainties
          auto flat_uncertainty =
              uncertainty_values.reshape({-1}).to(torch::kCPU);
          // If variance pointer provided, set it to the mean of per-atom
          // uncertainties
          if (variance != nullptr && flat_uncertainty.numel() > 0) {
            try {
              // TODO(rg): maybe allow the user to set the variance computation
              auto mean_unc = flat_uncertainty.to(torch::kFloat64).mean();
              *variance = mean_unc.item<double>();
            } catch (...) {
              // If mean computation fails, leave variance untouched.
              m_log->debug("[MetatomicPotential] Failed to compute mean "
                           "uncertainty for variance.");
            }
          }

          // Compare with threshold
          auto atoms_above_threshold =
              flat_uncertainty > this->uncertainty_threshold_;
          if (torch::any(atoms_above_threshold).item<bool>()) {
            // Get the sample "atom" column and index it by the boolean mask
            // Samples holder utilities may vary; follow the typical metatensor
            // API
            auto samples = uncertainty_block->samples();
            // samples->column("atom") should return a tensor with atom indices
            auto atom_indices_all = samples->column("atom").to(torch::kCPU);
            auto atom_indices_above =
                atom_indices_all.index({atoms_above_threshold});

            // Build a short human-readable list of offending atom indices
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

            m_log->warn(
                "[MetatomicPotential] The uncertainty on atomic energies for "
                "{} "
                "are larger than the threshold of {}. Be careful when "
                "analyzing "
                "the results, and consider retraining the model to better "
                "describe these configurations.",
                ss.str(), this->uncertainty_threshold_);
          }
        }
      } catch (const std::exception &e) {
        // Don't fail the run for uncertainty-check failures; log and continue.
        m_log->warn(
            "[MetatomicPotential] Failed to check energy_uncertainty: {}",
            e.what());
      }
    }

  } catch (const std::exception &e) {
    m_log->error("[MetatomicPotential] Model evaluation failed: {}", e.what());
    throw;
  }

  // 5. Extract energy and compute gradients (forces)
  auto energy_block =
      metatensor_torch::TensorMapHolder::block_by_id(output_map, 0);
  auto energy_tensor = energy_block->values(); // This should be a [1, 1] tensor

  if (energy_tensor.sizes().vec() != std::vector<int64_t>{1, 1}) {
    throw std::runtime_error(
        "Model did not return a single scalar energy value.");
  }

  // Set the energy value
  *energy = energy_tensor.item<double>();

  // Compute gradients w.r.t positions
  energy_tensor.backward();
  auto positions_grad = system->positions().grad();

  // 6. Copy forces back to the output array
  // Forces are the negative gradient of the potential energy
  auto forces_tensor = -positions_grad.to(torch::kCPU).to(torch::kFloat64);

  std::memcpy(forces, forces_tensor.contiguous().data_ptr<double>(),
              nAtoms * 3 * sizeof(double));
}

// --- MetatomicPotential::computeNeighbors (helper) ---

metatensor_torch::TensorBlock MetatomicPotential::computeNeighbors(
    metatomic_torch::NeighborListOptions request, long nAtoms,
    const double *positions, const double *box, bool periodic) {

  auto cutoff = request->engine_cutoff(m_params->metatomic_options.length_unit);

  VesinOptions options;
  options.cutoff = cutoff;
  options.full = request->full_list();
  options.return_shifts = true;
  options.return_distances = false; // we don't need distances
  options.return_vectors = true;    // metatomic uses vectors for autograd

  VesinNeighborList *vesin_neighbor_list = new VesinNeighborList();

  VesinDevice cpu{VesinCPU, 0};
  const char *error_message = nullptr;
  int status = vesin_neighbors(
      reinterpret_cast<const double (*)[3]>(positions),
      static_cast<size_t>(nAtoms), reinterpret_cast<const double (*)[3]>(box),
      &periodic, cpu, options, vesin_neighbor_list, &error_message);

  if (status != EXIT_SUCCESS) {
    std::string err_str = "vesin_neighbors failed";
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
