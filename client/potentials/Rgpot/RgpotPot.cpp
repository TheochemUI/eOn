/*
** This file is part of eOn.
*/
#include "eon/potentials/Rgpot/RgpotPot.h"
#include "eon/Parameters.h"
#include "eon/PotRegistry.h"
#include "eon/potentials/Rgpot/RGPotEngine.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>

RgpotPot::RgpotPot(const eonc::Parameters &p)
    : eonc::Potential(eonc::PotType::RGPOT, p) {
  RGPotEngineOptions opt;
  const auto &o = p.rgpot_options();
  opt.backend = o.backend;
  opt.basis = o.basis;
  opt.theory = o.theory;
  opt.scf_type = o.scf_type;
  opt.functional = o.functional;
  opt.cutoff_ry = o.cutoff_ry;
  opt.charge = o.charge;
  opt.multiplicity = o.multiplicity;
  opt.engine_path = o.engine_path;
  opt.engine_library = o.engine_library;
  opt.engine_root = o.engine_root;
  opt.title = o.title;
  opt.memory_mb = o.memory_mb;
  opt.scratch_dir = o.scratch_dir;
  opt.input_block = o.input_block;
  opt.permanent_dir = o.permanent_dir;
  opt.params_path = o.params_path;
  opt.ranks_per_image = o.ranks_per_image;
  opt.model_path = o.model_path;
  opt.device = o.device;
  opt.length_unit = o.length_unit;
  opt.extensions_directory = o.extensions_directory;
  opt.check_consistency = o.check_consistency;
  opt.uncertainty_threshold = o.uncertainty_threshold;
  opt.torch_determinism_strict = o.torch_determinism_strict;
  opt.xtb_paramset = o.xtb_paramset;
  opt.xtb_accuracy = o.xtb_accuracy;
  opt.xtb_electronic_temperature = o.xtb_electronic_temperature;
  opt.xtb_max_iterations = o.xtb_max_iterations;
  opt.xtb_charge = o.xtb_charge;
  opt.xtb_uhf = o.xtb_uhf;

  // Env overrides (CI / benchmarks)
  if (const char *e = std::getenv("RGPOT_BACKEND"))
    opt.backend = e;
  if (const char *e = std::getenv("RGPOT_NWCHEM_BASIS"))
    opt.basis = e;
  if (const char *e = std::getenv("RGPOT_NWCHEM_THEORY"))
    opt.theory = e;
  if (const char *e = std::getenv("RGPOT_NWCHEM_SCF_TYPE"))
    opt.scf_type = e;
  if (const char *e = std::getenv("RGPOT_PARAMS_PATH"))
    opt.params_path = e;
  // Engine-path env overrides are backend-scoped: NWCHEMC_LIBRARY must not
  // leak into a cpmdc configure (CPMDPot resolves CPMDC_LIBRARY itself).
  std::string backend_lc = opt.backend;
  std::ranges::transform(backend_lc, backend_lc.begin(),
                         [](unsigned char c) { return std::tolower(c); });
  if (backend_lc.rfind("nwchem", 0) == 0) {
    if (const char *e = std::getenv("NWCHEMC_LIBRARY"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("RGPOT_NWCHEMC_ENGINE"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("RGPOT_NWCHEM_ENGINE"))
      opt.engine_path = e;
  } else if (backend_lc.rfind("cpmd", 0) == 0) {
    if (const char *e = std::getenv("CPMDC_LIBRARY"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("RGPOT_CPMDC_ENGINE"))
      opt.engine_path = e;
  } else if (backend_lc.rfind("meta", 0) == 0 || backend_lc == "mta") {
    if (const char *e = std::getenv("RGPOT_METATOMIC_ENGINE"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("METATOMIC_ENGINE"))
      opt.engine_path = e;
    if (const char *e = std::getenv("RGPOT_METATOMIC_MODEL"))
      opt.model_path = e;
  } else if (backend_lc == "xtb" || backend_lc == "xtbpot" ||
             backend_lc == "gfn" || backend_lc == "gfnxtb") {
    if (const char *e = std::getenv("RGPOT_XTB_ENGINE"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("XTB_ENGINE"))
      opt.engine_path = e;
    if (opt.xtb_paramset.empty() || opt.xtb_paramset == "GFN2xTB") {
      if (!p.xtb_options().paramset.empty())
        opt.xtb_paramset = p.xtb_options().paramset;
    }
  }

  // Dual-read [Metatomic] when RGPOT backend is metatomic
  if ((backend_lc.rfind("meta", 0) == 0 || backend_lc == "mta") &&
      opt.model_path.empty())
    opt.model_path = p.metatomic_options().model_path;
  if ((backend_lc.rfind("meta", 0) == 0 || backend_lc == "mta") &&
      opt.device == "cpu" && !p.metatomic_options().device.empty())
    opt.device = p.metatomic_options().device;

  impl_ = std::make_unique<RGPotEngine>(opt);
  backend_ = impl_->backend();
  driver_ = impl_->worldRank() == 0;
  std::cout
      << "RgpotPot: in-process rgpot backend=" << backend_
      << " (dlopen: libnwchemc/libcpmdc/libmetatomic_engine/libxtb_engine)"
      << std::endl;
  // Registered first, so MPI_Finalize runs after the driver's stop handler.
  impl_->finalizeMpiAtExit();
  if (!driver_)
    serveWorker();
}

RgpotPot::~RgpotPot() { sendStop(); }

bool RgpotPot::engineAvailable() const { return impl_ && impl_->available(); }

namespace {
// Request header broadcast from the driver: kind, atoms, systems.
enum : std::int64_t { kStop = 0, kSingle = 1, kBatch = 2 };

// The driver's live potential, so an exit that skips the destructor still
// releases the workers before MPI_Finalize.
RgpotPot *g_driver = nullptr;
} // namespace

void RgpotPot::sendStop() {
  if (!impl_ || !driver_ || stopped_ || impl_->calculatorWorld() <= 1)
    return;
  std::int64_t hdr[3] = {kStop, 0, 0};
  impl_->broadcastFromDriver(hdr, sizeof(hdr));
  stopped_ = true;
  if (g_driver == this)
    g_driver = nullptr;
}

void RgpotPot::releaseWorkersAtExit() {
  if (g_driver || impl_->calculatorWorld() <= 1)
    return;
  // Runs before the MPI_Finalize handler registered in the constructor.
  g_driver = this;
  std::atexit([] {
    if (g_driver)
      g_driver->sendStop();
  });
}

void RgpotPot::computeSingle(long N, const double *R, const int *atomicNrs,
                             double *F, double *U, const double *box) {
  // Group 0 computes; its first rank's result reaches every rank.
  if (impl_->calculatorIndex() == 0)
    impl_->force(N, R, atomicNrs, F, U, box);
  impl_->shareResult(0, N, F, U);
}

void RgpotPot::computeBatch(long nSystems, long nAtoms,
                            const double *const *positions,
                            const int *const *atomicNrs, double *const *forces,
                            double *energies, const double *const *boxes) {
  const int groups = impl_->calculatorGroups();
  const int mine = impl_->calculatorIndex();
  for (long j = 0; j < nSystems; j++) {
    if (groups <= 1 || static_cast<int>(j % groups) == mine)
      impl_->force(nAtoms, positions[j], atomicNrs[j], forces[j], &energies[j],
                   boxes[j]);
  }
  if (impl_->calculatorWorld() > 1)
    for (long j = 0; j < nSystems; j++)
      impl_->shareResult(static_cast<int>(j % groups), nAtoms, forces[j],
                         &energies[j]);
}

void RgpotPot::serveWorker() {
  for (;;) {
    std::int64_t hdr[3] = {kStop, 0, 0};
    impl_->broadcastFromDriver(hdr, sizeof(hdr));
    if (hdr[0] == kStop)
      std::exit(0);
    const long n = static_cast<long>(hdr[1]);
    const long m = hdr[0] == kBatch ? static_cast<long>(hdr[2]) : 1;
    std::vector<double> R(static_cast<size_t>(3 * n * m));
    std::vector<int> Z(static_cast<size_t>(n * m));
    std::vector<double> box(static_cast<size_t>(9 * m));
    impl_->broadcastFromDriver(R.data(), R.size() * sizeof(double));
    impl_->broadcastFromDriver(Z.data(), Z.size() * sizeof(int));
    impl_->broadcastFromDriver(box.data(), box.size() * sizeof(double));
    std::vector<double> F(R.size()), U(static_cast<size_t>(m));
    if (hdr[0] == kSingle) {
      computeSingle(n, R.data(), Z.data(), F.data(), U.data(), box.data());
      continue;
    }
    std::vector<const double *> pos(static_cast<size_t>(m)), bx(pos.size());
    std::vector<const int *> nrs(pos.size());
    std::vector<double *> frc(pos.size());
    for (long j = 0; j < m; j++) {
      pos[static_cast<size_t>(j)] = R.data() + 3 * n * j;
      nrs[static_cast<size_t>(j)] = Z.data() + n * j;
      frc[static_cast<size_t>(j)] = F.data() + 3 * n * j;
      bx[static_cast<size_t>(j)] = box.data() + 9 * j;
    }
    computeBatch(m, n, pos.data(), nrs.data(), frc.data(), U.data(), bx.data());
  }
}

void RgpotPot::force(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, double *variance, const double *box) {
  if (variance)
    *variance = 0.0;
  if (impl_->calculatorWorld() <= 1) {
    impl_->force(N, R, atomicNrs, F, U, box);
    return;
  }
  std::int64_t hdr[3] = {kSingle, N, 1};
  impl_->broadcastFromDriver(hdr, sizeof(hdr));
  impl_->broadcastFromDriver(const_cast<double *>(R), 3 * N * sizeof(double));
  impl_->broadcastFromDriver(const_cast<int *>(atomicNrs), N * sizeof(int));
  impl_->broadcastFromDriver(const_cast<double *>(box), 9 * sizeof(double));
  computeSingle(N, R, atomicNrs, F, U, box);
  releaseWorkersAtExit();
}

bool RgpotPot::supportsBatchEvaluation() const noexcept {
  return impl_ && impl_->calculatorGroups() > 1;
}

void RgpotPot::forceBatch(long nSystems, long nAtoms,
                          const double *const *positions,
                          const int *const *atomicNrs, double *const *forces,
                          double *energies, double *variances,
                          const double *const *boxes) {
  if (impl_->calculatorWorld() > 1) {
    std::int64_t hdr[3] = {kBatch, nAtoms, nSystems};
    impl_->broadcastFromDriver(hdr, sizeof(hdr));
    std::vector<double> R(static_cast<size_t>(3 * nAtoms * nSystems));
    std::vector<int> Z(static_cast<size_t>(nAtoms * nSystems));
    std::vector<double> box(static_cast<size_t>(9 * nSystems));
    for (long j = 0; j < nSystems; j++) {
      std::copy(positions[j], positions[j] + 3 * nAtoms,
                R.begin() + 3 * nAtoms * j);
      std::copy(atomicNrs[j], atomicNrs[j] + nAtoms, Z.begin() + nAtoms * j);
      std::copy(boxes[j], boxes[j] + 9, box.begin() + 9 * j);
    }
    impl_->broadcastFromDriver(R.data(), R.size() * sizeof(double));
    impl_->broadcastFromDriver(Z.data(), Z.size() * sizeof(int));
    impl_->broadcastFromDriver(box.data(), box.size() * sizeof(double));
  }
  computeBatch(nSystems, nAtoms, positions, atomicNrs, forces, energies, boxes);
  releaseWorkersAtExit();
  for (long j = 0; j < nSystems; j++) {
    if (variances)
      variances[j] = 0.0;
    forceCallCounter++;
    eonc::PotRegistry::get().on_force_call(ptype);
  }
}
