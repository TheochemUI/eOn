/*
** This file is part of eOn.
*/
#include "eon/potentials/Rgpot/RgpotPot.h"
#include "eon/Parameters.h"
#include "eon/PotRegistry.h"
#include "eon/potentials/Rgpot/RGPotEngine.h"

#include <algorithm>
#include <cctype>
#include <chrono>
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
  opt.task_name = o.task_name;

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
      << " (dlopen: "
         "libnwchemc/libcpmdc/librgpot_metatomic_engine/librgpot_xtb_engine)"
      << std::endl;
  // Finalize is registered first. The grouped-exit handler is next, and
  // the stop handler is last, so exit runs stop, then Finalize, then _Exit.
  impl_->finalizeMpiAtExit();
  impl_->armGroupedExit();
  releaseWorkersAtExit();
  if (!driver_)
    serveWorker();
}

RgpotPot::~RgpotPot() {
  const bool grouped = impl_ && impl_->calculatorWorld() > 1;
  stopAndDrop();
  // ~CPMDPot dlcloses libcpmdc. On several ranks that runs Fortran
  // destructors while another rank may already be in MPI_Finalize.
  // The process reclaims the engine at _Exit.
  if (grouped)
    (void)impl_.release();
}

bool RgpotPot::engineAvailable() const { return impl_ && impl_->available(); }

namespace {
// Request header broadcast from the driver: kind, atoms, systems, and the
// group a single request runs on.
// kDown is the second broadcast: every rank has dropped CPMD, and
// workers may enter MPI_Finalize.
enum : std::int64_t { kStop = 0, kSingle = 1, kBatch = 2, kDown = 3 };

// The driver's live potential, so an exit that skips the destructor still
// releases the workers before MPI_Finalize.
RgpotPot *g_driver = nullptr;
} // namespace

void RgpotPot::sendStop() {
  if (!impl_ || !driver_ || stopped_ || impl_->calculatorWorld() <= 1)
    return;
  std::int64_t hdr[4] = {kStop, 0, 0, 0};
  impl_->broadcastFromDriver(hdr, sizeof(hdr));
  stopped_ = true;
  if (g_driver == this)
    g_driver = nullptr;
}

void RgpotPot::exchangeGroupUse() {
  const int groups = impl_->calculatorGroups();
  use_.busy.assign(static_cast<size_t>(groups), 0.0);
  use_.systems.assign(static_cast<size_t>(groups), 0.0);
  for (int g = 0; g < groups; g++) {
    double tally[3] = {busy_, systemsDone_, 0.0};
    double unused = 0.0;
    std::string error;
    (void)impl_->shareResult(g, 1, tally, &unused, true, error);
    use_.busy[static_cast<size_t>(g)] = tally[0];
    use_.systems[static_cast<size_t>(g)] = tally[1];
  }
}

void RgpotPot::stopAndDrop() {
  if (!impl_)
    return;
  const bool grouped = driver_ && impl_->calculatorWorld() > 1;
  if (grouped && !stopped_) {
    sendStop();
    try {
      exchangeGroupUse();
      std::cout << "RgpotPot: " << use_.table() << std::flush;
    } catch (const std::exception &ex) {
      std::cerr << "RgpotPot: no calculator-group summary: " << ex.what()
                << std::endl;
    }
  }
  if (!dropped_) {
    impl_->shutdownModule();
    dropped_ = true;
  }
  if (impl_->calculatorWorld() > 1)
    impl_->armGroupedExit();
  if (grouped && !acked_) {
    std::int64_t hdr[4] = {kDown, 0, 0, 0};
    impl_->broadcastFromDriver(hdr, sizeof(hdr));
    acked_ = true;
  }
}

void RgpotPot::releaseWorkersAtExit() {
  if (g_driver || !impl_ || impl_->calculatorWorld() <= 1)
    return;
  // Runs before the MPI_Finalize handler registered in the constructor.
  // After a failed engine call the workers can sit in a collective that
  // the stop broadcast never meets; the grouped exit handler aborts the
  // world instead.
  g_driver = this;
  std::atexit([] {
    if (g_driver && !RGPotEngine::mpiAbortRequested())
      g_driver->stopAndDrop();
  });
}

namespace {
// Runs one evaluation and keeps its error instead of throwing, so the
// rank still joins every broadcast that follows.
bool try_force(const RGPotEngine &engine, long N, const double *R,
               const int *atomicNrs, double *F, double *U, const double *box,
               std::string &error) {
  try {
    engine.force(N, R, atomicNrs, F, U, box);
    return true;
  } catch (const std::exception &ex) {
    error = ex.what();
    return false;
  }
}

// The orbital key of system j of a batch: its owner hint (a NEB image, a
// ring bead) when it has one, else its batch position, kept apart from
// the hints. A single evaluation has its own key. cpmdc keeps converged
// orbitals per key, so a calculator that evaluates several systems in
// turn starts each SCF from that system's own previous orbitals.
constexpr std::int64_t kSingleKey = -1;
std::int64_t orbital_key(const std::int64_t *owners, long j) {
  if (owners && owners[j] >= 0)
    return owners[j];
  return -2 - static_cast<std::int64_t>(j);
}

[[noreturn]] void raise_failure(long system, int owner,
                                const std::string &error) {
  std::string msg = "RGPOT: calculator " + std::to_string(owner) +
                    " failed on system " + std::to_string(system);
  if (!error.empty())
    msg += ": " + error;
  throw std::runtime_error(msg);
}
} // namespace

bool RgpotPot::evaluate(long N, const double *R, const int *atomicNrs,
                        double *F, double *U, const double *box,
                        std::string &error) {
  const auto t0 = std::chrono::steady_clock::now();
  const bool ok = try_force(*impl_, N, R, atomicNrs, F, U, box, error);
  busy_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
               .count();
  systemsDone_ += 1.0;
  return ok;
}

void RgpotPot::computeSingle(long N, const double *R, const int *atomicNrs,
                             double *F, double *U, const double *box,
                             int group) {
  // Every calculator evaluates the structure, because the engine agrees
  // on errors across MPI_COMM_WORLD after each call and a calculator that
  // skipped the call would never join that agreement. The first rank of
  // `group`, the one whose stored orbitals are nearest, sends its result
  // to every rank.
  std::string error;
  // Only that group's call counts toward the calculator usage report.
  impl_->selectOrbitals(kSingleKey);
  const bool ok = impl_->calculatorIndex() == group
                      ? evaluate(N, R, atomicNrs, F, U, box, error)
                      : try_force(*impl_, N, R, atomicNrs, F, U, box, error);
  if (!impl_->shareResult(group, N, F, U, ok, error))
    raise_failure(0, group, error);
}

void RgpotPot::computeBatch(long nSystems, long nAtoms,
                            const double *const *positions,
                            const int *const *atomicNrs, double *const *forces,
                            double *energies, const double *const *boxes,
                            const std::int64_t *owners,
                            const std::int64_t *route) {
  const int groups = impl_->calculatorGroups();
  const int mine = impl_->calculatorIndex();
  auto ownerOf = [&](long j) { return static_cast<int>(route[j]); };
  if (impl_->calculatorWorld() <= 1) {
    for (long j = 0; j < nSystems; j++) {
      impl_->selectOrbitals(orbital_key(owners, j));
      impl_->force(nAtoms, positions[j], atomicNrs[j], forces[j], &energies[j],
                   boxes[j]);
    }
    return;
  }
  std::vector<char> ok(static_cast<size_t>(nSystems), 1);
  std::vector<std::string> errors(static_cast<size_t>(nSystems));
  std::vector<long> owned(static_cast<size_t>(groups), 0);
  long last = -1;
  for (long j = 0; j < nSystems; j++) {
    owned[static_cast<size_t>(ownerOf(j))]++;
    if (ownerOf(j) == mine) {
      last = j;
      impl_->selectOrbitals(orbital_key(owners, j));
      ok[static_cast<size_t>(j)] =
          evaluate(nAtoms, positions[j], atomicNrs[j], forces[j], &energies[j],
                   boxes[j], errors[static_cast<size_t>(j)]);
    }
  }
  // Every calculator makes the same number of engine calls: the engine
  // agrees on errors across MPI_COMM_WORLD after each one. A calculator
  // that owns fewer systems repeats one into scratch buffers, its last
  // own system or system 0 when it owns none. The repeat runs under that
  // system's key, so a repeat of the last own system is the session's
  // stored result and costs no SCF.
  const long most = *std::max_element(owned.begin(), owned.end());
  const long pad = most - owned[static_cast<size_t>(mine)];
  if (pad > 0) {
    const long j = last >= 0 ? last : 0;
    impl_->selectOrbitals(orbital_key(owners, j));
    std::vector<double> scratchF(static_cast<size_t>(3 * nAtoms));
    double scratchU = 0.0;
    std::string scratchError;
    for (long k = 0; k < pad; k++)
      try_force(*impl_, nAtoms, positions[j], atomicNrs[j], scratchF.data(),
                &scratchU, boxes[j], scratchError);
  }
  // Every share runs before any rank raises, so all ranks leave together.
  long failed = -1;
  for (long j = 0; j < nSystems; j++) {
    const int owner = ownerOf(j);
    if (!impl_->shareResult(owner, nAtoms, forces[j], &energies[j],
                            ok[static_cast<size_t>(j)] != 0,
                            errors[static_cast<size_t>(j)]) &&
        failed < 0)
      failed = j;
  }
  if (failed >= 0)
    raise_failure(failed, ownerOf(failed), errors[static_cast<size_t>(failed)]);
}

void RgpotPot::serveWorker() {
  for (;;) {
    std::int64_t hdr[4] = {kStop, 0, 0, 0};
    impl_->broadcastFromDriver(hdr, sizeof(hdr));
    if (hdr[0] == kStop) {
      // Drop CPMD, wait until the driver has dropped it too, then exit 0.
      // MPI_Finalize is collective and runs from the exit handler.
      try {
        exchangeGroupUse();
      } catch (const std::exception &) {
      }
      if (!dropped_) {
        impl_->shutdownModule();
        dropped_ = true;
      }
      std::int64_t ack[4] = {kDown, 0, 0, 0};
      impl_->broadcastFromDriver(ack, sizeof(ack));
      impl_->armGroupedExit();
      std::exit(0);
    }
    const long n = static_cast<long>(hdr[1]);
    const long m = hdr[0] == kBatch ? static_cast<long>(hdr[2]) : 1;
    std::vector<double> R(static_cast<size_t>(3 * n * m));
    std::vector<int> Z(static_cast<size_t>(n * m));
    std::vector<double> box(static_cast<size_t>(9 * m));
    impl_->broadcastFromDriver(R.data(), R.size() * sizeof(double));
    impl_->broadcastFromDriver(Z.data(), Z.size() * sizeof(int));
    impl_->broadcastFromDriver(box.data(), box.size() * sizeof(double));
    std::vector<std::int64_t> ids(static_cast<size_t>(m), -1);
    std::vector<std::int64_t> route(static_cast<size_t>(m), 0);
    if (hdr[0] == kBatch) {
      impl_->broadcastFromDriver(ids.data(), ids.size() * sizeof(std::int64_t));
      impl_->broadcastFromDriver(route.data(),
                                 route.size() * sizeof(std::int64_t));
    }
    std::vector<double> F(R.size()), U(static_cast<size_t>(m));
    // A failed request raised on the driver too; the driver decides
    // whether the job goes on, so a worker keeps serving.
    if (hdr[0] == kSingle) {
      try {
        computeSingle(n, R.data(), Z.data(), F.data(), U.data(), box.data(),
                      static_cast<int>(hdr[3]));
      } catch (const std::runtime_error &) {
      }
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
    try {
      computeBatch(m, n, pos.data(), nrs.data(), frc.data(), U.data(),
                   bx.data(), ids.data(), route.data());
    } catch (const std::runtime_error &) {
    }
  }
}

void RgpotPot::force(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, double *variance, const double *box) {
  if (variance)
    *variance = 0.0;
  if (impl_->calculatorWorld() <= 1) {
    impl_->selectOrbitals(kSingleKey);
    impl_->force(N, R, atomicNrs, F, U, box);
    return;
  }
  const auto t0 = std::chrono::steady_clock::now();
  if (!schedule_)
    schedule_ =
        std::make_unique<eonc::GroupSchedule>(impl_->calculatorGroups());
  // The group whose last geometry is nearest holds the warmest orbitals.
  const int group = schedule_->nearest(R, 3 * N);
  schedule_->settle(kSingleKey, group);
  schedule_->record(kSingleKey, R, 3 * N);
  std::int64_t hdr[4] = {kSingle, N, 1, group};
  impl_->broadcastFromDriver(hdr, sizeof(hdr));
  impl_->broadcastFromDriver(const_cast<double *>(R), 3 * N * sizeof(double));
  impl_->broadcastFromDriver(const_cast<int *>(atomicNrs), N * sizeof(int));
  impl_->broadcastFromDriver(const_cast<double *>(box), 9 * sizeof(double));
  // Before the evaluation, so a failure that ends the job still stops
  // the workers.
  releaseWorkersAtExit();
  computeSingle(N, R, atomicNrs, F, U, box, group);
  const double dt =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
          .count();
  use_.wall += dt;
  use_.singleWall += dt;
  use_.singles++;
}

bool RgpotPot::supportsBatchEvaluation() const noexcept {
  // One calculator gains from a batch too when the engine keeps orbitals
  // per key: the batch carries each image's or bead's identity, a single
  // force() does not.
  return impl_ &&
         (impl_->calculatorGroups() > 1 || impl_->keepsOrbitalsPerKey());
}

void RgpotPot::forceBatch(long nSystems, long nAtoms,
                          const double *const *positions,
                          const int *const *atomicNrs, double *const *forces,
                          double *energies, double *variances,
                          const double *const *boxes) {
  forceBatchOwned(nSystems, nAtoms, positions, atomicNrs, forces, energies,
                  variances, boxes, nullptr);
}

void RgpotPot::forceBatchOwned(long nSystems, long nAtoms,
                               const double *const *positions,
                               const int *const *atomicNrs,
                               double *const *forces, double *energies,
                               double *variances, const double *const *boxes,
                               const long *owners) {
  std::vector<std::int64_t> ids(static_cast<size_t>(nSystems), -1);
  if (owners) {
    for (long j = 0; j < nSystems; j++)
      ids[static_cast<size_t>(j)] = owners[j];
  }
  const bool grouped = impl_->calculatorWorld() > 1;
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<std::int64_t> route(static_cast<size_t>(nSystems), 0);
  if (grouped) {
    if (!schedule_)
      schedule_ =
          std::make_unique<eonc::GroupSchedule>(impl_->calculatorGroups());
    std::vector<std::int64_t> keys(static_cast<size_t>(nSystems));
    for (long j = 0; j < nSystems; j++)
      keys[static_cast<size_t>(j)] = orbital_key(ids.data(), j);
    const std::vector<int> groupOf = schedule_->assign(keys);
    for (long j = 0; j < nSystems; j++) {
      route[static_cast<size_t>(j)] = groupOf[static_cast<size_t>(j)];
      schedule_->record(keys[static_cast<size_t>(j)], positions[j], 3 * nAtoms);
    }
    std::int64_t hdr[4] = {kBatch, nAtoms, nSystems, 0};
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
    impl_->broadcastFromDriver(ids.data(), ids.size() * sizeof(std::int64_t));
    impl_->broadcastFromDriver(route.data(),
                               route.size() * sizeof(std::int64_t));
  }
  releaseWorkersAtExit();
  computeBatch(nSystems, nAtoms, positions, atomicNrs, forces, energies, boxes,
               ids.data(), route.data());
  if (grouped) {
    use_.wall +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
            .count();
    use_.batches++;
  }
  for (long j = 0; j < nSystems; j++) {
    if (variances)
      variances[j] = 0.0;
    forceCallCounter++;
    eonc::PotRegistry::get().on_force_call(ptype);
  }
}
