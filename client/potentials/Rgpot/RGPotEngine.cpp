// Isolated TU: rgpot only (no eOn Potential.h) — avoids Cap'n Proto Potential
// clash.
// on_exit is a glibc extension. -std=c++20 hides it unless this is set
// before the libc headers.
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "eon/potentials/Rgpot/RGPotEngine.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#if !defined(_WIN32)
#include <dlfcn.h>
#endif
#if defined(__linux__) && defined(EON_RGPOT_MPI)
#include <mpi.h>
#include <stdlib.h>
#endif

#include <capnp/message.h>
#include <capnp/serialize.h>

#include "eon/potentials/Rgpot/MetatomicEngineLoader.h"
#include "eon/potentials/Rgpot/XTBEngineLoader.h"
#include "rgpot/CPMDPot/CPMDPot.hpp"
#include "rgpot/CalculatorGroup.hpp"
#include "rgpot/NWChemPot/NWChemPot.hpp"
#include "rgpot/rpc/Potentials.capnp.h"

using rgpot::types::AtomMatrix;

#if defined(__linux__) && defined(EON_RGPOT_MPI)
extern "C" void eon_rgpot_hard_exit(int status, void *) {
  int inited = 0;
  int finalized = 0;
  MPI_Initialized(&inited);
  MPI_Finalized(&finalized);
  if (inited && !finalized)
    MPI_Finalize();
  std::fflush(nullptr);
  std::_Exit(status);
}
#endif

namespace {

// Serialized Cap'n Proto message (standard segment framing, word-aligned) as
// written by `capnp encode` or messageToFlatArray.
std::vector<::capnp::word> read_params_file(const std::string &path) {
  std::ifstream in(path, std::ios::binary | std::ios::ate);
  if (!in)
    throw std::runtime_error("RGPOT: cannot open params_path: " + path);
  const std::streamsize bytes = in.tellg();
  if (bytes <= 0 || (static_cast<size_t>(bytes) % sizeof(::capnp::word)) != 0)
    throw std::runtime_error(
        "RGPOT: params_path is not a capnp flat message: " + path);
  std::vector<::capnp::word> words(static_cast<size_t>(bytes) /
                                   sizeof(::capnp::word));
  in.seekg(0);
  in.read(reinterpret_cast<char *>(words.data()), bytes);
  if (!in)
    throw std::runtime_error("RGPOT: short read on params_path: " + path);
  return words;
}

std::string to_lower(std::string s) {
  for (char &c : s)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::array<std::array<double, 3>, 3> box_from_row_major(const double *box) {
  std::array<std::array<double, 3>, 3> out{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out[static_cast<size_t>(i)][static_cast<size_t>(j)] = box[i * 3 + j];
  return out;
}

bool looks_like_dft_xc(const std::string &s) {
  if (s.empty())
    return false;
  static const char *k[] = {"b3lyp", "blyp", "pbe",    "pw91",    "bp86",
                            "hcth",  "ft97", "hfexch", "xperpbe", nullptr};
  for (int i = 0; k[i]; ++i) {
    if (s.size() >= std::char_traits<char>::length(k[i]) &&
        s.compare(0, std::char_traits<char>::length(k[i]), k[i]) == 0)
      return true;
  }
  return false;
}

/** Map eOn/native XTB paramset names to RGPOT_XTB_METHOD_* ABI codes. */
int xtb_method_from_paramset(const std::string &paramset) {
  const std::string p = to_lower(paramset);
  if (p == "gfnff" || p == "gfn-ff")
    return RGPOT_XTB_METHOD_GFNFF;
  if (p == "gfn0xtb" || p == "gfn0" || p == "gfn0-xtb")
    return RGPOT_XTB_METHOD_GFN0;
  if (p == "gfn1xtb" || p == "gfn1" || p == "gfn1-xtb")
    return RGPOT_XTB_METHOD_GFN1;
  if (p == "gfn2xtb" || p == "gfn2" || p == "gfn2-xtb" || p.empty())
    return RGPOT_XTB_METHOD_GFN2;
  throw std::runtime_error(
      "RGPOT(xtb): paramset must be GFNFF, GFN0xTB, GFN1xTB, or GFN2xTB "
      "(got '" +
      paramset + "')");
}

// Ranks mpirun assigned. 1 when this process was not started under mpirun,
// so a one-process test does not call MPI_Init.
int mpi_world_hint() {
  const char *s = std::getenv("OMPI_COMM_WORLD_SIZE");
  if (s == nullptr || s[0] == '\0')
    s = std::getenv("PMI_SIZE");
  if (s == nullptr || s[0] == '\0')
    return 1;
  return std::atoi(s);
}

void pin_cpmd_library(const std::string &path) {
#if !defined(_WIN32)
  // A later dlclose must not run Fortran destructors. The mapping stays
  // until process exit, and the grouped exit path does not run them.
  const char *name = path.empty() ? "libcpmdc.so" : path.c_str();
  dlopen(name, RTLD_NOW | RTLD_NOLOAD | RTLD_GLOBAL | RTLD_NODELETE);
#else
  (void)path;
#endif
}

#if defined(__linux__) && defined(EON_RGPOT_MPI)
// Every rank calls this. A rank that failed locally still enters, so nobody
// is left in MPI_Comm_split. The shared message is the lowest failing rank's
// text. Returns false when any rank failed.
bool agree_construction(std::string &message) {
  int inited = 0;
  MPI_Initialized(&inited);
  if (!inited)
    MPI_Init(nullptr, nullptr);
  ::rgpot::finalizeMpiAtExit();
  const int ok = message.empty() ? 1 : 0;
  int all_ok = 0;
  MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (all_ok)
    return true;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int mine = ok ? size : rank;
  int owner = 0;
  MPI_Allreduce(&mine, &owner, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  int len = 0;
  if (rank == owner)
    len = static_cast<int>(
        std::min(message.size(), static_cast<std::size_t>(4095)));
  MPI_Bcast(&len, 1, MPI_INT, owner, MPI_COMM_WORLD);
  std::vector<char> buf(static_cast<std::size_t>(len) + 1, '\0');
  if (rank == owner && len > 0)
    std::memcpy(buf.data(), message.data(), static_cast<std::size_t>(len));
  if (len > 0)
    MPI_Bcast(buf.data(), len, MPI_CHAR, owner, MPI_COMM_WORLD);
  message.assign(buf.data(), static_cast<std::size_t>(len));
  return false;
}
#endif

} // namespace

struct RGPotEngine::Impl {
  enum class Backend { Nwchemc, Cpmdc, Metatomic, Xtb };
  Backend backend{Backend::Nwchemc};
  std::unique_ptr<rgpot::NWChemPot> nwchem;
  std::unique_ptr<rgpot::CPMDPot> cpmd;
  std::unique_ptr<MetatomicEngineLoader> metatomic;
  std::unique_ptr<XTBEngineLoader> xtb;
  // Calculator groups (cpmdc, ranks_per_image > 0). One group when off.
  int groups{1};
  int group{0};
  int world{1};
  bool module_down{false};
};

RGPotEngine::RGPotEngine(const RGPotEngineOptions &opt)
    : impl_(std::make_unique<Impl>()) {
  backend_ = to_lower(opt.backend);
  if (backend_ == "nwchem" || backend_ == "nwchemc" ||
      backend_ == "nwchempot") {
    backend_ = "nwchemc";
    impl_->backend = Impl::Backend::Nwchemc;
    ::capnp::MallocMessageBuilder msg;
    auto params = msg.initRoot<::NWChemParams>();
    params.setBasis(opt.basis);
    params.setTheory(opt.theory);
    params.setScfType(opt.scf_type);
    params.setCharge(opt.charge);
    params.setMultiplicity(opt.multiplicity);
    if (!opt.engine_path.empty())
      params.setEnginePath(opt.engine_path);
    else if (!opt.engine_library.empty())
      params.setEnginePath(opt.engine_library);
    if (!opt.engine_root.empty())
      params.setNwchemRoot(opt.engine_root);
    if (!opt.title.empty())
      params.setTitle(opt.title);
    if (opt.memory_mb > 0)
      params.setMemoryMb(static_cast<uint32_t>(opt.memory_mb));
    if (!opt.scratch_dir.empty())
      params.setScratchDir(opt.scratch_dir);
    // DFT XC as inputBlocks when theory=dft and scfType is an XC label
    std::string block = opt.input_block;
    if (block.empty()) {
      if (const char *env = std::getenv("RGPOT_NWCHEM_INPUT_BLOCK"))
        block = env;
    }
    if (block.empty() && (opt.theory == "dft" || opt.theory == "DFT") &&
        looks_like_dft_xc(opt.scf_type)) {
      block = "dft\n  xc " + opt.scf_type + "\n  mult " +
              std::to_string(opt.multiplicity) + "\nend";
    } else if (block.empty() && looks_like_dft_xc(opt.theory)) {
      block = "dft\n  xc " + opt.theory + "\n  mult " +
              std::to_string(opt.multiplicity) + "\nend";
    }
    if (!block.empty()) {
      auto blocks = params.initInputBlocks(1);
      blocks.set(0, block);
    }
    impl_->nwchem = std::make_unique<rgpot::NWChemPot>(params.asReader());
    if (!impl_->nwchem->available())
      throw std::runtime_error(
          "RGPOT(nwchemc): engine not available (set NWCHEMC_LIBRARY / "
          "RGPOT_NWCHEMC_ENGINE or [RgpotPot] engine_path)");
  } else if (backend_ == "cpmd" || backend_ == "cpmdc" ||
             backend_ == "cpmdpot") {
    backend_ = "cpmdc";
    impl_->backend = Impl::Backend::Cpmdc;
    std::string local_error;
    try {
      ::capnp::MallocMessageBuilder msg;
      ::CPMDParams::Builder params = msg.initRoot<::CPMDParams>();
      if (!opt.params_path.empty()) {
        // The whole CPMDParams message from disk: functional, cutoff, and the
        // inputSections (cell, symmetry, pseudopotential channels, optimizer).
        // The keys below only place the engine and its scratch.
        const auto words = read_params_file(opt.params_path);
        ::capnp::FlatArrayMessageReader reader(
            kj::arrayPtr(words.data(), words.size()));
        msg.setRoot(reader.getRoot<::CPMDParams>());
        params = msg.getRoot<::CPMDParams>();
      } else {
        params.setFunctional(opt.functional);
        params.setCutOffRy(opt.cutoff_ry);
        params.setCharge(opt.charge);
        params.setMultiplicity(opt.multiplicity);
        if (!opt.title.empty())
          params.setTitle(opt.title);
        if (opt.memory_mb > 0)
          params.setMemoryMb(static_cast<uint32_t>(opt.memory_mb));
      }
      if (!opt.engine_path.empty())
        params.setEnginePath(opt.engine_path);
      else if (!opt.engine_library.empty())
        params.setEnginePath(opt.engine_library);
      if (!opt.engine_root.empty())
        params.setCpmdRoot(opt.engine_root);
      if (!opt.scratch_dir.empty())
        params.setScratchDir(opt.scratch_dir);
      if (!opt.permanent_dir.empty())
        params.setPermanentDir(opt.permanent_dir);
      // Raw &SECTION text goes ahead of the sections cpmdc generates, so a
      // periodic &SYSTEM or a full &DFT here replaces the isolated cold deck.
      std::string block = opt.input_block;
      if (block.empty()) {
        if (const char *env = std::getenv("RGPOT_CPMD_INPUT_BLOCK"))
          block = env;
      }
      if (!block.empty()) {
        auto blocks = params.initInputBlocks(1);
        blocks.set(0, block);
      }
      impl_->cpmd = std::make_unique<rgpot::CPMDPot>(params.asReader());
      pin_cpmd_library(opt.engine_path.empty() ? opt.engine_library
                                               : opt.engine_path);
      if (!impl_->cpmd->available())
        throw std::runtime_error(
            "RGPOT(cpmdc): engine not available (set CPMDC_LIBRARY / "
            "RGPOT_CPMDC_ENGINE or [RgpotPot] engine_path)");
    } catch (const std::exception &ex) {
      local_error = ex.what();
    }
    // Several ranks: agree before MPI_Comm_split. One rank that cannot
    // read params_path must not leave the others inside the split.
    if (mpi_world_hint() <= 1 && !local_error.empty())
      throw std::runtime_error(local_error);
#if defined(__linux__) && defined(EON_RGPOT_MPI)
    if (::rgpot::calculatorsUseMpi() && mpi_world_hint() > 1 &&
        !agree_construction(local_error)) {
      armGroupedExit();
      throw std::runtime_error(
          local_error.empty()
              ? std::string(
                    "RGPOT(cpmdc): a rank failed before the calculator split")
              : local_error);
    }
#endif
    if (!local_error.empty())
      throw std::runtime_error(local_error);
    if (::rgpot::calculatorsUseMpi()) {
      // Collective on MPI_COMM_WORLD, before the first force: the engine
      // installs the group communicator ahead of CPMD's mp_start.
      // ranks_per_image = 0 is one calculator on the whole world.
      const rgpot::CalculatorGroup g =
          ::rgpot::bindCalculators(opt.ranks_per_image);
      ::rgpot::finalizeMpiAtExit();
      if (mpi_world_hint() > 1)
        armGroupedExit();
      if (g.index < 0)
        throw std::runtime_error(
            "RGPOT(cpmdc): ranks_per_image=" +
            std::to_string(opt.ranks_per_image) +
            " does not divide the MPI world into calculator groups");
      impl_->groups = ::rgpot::calculatorCount();
      impl_->group = g.index;
      impl_->world = ::rgpot::calculatorWorldSize();
    } else if (opt.ranks_per_image > 0) {
      throw std::runtime_error(
          "RGPOT(cpmdc): ranks_per_image needs rgpot built with MPI "
          "(-Drgpot:with_mpi=enabled)");
    }
  } else if (backend_ == "metatomic" || backend_ == "mta" ||
             backend_ == "metatomicpot") {
    backend_ = "metatomic";
    impl_->backend = Impl::Backend::Metatomic;
    MetatomicEngineOptions mopt;
    mopt.model_path = opt.model_path;
    mopt.device = opt.device;
    mopt.length_unit = opt.length_unit;
    mopt.extensions_directory = opt.extensions_directory;
    mopt.check_consistency = opt.check_consistency;
    mopt.uncertainty_threshold = opt.uncertainty_threshold;
    mopt.engine_path =
        !opt.engine_path.empty() ? opt.engine_path : opt.engine_library;
    mopt.torch_determinism_strict = opt.torch_determinism_strict;
    impl_->metatomic = std::make_unique<MetatomicEngineLoader>(mopt);
    if (!impl_->metatomic->available())
      throw std::runtime_error(
          "RGPOT(metatomic): engine not available (set RGPOT_METATOMIC_ENGINE "
          "or [RgpotPot] engine_path to libmetatomic_engine.so)");
  } else if (backend_ == "xtb" || backend_ == "xtbpot" || backend_ == "gfn" ||
             backend_ == "gfnxtb") {
    backend_ = "xtb";
    impl_->backend = Impl::Backend::Xtb;
    XTBEngineOptions xopt;
    xopt.method = xtb_method_from_paramset(opt.xtb_paramset);
    xopt.accuracy = opt.xtb_accuracy;
    xopt.electronic_temperature = opt.xtb_electronic_temperature;
    xopt.max_iterations = opt.xtb_max_iterations;
    xopt.charge = opt.xtb_charge;
    xopt.uhf = opt.xtb_uhf;
    xopt.engine_path =
        !opt.engine_path.empty() ? opt.engine_path : opt.engine_library;
    impl_->xtb = std::make_unique<XTBEngineLoader>(xopt);
    if (!impl_->xtb->available())
      throw std::runtime_error(
          "RGPOT(xtb): engine not available (set RGPOT_XTB_ENGINE or "
          "[RgpotPot] engine_path to libxtb_engine.so)");
  } else {
    throw std::runtime_error("RGPOT: unknown backend '" + opt.backend +
                             "' (expected nwchemc, cpmdc, metatomic, or xtb)");
  }
}

RGPotEngine::~RGPotEngine() = default;

int RGPotEngine::calculatorGroups() const noexcept {
  return impl_ ? impl_->groups : 1;
}

int RGPotEngine::calculatorIndex() const noexcept {
  return impl_ ? impl_->group : 0;
}

int RGPotEngine::calculatorWorld() const noexcept {
  return impl_ ? impl_->world : 1;
}

int RGPotEngine::worldRank() const noexcept {
  if (!impl_ || impl_->world <= 1)
    return 0;
  const rgpot::CalculatorGroup &g = ::rgpot::thisCalculator();
  return g.index * g.ranks + g.rank_in_group;
}

void RGPotEngine::finalizeMpiAtExit() const {
  if (impl_ && impl_->world > 1)
    ::rgpot::finalizeMpiAtExit();
}

void RGPotEngine::armGroupedExit() const {
#if defined(__linux__) && defined(EON_RGPOT_MPI)
  if (mpi_world_hint() <= 1)
    return;
  static std::once_flag once;
  std::call_once(once, [] { ::on_exit(eon_rgpot_hard_exit, nullptr); });
#else
  (void)this;
#endif
}

void RGPotEngine::shutdownModule() noexcept {
  if (!impl_ || impl_->module_down)
    return;
  impl_->module_down = true;
#if !defined(_WIN32)
  if (auto *fin =
          reinterpret_cast<void (*)()>(dlsym(RTLD_DEFAULT, "cpmdc_finalize")))
    fin();
  pin_cpmd_library({});
#endif
}

void RGPotEngine::broadcastFromDriver(void *data, std::size_t bytes) const {
  if (!impl_ || impl_->world <= 1 || bytes == 0)
    return;
  // World rank 0 is the first rank of calculator 0.
  if (::rgpot::shareFromCalculator(0, data, bytes) == 0)
    throw std::runtime_error("RGPOT: could not broadcast from the driver rank");
}

bool RGPotEngine::shareResult(int owner, long N, double *F, double *U, bool ok,
                              std::string &error) const {
  if (!impl_ || impl_->world <= 1) {
    (void)error;
    return ok;
  }
  // Slot 0 carries the owner's status, so a failed evaluation reaches
  // every rank through the same broadcast as a good one.
  std::vector<double> buf(static_cast<size_t>(3 * N + 2));
  if (owner == impl_->group) {
    buf[0] = ok ? 1.0 : 0.0;
    buf[1] = *U;
    std::copy(F, F + 3 * N, buf.begin() + 2);
  }
  if (::rgpot::shareFromCalculator(owner, buf.data(),
                                   buf.size() * sizeof(double)) == 0)
    throw std::runtime_error("RGPOT: could not share a calculator result");
  *U = buf[1];
  std::copy(buf.begin() + 2, buf.end(), F);
  // The owner's engine text travels with the status. Other ranks hold an
  // empty string until this broadcast.
  std::array<char, 512> msg{};
  if (owner == impl_->group && !error.empty()) {
    const auto n = std::min(error.size(), msg.size() - 1);
    std::memcpy(msg.data(), error.data(), n);
  }
  if (::rgpot::shareFromCalculator(owner, msg.data(), msg.size()) == 0)
    throw std::runtime_error("RGPOT: could not share a calculator error");
  if (msg[0] != '\0')
    error.assign(msg.data());
  return buf[0] == 1.0;
}

bool RGPotEngine::available() const {
  if (!impl_)
    return false;
  if (impl_->backend == Impl::Backend::Nwchemc && impl_->nwchem)
    return impl_->nwchem->available();
  if (impl_->backend == Impl::Backend::Cpmdc && impl_->cpmd)
    return impl_->cpmd->available();
  if (impl_->backend == Impl::Backend::Metatomic && impl_->metatomic)
    return impl_->metatomic->available();
  if (impl_->backend == Impl::Backend::Xtb && impl_->xtb)
    return impl_->xtb->available();
  return false;
}

void RGPotEngine::force(long N, const double *R, const int *atomicNrs,
                        double *F, double *U, const double *box) const {
  // Set on one rank by the multi-rank failure test. The text is the
  // engine error that rank has to deliver to rank 0.
  if (const char *fail = std::getenv("RGPOT_FORCE_FAIL")) {
    if (fail[0] != '\0')
      throw std::runtime_error(fail);
  }
  if (N <= 0)
    throw std::runtime_error("RGPotEngine::force called with N <= 0");

  AtomMatrix positions(static_cast<int>(N), 3);
  for (long i = 0; i < N; ++i) {
    const int ii = static_cast<int>(i);
    positions(ii, 0) = R[3 * i + 0];
    positions(ii, 1) = R[3 * i + 1];
    positions(ii, 2) = R[3 * i + 2];
  }
  std::vector<int> atmtypes(atomicNrs, atomicNrs + N);
  const auto cell = box_from_row_major(box);

  if (impl_->backend == Impl::Backend::Metatomic) {
    impl_->metatomic->force(N, R, atomicNrs, F, U, nullptr, box);
    return;
  }
  if (impl_->backend == Impl::Backend::Xtb) {
    impl_->xtb->force(N, R, atomicNrs, F, U, nullptr, box);
    return;
  }

  // rgpot >= 2.5.0: operator() returns (energy, forces, variance).
  std::tuple<double, AtomMatrix, double> result;
  if (impl_->backend == Impl::Backend::Nwchemc)
    result = (*impl_->nwchem)(positions, atmtypes, cell);
  else
    result = (*impl_->cpmd)(positions, atmtypes, cell);

  *U = std::get<0>(result);
  const auto &forces = std::get<1>(result);
  for (long i = 0; i < N; ++i) {
    const int ii = static_cast<int>(i);
    F[3 * i + 0] = forces(ii, 0);
    F[3 * i + 1] = forces(ii, 1);
    F[3 * i + 2] = forces(ii, 2);
  }
}
