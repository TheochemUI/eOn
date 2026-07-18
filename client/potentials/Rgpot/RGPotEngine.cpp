// Isolated TU: rgpot only (no eOn Potential.h) — avoids Cap'n Proto Potential
// clash.
#include "RGPotEngine.h"

#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include <capnp/message.h>

#include "rgpot/CPMDPot/CPMDPot.hpp"
#include "rgpot/NWChemPot/NWChemPot.hpp"
#include "MetatomicEngineLoader.h"
#include "XTBEngineLoader.h"
#include "rgpot/rpc/Potentials.capnp.h"

namespace {

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

} // namespace

struct RGPotEngine::Impl {
  enum class Backend { Nwchemc, Cpmdc, Metatomic, Xtb };
  Backend backend{Backend::Nwchemc};
  std::unique_ptr<rgpot::NWChemPot> nwchem;
  std::unique_ptr<rgpot::CPMDPot> cpmd;
  std::unique_ptr<MetatomicEngineLoader> metatomic;
  std::unique_ptr<XTBEngineLoader> xtb;
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
    ::capnp::MallocMessageBuilder msg;
    auto params = msg.initRoot<::CPMDParams>();
    params.setFunctional(opt.functional);
    params.setCutOffRy(opt.cutoff_ry);
    params.setCharge(opt.charge);
    params.setMultiplicity(opt.multiplicity);
    if (!opt.engine_path.empty())
      params.setEnginePath(opt.engine_path);
    else if (!opt.engine_library.empty())
      params.setEnginePath(opt.engine_library);
    if (!opt.engine_root.empty())
      params.setCpmdRoot(opt.engine_root);
    if (!opt.title.empty())
      params.setTitle(opt.title);
    if (opt.memory_mb > 0)
      params.setMemoryMb(static_cast<uint32_t>(opt.memory_mb));
    if (!opt.scratch_dir.empty())
      params.setScratchDir(opt.scratch_dir);
    impl_->cpmd = std::make_unique<rgpot::CPMDPot>(params.asReader());
    if (!impl_->cpmd->available())
      throw std::runtime_error(
          "RGPOT(cpmdc): engine not available (set CPMDC_LIBRARY / "
          "RGPOT_CPMDC_ENGINE or [RgpotPot] engine_path)");
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
    throw std::runtime_error(
        "RGPOT: unknown backend '" + opt.backend +
        "' (expected nwchemc, cpmdc, metatomic, or xtb)");
  }
}

RGPotEngine::~RGPotEngine() = default;

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

  std::pair<double, AtomMatrix> result;
  if (impl_->backend == Impl::Backend::Nwchemc)
    result = (*impl_->nwchem)(positions, atmtypes, cell);
  else
    result = (*impl_->cpmd)(positions, atmtypes, cell);

  *U = result.first;
  const auto &forces = result.second;
  for (long i = 0; i < N; ++i) {
    const int ii = static_cast<int>(i);
    F[3 * i + 0] = forces(ii, 0);
    F[3 * i + 1] = forces(ii, 1);
    F[3 * i + 2] = forces(ii, 2);
  }
}
