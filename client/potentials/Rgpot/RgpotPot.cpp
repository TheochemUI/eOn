/*
** This file is part of eOn.
*/
#include "RgpotPot.h"
#include "RGPotEngine.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <string>

RgpotPot::RgpotPot(const Parameters &p)
    : Potential(PotType::RGPOT, p) {
  RGPotEngineOptions opt;
  const auto &o = p.rgpot_options;
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
  // Engine-path env overrides are backend-scoped: NWCHEMC_LIBRARY must not
  // leak into a cpmdc configure (CPMDPot resolves CPMDC_LIBRARY itself).
  std::string backend_lc = opt.backend;
  std::transform(backend_lc.begin(), backend_lc.end(), backend_lc.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (backend_lc.rfind("nwchem", 0) == 0) {
    if (const char *e = std::getenv("NWCHEMC_LIBRARY"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("RGPOT_NWCHEMC_ENGINE"))
      opt.engine_path = e;
  } else if (backend_lc.rfind("cpmd", 0) == 0) {
    if (const char *e = std::getenv("CPMDC_LIBRARY"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("RGPOT_CPMDC_ENGINE"))
      opt.engine_path = e;
  } else if (backend_lc == "xtb" || backend_lc == "xtbpot" ||
             backend_lc == "gfn" || backend_lc == "gfnxtb") {
    if (const char *e = std::getenv("RGPOT_XTB_ENGINE"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("XTB_ENGINE"))
      opt.engine_path = e;
    // Prefer dedicated [XTBPot]/xtb_options when still default-ish from RGPOT
    if (opt.xtb_paramset.empty() || opt.xtb_paramset == "GFN2xTB") {
      if (!p.xtb_options.paramset.empty())
        opt.xtb_paramset = p.xtb_options.paramset;
    }
  }

  impl_ = std::make_unique<RGPotEngine>(opt);
  backend_ = impl_->backend();
  std::cout << "RgpotPot: in-process rgpot backend=" << backend_
            << " (dlopen engines: libnwchemc/libcpmdc/libxtb_engine)"
            << std::endl;
}

RgpotPot::~RgpotPot() = default;

bool RgpotPot::engineAvailable() const { return impl_ && impl_->available(); }

void RgpotPot::force(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, double *variance, const double *box) {
  (void)variance;
  impl_->force(N, R, atomicNrs, F, U, box);
}
