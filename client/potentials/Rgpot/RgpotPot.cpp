/*
** This file is part of eOn.
*/
#include "RgpotPot.h"
#include "RGPotEngine.h"

#include <cstdlib>
#include <iostream>

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

  // Env overrides (CI / benchmarks)
  if (const char *e = std::getenv("RGPOT_BACKEND"))
    opt.backend = e;
  if (const char *e = std::getenv("RGPOT_NWCHEM_BASIS"))
    opt.basis = e;
  if (const char *e = std::getenv("RGPOT_NWCHEM_THEORY"))
    opt.theory = e;
  if (const char *e = std::getenv("RGPOT_NWCHEM_SCF_TYPE"))
    opt.scf_type = e;
  if (const char *e = std::getenv("NWCHEMC_LIBRARY"))
    opt.engine_path = e;
  else if (const char *e = std::getenv("RGPOT_NWCHEMC_ENGINE"))
    opt.engine_path = e;

  impl_ = std::make_unique<RGPotEngine>(opt);
  backend_ = impl_->backend();
  std::cout << "RgpotPot: in-process rgpot backend=" << backend_
            << " (dlopen libnwchemc/libcpmdc, no potserv RPC)" << std::endl;
}

RgpotPot::~RgpotPot() = default;

bool RgpotPot::engineAvailable() const { return impl_ && impl_->available(); }

void RgpotPot::force(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, double *variance, const double *box) {
  (void)variance;
  impl_->force(N, R, atomicNrs, F, U, box);
}
