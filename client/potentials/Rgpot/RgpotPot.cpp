#include "RgpotPot.h"
#include "RgpotClientRpc.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

RgpotPot::RgpotPot(const Parameters &p)
    : Potential(PotType::RGPOT, p) {
  const auto &o = p.rgpot_options;
  host_ = o.host;
  port_ = o.port;
  backend_ = o.backend;
  nwchem_basis_ = o.nwchem_basis;
  nwchem_theory_ = o.nwchem_theory;
  nwchem_scf_type_ = o.nwchem_scf_type;
  nwchem_charge_ = o.nwchem_charge;
  nwchem_multiplicity_ = o.nwchem_multiplicity;
  cpmd_functional_ = o.cpmd_functional;
  cpmd_task_ = o.cpmd_task;
  cpmd_cut_off_ry_ = o.cpmd_cut_off_ry;
  cpmd_charge_ = o.cpmd_charge;
  cpmd_multiplicity_ = o.cpmd_multiplicity;

  // Env overrides for CI / scripts
  if (const char *e = std::getenv("RGPOT_POTSERV_HOST")) {
    host_ = e;
  }
  if (const char *e = std::getenv("RGPOT_POTSERV_PORT")) {
    port_ = std::atoi(e);
  }
  if (const char *e = std::getenv("RGPOT_BACKEND")) {
    backend_ = e;
  }

  std::cout << "RgpotPot: client to " << host_ << ":" << port_
            << " backend=" << backend_ << std::endl;
  client_holder_ = rgpotClientCreate(host_, port_);
}

RgpotPot::~RgpotPot() {
  if (client_holder_) {
    rgpotClientDestroy(client_holder_);
    client_holder_ = nullptr;
  }
}

void RgpotPot::ensureConfigured() {
  if (configured_) {
    return;
  }
  RgpotConfigureSpec spec;
  spec.backend = backend_;
  spec.nwchem_basis = nwchem_basis_;
  spec.nwchem_theory = nwchem_theory_;
  spec.nwchem_scf_type = nwchem_scf_type_;
  spec.nwchem_charge = nwchem_charge_;
  spec.nwchem_multiplicity = nwchem_multiplicity_;
  spec.cpmd_functional = cpmd_functional_;
  spec.cpmd_task = cpmd_task_;
  spec.cpmd_cut_off_ry = cpmd_cut_off_ry_;
  spec.cpmd_charge = cpmd_charge_;
  spec.cpmd_multiplicity = cpmd_multiplicity_;

  std::string msg;
  const bool ok = rgpotClientConfigure(client_holder_, spec, &msg);
  if (!ok) {
    throw std::runtime_error("RgpotPot configure failed: " + msg);
  }
  std::cout << "RgpotPot: configured ok (" << msg << ")" << std::endl;
  configured_ = true;
}

void RgpotPot::force(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, double *variance, const double *box) {
  (void)variance;
  ensureConfigured();
  double box9[9];
  // eOn box is 3x3 row-major in R
  for (int i = 0; i < 9; ++i) {
    box9[i] = box[i];
  }
  rgpotClientCalculate(client_holder_, N, R, atomicNrs, box9, U, F);
}
