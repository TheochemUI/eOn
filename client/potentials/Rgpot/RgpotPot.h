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
#pragma once

#include "../../Potential.h"
#include <memory>
#include <string>

/**
 * @brief Cap'n Proto client potential talking to rgpot potserv (NWChem / CPMD).
 *
 * Does not include Potentials.capnp.h (name collision with Potential). RPC lives
 * in RgpotClientRpc.cpp. Units on the wire: angstrom / eV (eOn force() convention).
 */
class RgpotPot : public Potential {
public:
  explicit RgpotPot(const Parameters &p);
  ~RgpotPot() override;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  void ensureConfigured();

  std::string host_;
  int port_{12345};
  std::string backend_{"NWChem"};

  std::string nwchem_basis_{"sto-3g"};
  std::string nwchem_theory_{"scf"};
  std::string nwchem_scf_type_{"rhf"};
  int nwchem_charge_{0};
  int nwchem_multiplicity_{1};

  std::string cpmd_functional_{"BLYP"};
  std::string cpmd_task_{"gradient"};
  double cpmd_cut_off_ry_{70.0};
  int cpmd_charge_{0};
  int cpmd_multiplicity_{1};

  bool configured_{false};
  // Opaque client held as void* so header stays free of Cap'n Proto types.
  void *client_holder_{nullptr};
};
