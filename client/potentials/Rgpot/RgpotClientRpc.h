/*
** Cap'n Proto RPC helpers for RgpotPot (no eOn Potential.h include).
*/
#pragma once

#include <string>
#include <vector>

struct RgpotConfigureSpec {
  std::string backend; // "NWChem" or "CPMD"
  std::string nwchem_basis;
  std::string nwchem_theory;
  std::string nwchem_scf_type;
  int nwchem_charge{0};
  int nwchem_multiplicity{1};
  /// Optional NWChem input blocks (e.g. "dft\\n  xc b3lyp\\nend") for reliable XC
  /// when host libnwchemc does not promote scfType alone.
  std::vector<std::string> nwchem_input_blocks;
  std::string cpmd_functional;
  std::string cpmd_task;
  double cpmd_cut_off_ry{70.0};
  int cpmd_charge{0};
  int cpmd_multiplicity{1};
};

// Opaque EzRpcClient + wait scope storage; destroy with rgpotClientDestroy.
void *rgpotClientCreate(const std::string &host, int port);
void rgpotClientDestroy(void *holder);

// Returns true if configure() reported ok.
bool rgpotClientConfigure(void *holder, const RgpotConfigureSpec &spec,
                          std::string *message_out);

// Positions/box in Angstrom; energy eV; forces eV/Angstrom; box is 9 doubles row-major.
void rgpotClientCalculate(void *holder, long n_atoms, const double *positions,
                          const int *atomic_nrs, const double *box_9,
                          double *energy_out, double *forces_out);
