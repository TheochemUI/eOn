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
#include "XtbLoader.h"

#include <iostream>

namespace eonc {

XtbLoader &XtbLoader::instance() {
  static XtbLoader loader;
  return loader;
}

XtbLoader::XtbLoader() {
#ifdef _WIN32
  const char *names[] = {"xtb.dll", "libxtb.dll", "libxtb-6.dll", nullptr};
#elif defined(__APPLE__)
  const char *names[] = {"libxtb.dylib", "libxtb.6.dylib", nullptr};
#else
  const char *names[] = {"libxtb.so", "libxtb.so.6", nullptr};
#endif

  m_handle = dynlib::openFirst(names);
  if (!m_handle) {
    return; // Not found; require_loaded() raises if a caller asks for XTB.
  }

  new_environment =
      dynlib::loadSym<new_environment_fn>(m_handle, "xtb_newEnvironment");
  del_environment =
      dynlib::loadSym<del_environment_fn>(m_handle, "xtb_delEnvironment");
  check_environment =
      dynlib::loadSym<check_environment_fn>(m_handle, "xtb_checkEnvironment");
  get_error = dynlib::loadSym<get_error_fn>(m_handle, "xtb_getError");
  release_output =
      dynlib::loadSym<release_output_fn>(m_handle, "xtb_releaseOutput");
  set_verbosity =
      dynlib::loadSym<set_verbosity_fn>(m_handle, "xtb_setVerbosity");

  new_molecule = dynlib::loadSym<new_molecule_fn>(m_handle, "xtb_newMolecule");
  del_molecule = dynlib::loadSym<del_molecule_fn>(m_handle, "xtb_delMolecule");
  update_molecule =
      dynlib::loadSym<update_molecule_fn>(m_handle, "xtb_updateMolecule");

  new_calculator =
      dynlib::loadSym<new_calculator_fn>(m_handle, "xtb_newCalculator");
  del_calculator =
      dynlib::loadSym<del_calculator_fn>(m_handle, "xtb_delCalculator");
  load_gfnff = dynlib::loadSym<load_gfnff_fn>(m_handle, "xtb_loadGFNFF");
  load_gfn0 = dynlib::loadSym<load_gfn0_fn>(m_handle, "xtb_loadGFN0xTB");
  load_gfn1 = dynlib::loadSym<load_gfn1_fn>(m_handle, "xtb_loadGFN1xTB");
  load_gfn2 = dynlib::loadSym<load_gfn2_fn>(m_handle, "xtb_loadGFN2xTB");
  set_accuracy =
      dynlib::loadSym<set_accuracy_fn>(m_handle, "xtb_setAccuracy");
  set_max_iter = dynlib::loadSym<set_max_iter_fn>(m_handle, "xtb_setMaxIter");
  set_electronic_temp = dynlib::loadSym<set_electronic_temp_fn>(
      m_handle, "xtb_setElectronicTemp");

  singlepoint =
      dynlib::loadSym<singlepoint_fn>(m_handle, "xtb_singlepoint");

  new_results = dynlib::loadSym<new_results_fn>(m_handle, "xtb_newResults");
  del_results = dynlib::loadSym<del_results_fn>(m_handle, "xtb_delResults");
  get_energy = dynlib::loadSym<get_energy_fn>(m_handle, "xtb_getEnergy");
  get_gradient =
      dynlib::loadSym<get_gradient_fn>(m_handle, "xtb_getGradient");

  // Required minimum surface; the GFN parametrisation loaders are
  // checked at use site so users can run on a libxtb missing one of
  // the four (e.g. distributors stripping GFNFF).
  if (!new_environment || !del_environment || !check_environment ||
      !get_error || !release_output || !set_verbosity || !new_molecule ||
      !del_molecule || !update_molecule || !new_calculator ||
      !del_calculator || !singlepoint || !new_results || !del_results ||
      !get_energy || !get_gradient || !set_accuracy || !set_max_iter ||
      !set_electronic_temp) {
    std::cerr << "[XTB] libxtb loaded but lacks required symbols\n";
    dynlib::close(m_handle);
    m_handle = {};
    return;
  }

  m_loaded = true;
}

XtbLoader::~XtbLoader() { dynlib::close(m_handle); }

void XtbLoader::require_loaded() const {
  if (!m_loaded) {
    throw std::runtime_error(
        "XTB potential requested but libxtb not found.\n"
        "Install via: conda install -c conda-forge xtb\n"
        "Or ensure libxtb is in your library search path "
        "(LD_LIBRARY_PATH / DYLD_LIBRARY_PATH / PATH).");
  }
}

} // namespace eonc
