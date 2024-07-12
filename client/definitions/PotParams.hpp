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
#include "../EnvHelpers.hpp"
#include <string>
using namespace std::string_literals;

namespace eonc::def {
using namespace helper_functions;
struct XTBParams {
  double acc{1.0};
  double elec_temperature{0.0};
  size_t maxiter{250};
  std::string paramset{"GNFF"s};
};

struct ASEOrcaParams {
  std::string orca_path;
  std::string simpleinput;
  std::string orca_nproc;
  ASEOrcaParams()
      : orca_path{get_value_from_env_or_param("ORCA_COMMAND", PDef(""s, ""s),
                                              "", false)},
        simpleinput{get_value_from_env_or_param(
            "ORCA_SIMPLEINPUT", PDef("ENGRAD HF-3c", "ENGRAD HF-3c"),
            "Using ENGRAD HF-3c as a default input, set simpleinput or the "
            "environment variable ORCA_SIMPLEINPUT.\n")},
        orca_nproc{get_value_from_env_or_param("ORCA_NPROC", PDef("1"s, "auto"),
                                               "", false)} {}
};

struct AMS_ENV {
  std::string amshome;
  std::string scm_tmpdir;
  std::string scmlicense;
  std::string scm_pythondir;
  std::string amsbin;
  std::string amsresources;

  AMS_ENV() {
    amshome = get_value_from_env_or_param("AMSHOME", PDef("", ""), "", false);
    scm_tmpdir =
        get_value_from_env_or_param("SCM_TMPDIR", PDef("", "/tmp"), "", false);
    scmlicense =
        get_value_from_env_or_param("SCMLICENSE", PDef("", ""), "", false);
    scm_pythondir =
        get_value_from_env_or_param("SCM_PYTHONDIR", PDef("", ""), "", false);
    amsbin = get_value_from_env_or_param("AMSBIN", PDef("", ""), "", false);
    amsresources =
        get_value_from_env_or_param("AMSRESOURCES", PDef("", ""), "", false);
  }
};

struct AMSParams {
  std::string engine;     // MOPAC, ADF, BAND, REAXFF, FORCEFIELD
  std::string forcefield; // OPt.ff etc. (REAXFF)
  std::string model;      // Model hamiltonian (MOPAC)
  std::string resources;  // DFTB
  std::string xc;         // Exchange (BAND, ADF)
  std::string basis;      // Basis (BAND, ADF)
  AMS_ENV amsenv;
  AMSParams()
      : engine{""s},
        forcefield{""s},
        model{""s},
        resources{""s},
        xc{""s},
        basis{""s},
        amsenv{AMS_ENV()} {}
};

// [CatLearn] //
// No reasonable default for catl_path
struct CatLearnParams {
  std::string path;
  std::string model;
  std::string prior;
  bool use_deriv;
  bool use_fingerprint;
  bool parallel;
  CatLearnParams()
      : path{""s},
        model{"gp"s},
        prior{"median"s},
        use_deriv{true},
        use_fingerprint{false},
        parallel{true} {}
};

} // namespace eonc::def
