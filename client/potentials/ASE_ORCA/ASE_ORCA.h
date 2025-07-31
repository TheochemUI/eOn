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
#define PYBIND11_DETAILED_ERROR_MESSAGES

#include "../../EnvHelpers.hpp"
#include "../../Potential.h"

#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>

using namespace std::string_literals;
namespace eonc {

namespace py = pybind11;
using namespace pybind11::literals; // to bring in the `_a` literal

class ASEOrcaPot final : public Potential<ASEOrcaPot> {
public:
  struct Params final {
    std::string orca_path{def::get_value_from_env_or_param(
        "ORCA_COMMAND", def::PDef(""s, ""s), "", false)};
    std::string simpleinput{get_value_from_env_or_param(
        "ORCA_SIMPLEINPUT", def::PDef("ENGRAD HF-3c", "ENGRAD HF-3c"),
        "Using ENGRAD HF-3c as a default input, set simpleinput or the "
        "environment variable ORCA_SIMPLEINPUT.\n")};
    std::string orca_nproc{def::get_value_from_env_or_param(
        "ORCA_NPROC", def::PDef("1"s, "auto"), "", false)};
  };

private:
  py::object _calc;
  py::object _ase;

public:
  ASEOrcaPot(const ASEOrcaPot::Params &);
  void setParameters(const ASEOrcaPot::Params &);

  // Functions
  void forceImpl(const ForceInput &, ForceOut *) override final;
};

} // namespace eonc
