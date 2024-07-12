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

#include "C_Structs.h"
#include "Eigen.h"
#include "thirdparty/toml.hpp"
#include <memory>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

namespace eonc {

class PotBase {
protected:
  std::shared_ptr<spdlog::logger> m_log;
  // Does not take into account the fixed / free atoms
  // Callers responsibility that ForceOut has enough space!!!
  // Protected since this is really only to be implemented... callers use get_ef
  void virtual force(const ForceInput &params, ForceOut *efvd) = 0;

public:
  // Safer, saner returns, and also allocates memory for force()
  // TODO(rg):: A variant return can unify SurrogatePotential and Potential
  std::tuple<double, AtomMatrix>
  get_ef(const AtomMatrix pos, const Vector<size_t> atmnrs, const Matrix3S box);
};

template <typename T> class Potential : public PotBase {
public:
  void force(const ForceInput &params, ForceOut *efvd) override {
    return static_cast<T *>(this)->force(params, efvd);
  }
};

std::shared_ptr<PotBase> makePotential(const toml::table &config);
} // namespace eonc
