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

#include "Eigen.h"
#include "thirdparty/toml.hpp"
#include <memory>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  // pointer to number of atoms, pointer to array of positions
  // address to supercell size
  const size_t nAtoms;
  const double *pos;
  const size_t *atmnrs;
  const double *box;
} ForceInput;

typedef struct {
  // pointer to array of forces
  double *F;
  // Internal energy
  double energy;
  // Variance here is 0 when not needed and that's OK
  double variance;
} ForceOut;

#ifdef __cplusplus
}
#endif

namespace eonc {

class PotBase {
protected:
  std::shared_ptr<spdlog::logger> m_log;
  void initializeLogger() {
    if (!spdlog::get("_potcalls")) {
      // Create logger if it doesn't exist
      m_log = spdlog::basic_logger_mt("_potcalls", "_potcalls.log", true);
      m_log->set_pattern("[%l] [%Y-%m-%d %H:%M:%S] %v");
    } else {
      // Use existing logger
      m_log = spdlog::get("_potcalls");
    }
    // TODO(rg): Fix this
    // if (m_log) {
    //   m_log->trace("[{}] created", typeid(this).name());
    // }
  }

  // Does not take into account the fixed / free atoms
  // Callers responsibility that ForceOut has enough space!!!
  // Protected since this is really only to be implemented... callers use get_ef
  void virtual force(const ForceInput &params, ForceOut *efvd) = 0;

public:
  PotBase() {
    SPDLOG_TRACE("CREATED WITH {}", forceCallCounter);
    initializeLogger();
  }
  size_t forceCallCounter{0};
  virtual ~PotBase() {
    SPDLOG_TRACE("DESTROYED AFTER {}\n", forceCallCounter);
    // TODO(rg): Fix this
    // if (m_log) {
    //   m_log->trace("[{}] destroyed after {} calls", typeid(this).name(),
    //                forceCallCounter);
    // } else {
    //   std::cerr << "Logger is not initialized\n";
    // }
  }

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
