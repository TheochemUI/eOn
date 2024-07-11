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
#include "Parameters.h"
#include <iostream>
#include <memory>

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

  // Does not take into account the fixed / free atoms
  // Variance here is null when not needed and that's OK
  void virtual force(long nAtoms, const double *positions, const int *atomicNrs,
                     double *forces, double *energy, double *variance,
                     const double *box) = 0;

  std::tuple<double, AtomMatrix>
  get_ef(const AtomMatrix pos, const Vector<int> atmnrs, const Matrix3S box);
};

template <typename T> class Potential : public PotBase {
public:
  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    return static_cast<T *>(this)->force(nAtoms, positions, atomicNrs, forces,
                                         energy, variance, box);
  }
};

std::shared_ptr<PotBase> makePotential(const toml::table &config);
} // namespace eonc
