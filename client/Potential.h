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

#include "client/C_Structs.h"
#include "client/Eigen.h"
#include "client/potentials/PotHelpers.hpp"
#include "thirdparty/toml.hpp"
#include <memory>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

namespace eonc {

class PotBase {
protected:
  // Does not take into account the fixed / free atoms
  // Callers responsibility that ForceOut has enough space!!!
  // Protected since this is really only to be implemented... callers use get_ef
  virtual void force(const ForceInput &params, ForceOut *efvd) = 0;

public:
  virtual std::tuple<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                                const Vector<size_t> atmnrs,
                                                const Matrix3S box) = 0;
  virtual size_t getInstances() const = 0;
  virtual size_t getTotalForceCalls() const = 0;
};

template <typename T>
class Potential : public PotBase, public pot::registry<T> {
protected:
  void force(const ForceInput &params, ForceOut *efvd) override final {
    pot::registry<T>::incrementForceCalls();
    return static_cast<T *>(this)->forceImpl(params, efvd);
  }

public:
  // To be implemented by the child classes
  virtual void forceImpl(const ForceInput &params, ForceOut *efvd) = 0;
  // Safer, saner returns, and also allocates memory for force()
  // TODO(rg):: A variant return can unify SurrogatePotential and Potential
  std::tuple<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                        const Vector<size_t> atmnrs,
                                        const Matrix3S box) override final {
    size_t nAtoms{static_cast<size_t>(pos.rows())};
    // When not in debug mode the initial values are unchecked
    // So the initial data in efd matters!
    AtomMatrix forces{MatrixType::Zero(nAtoms, 3)};
    ForceOut efd{forces.data(), 0, 0};
    this->force({nAtoms, pos.data(), atmnrs.data(), box.data()}, &efd);
    return std::make_tuple(efd.energy, forces);
  };

  size_t getInstances() const override final { return pot::registry<T>::count; }
  size_t getTotalForceCalls() const override final {
    return pot::registry<T>::forceCalls;
  }
};

std::shared_ptr<PotBase> makePotential(const toml::table &config);
} // namespace eonc
