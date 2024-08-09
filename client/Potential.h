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
#include "client/thirdparty/toml.hpp"
#include <memory>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

#include <xxhash.h>

#include "client/potentials/PotentialCache.hpp"

namespace eonc {
class PotBase {
protected:
  eonc::cache::PotentialCache *pcache = nullptr;
  // Does not take into account the fixed / free atoms
  // Callers responsibility that ForceOut has enough space!!!
  // Protected since this is really only to be implemented... callers use get_ef
  virtual void force(const ForceInput &, ForceOut *) = 0;

public:
  virtual std::tuple<double, AtomMatrix>
  get_ef(const AtomMatrix &, const Vector<size_t> &, const Matrix3S &) = 0;
  virtual size_t getInstances() const = 0;
  virtual size_t getTotalForceCalls() const = 0;
  void set_cache(eonc::cache::PotentialCache *cc) { pcache = cc; }
};

template <typename T>
class Potential : public PotBase, public pot::registry<T> {
protected:
  void force(const ForceInput &params, ForceOut *efvd) override final {
    pot::registry<T>::incrementForceCalls();
    static_cast<T *>(this)->forceImpl(params, efvd);
    return;
  }

private:
  size_t computeHash(const AtomMatrix &pos,
                     const Vector<size_t> &atmnrs) const {
    XXH64_state_t *state = XXH64_createState();
    if (state == nullptr) {
      throw std::runtime_error("Failed to create XXH64 state");
    }

    // TODO(rg):: seed should be a user parameter
    XXH64_reset(state, 1995);

    XXH64_update(state, reinterpret_cast<const void *>(pos.data()),
                 sizeof(double) * pos.size());
    XXH64_update(state, reinterpret_cast<const void *>(atmnrs.data()),
                 sizeof(size_t) * atmnrs.size());
    const std::string potential_name = typeid(T).name();
    XXH64_update(state, potential_name.data(), potential_name.size());

    // Finalize the hash computation and free state
    size_t hash_result = XXH64_digest(state);
    XXH64_freeState(state);

    return hash_result;
  }

public:
  // Mostly for testing
  size_t getHash(const AtomMatrix &pos, const Vector<size_t> &atmnrs) const {
    return computeHash(pos, atmnrs);
  }
  // To be implemented by the child classes
  virtual void forceImpl(const ForceInput &params, ForceOut *efvd) = 0;
  // Safer, saner returns, and also allocates memory for force()
  // TODO(rg):: A variant return can unify SurrogatePotential and Potential
  std::tuple<double, AtomMatrix> get_ef(const AtomMatrix &pos,
                                        const Vector<size_t> &atmnrs,
                                        const Matrix3S &box) override final {
    const size_t nAtoms{static_cast<size_t>(pos.rows())};
    // When not in debug mode the initial values are unchecked
    // So the initial data in efd matters!
    AtomMatrix forces{MatrixType::Zero(nAtoms, 3)};
    ForceOut efd{forces.data(), 0, 0};
    if (this->pcache != nullptr) {
      SPDLOG_TRACE("Cache present");
      handle_cache(pos, atmnrs, box, efd, forces);
    } else {
      SPDLOG_TRACE("Cache not present");
      this->force({nAtoms, pos.data(), atmnrs.data(), box.data()}, &efd);
    }
    return std::make_tuple(efd.energy, forces);
  };

  void handle_cache(const AtomMatrix &pos, const Vector<size_t> &atmnrs,
                    const Matrix3S &box, ForceOut &efd, AtomMatrix &forces) {
    // TODO(rg):: This can probably be parsed faster
    const size_t nAtoms{static_cast<size_t>(pos.rows())};
    const size_t currentHash = getHash(pos, atmnrs);
    SPDLOG_TRACE("Current hash and key is {}", currentHash);
    const eonc::cache::KeyHash kv{currentHash};
    auto found_item = pcache->find(kv);
    if (found_item) {
      SPDLOG_TRACE("Cache hit");
      pcache->deserialize_hit(found_item, efd, forces);
      SPDLOG_TRACE("Found {}", efd.energy);
    } else {
      SPDLOG_TRACE("Cache miss");
      this->force({nAtoms, pos.data(), atmnrs.data(), box.data()}, &efd);
      pcache->add_serialized(kv, efd, forces);
      SPDLOG_TRACE("Not found in cache, so added {}", efd.energy);
    }
  }

  size_t getInstances() const override final { return pot::registry<T>::count; }
  size_t getTotalForceCalls() const override final {
    return pot::registry<T>::forceCalls;
  }
};

std::shared_ptr<PotBase> makePotential(const toml::table &config);
} // namespace eonc
