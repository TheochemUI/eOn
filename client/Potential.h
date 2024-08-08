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
#include "client/thirdparty/xxhash.hpp"
#include <memory>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

// TODO(rg) :: These defines need to be tested for and set in meson.build
#define HAVE_ALIGNED_ALLOC 1
#define HAVE_POSIX_MEMALIGN 1
#define CACHELOT_PLATFORM_BITS 64
#include <cachelot/cache.h>
#include <cachelot/common.h>

namespace eonc {
constexpr size_t cache_memory = 64 * cachelot::Megabyte;
constexpr size_t page_size = 4 * cachelot::Megabyte;
constexpr size_t hash_initial = 131072;

class PotBase {
protected:
  cachelot::cache::Cache *potCache;
  // Does not take into account the fixed / free atoms
  // Callers responsibility that ForceOut has enough space!!!
  // Protected since this is really only to be implemented... callers use get_ef
  virtual void force(const ForceInput &, ForceOut *) = 0;

public:
  virtual std::tuple<double, AtomMatrix>
  get_ef(const AtomMatrix &, const Vector<size_t> &, const Matrix3S &) = 0;
  virtual size_t getInstances() const = 0;
  virtual size_t getTotalForceCalls() const = 0;
  void set_cache(cachelot::cache::Cache *cc) { potCache = cc; }
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
    xxh::hash_state_t<32> hash_stream;
    for (auto idx = 0; idx < pos.size(); ++idx) {
      hash_stream.update(reinterpret_cast<const char *>(&pos.data()[idx]),
                         sizeof(pos.data()[idx]));
    }
    for (auto idx = 0; idx < atmnrs.size(); ++idx) {
      hash_stream.update(reinterpret_cast<const char *>(&atmnrs[idx]),
                         sizeof(atmnrs[idx]));
    }
    const std::string type_name = typeid(T).name();
    hash_stream.update(type_name);
    return hash_stream.digest();
  }

public:
  // To be implemented by the child classes
  virtual void forceImpl(const ForceInput &params, ForceOut *efvd) = 0;
  // Safer, saner returns, and also allocates memory for force()
  // TODO(rg):: A variant return can unify SurrogatePotential and Potential
  std::tuple<double, AtomMatrix> get_ef(const AtomMatrix &pos,
                                        const Vector<size_t> &atmnrs,
                                        const Matrix3S &box) override final {
    const size_t nAtoms{static_cast<size_t>(pos.rows())};
    const size_t currentHash = computeHash(pos, atmnrs);
    // When not in debug mode the initial values are unchecked
    // So the initial data in efd matters!
    AtomMatrix forces{MatrixType::Zero(nAtoms, 3)};
    ForceOut efd{forces.data(), 0, 0};
    if (this->potCache != nullptr) {
      SPDLOG_INFO("Cache present and hit");
      const auto chash = std::to_string(currentHash);
      cachelot::slice cache_key(chash.c_str(), chash.size());
      SPDLOG_INFO("Current Hash is {}, with key {}", currentHash,
                  cache_key.str());

      auto found_item = potCache->do_get(cache_key, currentHash);
      if (found_item) {
        // Deserialize the cached value
        auto buffer = found_item->value().str();
        std::memcpy(&efd.energy, buffer.c_str(), sizeof(double));
        std::memcpy(forces.data(), buffer.c_str() + sizeof(double),
                    forces.size() * sizeof(double));
        std::cout << "Found " << efd.energy << std::endl;
      } else {
        SPDLOG_INFO("Cache present and miss");
        this->force({nAtoms, pos.data(), atmnrs.data(), box.data()}, &efd);

        // Serialize the value
        size_t value_size = sizeof(double) + forces.size() * sizeof(double);
        std::vector<char> buffer(value_size);
        std::memcpy(buffer.data(), &efd.energy, sizeof(double));
        std::memcpy(buffer.data() + sizeof(double), forces.data(),
                    forces.size() * sizeof(double));

        cachelot::slice value_slice(buffer.data(), buffer.size());
        auto new_item =
            potCache->create_item(cache_key, currentHash, value_slice.length(),
                                  0, cachelot::cache::Item::infinite_TTL);
        new_item->assign_value(value_slice);
        bool isDone = potCache->do_add(new_item);
        SPDLOG_INFO("{} :: Not found, so added {}", isDone, efd.energy);
        if (not isDone) {
          throw std::runtime_error("Key collision for Potential cache");
        }
      }
    } else {
      SPDLOG_INFO("Cache not present");
      this->force({nAtoms, pos.data(), atmnrs.data(), box.data()}, &efd);
    }
    return std::make_tuple(efd.energy, forces);
  };

  size_t getInstances() const override final { return pot::registry<T>::count; }
  size_t getTotalForceCalls() const override final {
    return pot::registry<T>::forceCalls;
  }
};

std::shared_ptr<PotBase> makePotential(const toml::table &config);
} // namespace eonc
