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

// TODO(rg) :: These defines need to be tested for and set in meson.build
#define HAVE_ALIGNED_ALLOC 1
#define HAVE_POSIX_MEMALIGN 1
#define CACHELOT_PLATFORM_BITS 64
#include <cachelot/cache.h>
#include <cachelot/common.h>

namespace eonc::cache {
// TODO(rg) :: Make these parameters
constexpr size_t cache_memory = 64 * cachelot::Megabyte;
constexpr size_t page_size = 4 * cachelot::Megabyte;
constexpr size_t hash_initial = 131072;

struct KeyHash {
  size_t hash;
  cachelot::slice key;
  // TODO(rg): Slightly inefficient to do 2 string conversions for no good
  // reason
  KeyHash(size_t _hash)
      : hash{_hash},
        key(std::to_string(_hash).c_str(), std::to_string(_hash).size()) {}
  // .. but slightly ugly to use KeyHash(std::to_string(hash_in), hash_in)
  KeyHash(std::string _key, size_t _hash)
      : hash{_hash},
        key(_key.c_str(), _key.size()) {}
};

// Helper which doesn't compute the hashes, since those are obtained based on
// the Potential<T>
// TODO(rg) :: Consider reworking as a singleton like the logger
class PotentialCache {
private:
  cachelot::cache::Cache potCache;
  static std::atomic<PotentialCache *> instance;
  static std::mutex pcMutex;

  PotentialCache()
      : potCache(cachelot::cache::Cache::Create(
            eonc::cache::cache_memory, eonc::cache::page_size,
            eonc::cache::hash_initial, true)) {}

  ~PotentialCache() = default;

  PotentialCache(const PotentialCache &) = delete;
  PotentialCache &operator=(const PotentialCache &) = delete;

  void setup_cache();

public:
  static PotentialCache *getInstance();
  void deserialize_hit(cachelot::cache::ConstItemPtr &, ForceOut &,
                       AtomMatrix &);
  void add_serialized(const KeyHash &, const ForceOut &, const AtomMatrix &);
  cachelot::cache::ConstItemPtr find(const KeyHash &);
};

} // namespace eonc::cache
