/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#pragma once

#ifdef RGPOT_HAS_CACHE
#include "rgpot/PotentialCache.hpp"
#include <stdexcept>
#include <string>

namespace eonc {

class RgpotCacheProvider {
public:
  virtual ~RgpotCacheProvider() = default;
  virtual void setResultCache(rgpot::cache::PotentialCache *cache) = 0;
};

/// The potential borrows the database for the lifetime of this job.
class ScopedRgpotCache {
public:
  ScopedRgpotCache(RgpotCacheProvider &provider, const std::string &path)
      : cache_(path), provider_(provider) {
    checkOpen(cache_, path);
    provider_.setResultCache(&cache_);
  }
  ~ScopedRgpotCache() { provider_.setResultCache(nullptr); }
  ScopedRgpotCache(const ScopedRgpotCache &) = delete;
  ScopedRgpotCache &operator=(const ScopedRgpotCache &) = delete;

  auto counts() const { return cache_.counts(); }

private:
  template <class Cache>
  static void checkOpen(const Cache &cache, const std::string &path) {
    if constexpr (requires { cache.is_open(); }) {
      if (!cache.is_open()) {
        throw std::runtime_error("lineb: cannot open rgpot result cache " +
                                 path);
      }
    } else {
      throw std::runtime_error(
          "lineb: installed rgpot must expose PotentialCache::is_open to "
          "verify result_cache");
    }
  }
  rgpot::cache::PotentialCache cache_;
  RgpotCacheProvider &provider_;
};

} // namespace eonc
#endif
