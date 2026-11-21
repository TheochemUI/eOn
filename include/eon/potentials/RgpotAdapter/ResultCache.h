/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#pragma once

#ifdef RGPOT_HAS_CACHE
#include "rgpot/PotentialCache.hpp"
#include <memory>
#include <stdexcept>
#include <string>

namespace eonc {

class RgpotCacheProvider {
public:
  virtual ~RgpotCacheProvider() = default;
  virtual void setResultCache(rgpot::cache::PotentialCache *cache) = 0;
};

/// A dynamic rgpot engine owns its own database across the C ABI.
class RgpotEngineCacheProvider {
public:
  virtual ~RgpotEngineCacheProvider() = default;
  virtual void attachEngineCache(const std::string &path) = 0;
  virtual void detachEngineCache() noexcept = 0;
  virtual rgpot::cache::EvalCounts engineCacheCounts() const = 0;
};

/// The potential borrows the database for the lifetime of this job.
class ScopedRgpotCache {
public:
  ScopedRgpotCache(RgpotCacheProvider &provider, const std::string &path)
      : cache_(std::make_unique<rgpot::cache::PotentialCache>(path)),
        provider_(&provider) {
    checkOpen(*cache_, path);
    provider_->setResultCache(cache_.get());
  }
  ScopedRgpotCache(RgpotEngineCacheProvider &engine, const std::string &path)
      : engine_(&engine) {
    engine_->attachEngineCache(path);
  }
  ~ScopedRgpotCache() {
    if (engine_)
      engine_->detachEngineCache();
    else
      provider_->setResultCache(nullptr);
  }
  ScopedRgpotCache(const ScopedRgpotCache &) = delete;
  ScopedRgpotCache &operator=(const ScopedRgpotCache &) = delete;

  auto counts() const {
    return engine_ ? engine_->engineCacheCounts() : cache_->counts();
  }

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
  std::unique_ptr<rgpot::cache::PotentialCache> cache_;
  RgpotCacheProvider *provider_{};
  RgpotEngineCacheProvider *engine_{};
};

} // namespace eonc
#endif
