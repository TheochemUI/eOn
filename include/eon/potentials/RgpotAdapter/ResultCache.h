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

/// The potential borrows the database for the lifetime of this job.
class ScopedRgpotCache {
public:
  ScopedRgpotCache(RgpotCacheProvider &provider, const std::string &path)
      : provider_(provider) {
    rocksdb::Options options;
    options.create_if_missing = true;
#if RGPOT_ROCKSDB_UNIQUE_PTR_OPEN
    const auto status = rocksdb::DB::Open(options, path, &db_);
#else
    rocksdb::DB *raw = nullptr;
    const auto status = rocksdb::DB::Open(options, path, &raw);
    db_.reset(raw);
#endif
    if (!status.ok()) {
      throw std::runtime_error("lineb: cannot open rgpot result cache " + path +
                               ": " + status.ToString());
    }
    cache_.set_db(db_.get());
    provider_.setResultCache(&cache_);
  }
  ~ScopedRgpotCache() { provider_.setResultCache(nullptr); }
  ScopedRgpotCache(const ScopedRgpotCache &) = delete;
  ScopedRgpotCache &operator=(const ScopedRgpotCache &) = delete;

  auto counts() const { return cache_.counts(); }

private:
  std::unique_ptr<rocksdb::DB> db_;
  rgpot::cache::PotentialCache cache_;
  RgpotCacheProvider &provider_;
};

} // namespace eonc
#endif
