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
#include "client/potentials/PotentialCache.hpp"

namespace eonc::cache {
void PotentialCache::set_cache(cachelot::cache::Cache *cc) { potCache = cc; }
void PotentialCache::deserialize_hit(cachelot::cache::ConstItemPtr &hit,
                                     ForceOut &efd, AtomMatrix &forces) {
  // NOTE(rg) :: Assumes efd, and AtomMatrix are setup
  auto buffer = hit->value().str();
  std::memcpy(&efd.energy, buffer.c_str(), sizeof(double));
  std::memcpy(forces.data(), buffer.c_str() + sizeof(double),
              forces.size() * sizeof(double));
}
void PotentialCache::add_serialized(const KeyHash &kv, const ForceOut &efd,
                                    const AtomMatrix &forces) {
  // Serialize the value
  size_t value_size = sizeof(double) + forces.size() * sizeof(double);
  std::vector<char> buffer(value_size);
  std::memcpy(buffer.data(), &efd.energy, sizeof(double));
  std::memcpy(buffer.data() + sizeof(double), forces.data(),
              forces.size() * sizeof(double));

  cachelot::slice value_slice(buffer.data(), buffer.size());
  auto new_item = potCache->create_item(kv.key, kv.hash, value_slice.length(),
                                        0, cachelot::cache::Item::infinite_TTL);
  new_item->assign_value(value_slice);
  bool isDone = potCache->do_add(new_item);
  if (not isDone) {
    throw std::runtime_error("Key collision for Potential cache");
  }
}
cachelot::cache::ConstItemPtr PotentialCache::find(const KeyHash &kv) {
  return potCache->do_get(kv.key, kv.hash);
}
} // namespace eonc::cache
