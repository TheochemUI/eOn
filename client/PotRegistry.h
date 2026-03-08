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

#include "BaseStructures.h"
#include <array>
#include <atomic>
#include <chrono>
#include <mutex>
#include <string>
#include <vector>

namespace eonc {

class PotRegistry {
public:
  using Clock = std::chrono::system_clock;
  using TimePoint = Clock::time_point;

  struct InstanceRecord {
    uint64_t id;
    PotType type;
    TimePoint created_at;
    TimePoint destroyed_at;
    size_t force_calls;
  };

private:
  struct TypeStats {
    std::atomic<size_t> force_calls{0};
    std::atomic<size_t> created{0};
    std::atomic<size_t> alive{0};
  };
  std::array<TypeStats, magic_enum::enum_count<PotType>()> m_type_stats{};

  std::mutex m_records_mutex;
  std::vector<InstanceRecord> m_records;

  std::atomic<uint64_t> m_next_id{1};

public:
  static PotRegistry &get() noexcept;
  void reset();

  // Lifecycle events
  [[nodiscard]] uint64_t on_created(PotType t) noexcept;
  void on_destroyed(uint64_t id, PotType t, size_t force_calls,
                    TimePoint created_at);
  void on_force_call(PotType t) noexcept;

  // Per-type queries
  [[nodiscard]] size_t type_force_calls(PotType t) const noexcept;
  [[nodiscard]] size_t total_force_calls() const noexcept;
  [[nodiscard]] size_t type_alive(PotType t) const noexcept;

  // JSON output
  void write_summary(const std::string &path = "_potcalls.json") const;
};

} // namespace eonc

using eonc::PotRegistry;
