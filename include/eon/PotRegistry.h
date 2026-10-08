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
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace eonc {

/// Outlives the registry object. A Potential copies this pointer and
/// skips on_destroyed once the registry has started to tear down.
struct RegistryLifetime {
  std::atomic<bool> alive{true};
};

// ---------------------------------------------------------------------------
// IPotRegistry: injectable ABI. Production uses PotRegistry::get().
// ---------------------------------------------------------------------------
class IPotRegistry {
public:
  using Clock = std::chrono::system_clock;
  using TimePoint = Clock::time_point;

  virtual ~IPotRegistry() {
    if (lifetime_)
      lifetime_->alive.store(false, std::memory_order_release);
  }

  [[nodiscard]] std::shared_ptr<RegistryLifetime> lifetime() const {
    return lifetime_;
  }
  [[nodiscard]] virtual uint64_t on_created(PotType t) noexcept = 0;
  virtual void on_destroyed(uint64_t id, PotType t, size_t force_calls,
                            TimePoint created_at) = 0;
  virtual void on_force_call(PotType t) noexcept = 0;

  IPotRegistry(const IPotRegistry &) = delete;
  IPotRegistry &operator=(const IPotRegistry &) = delete;

protected:
  IPotRegistry() = default;

  std::shared_ptr<RegistryLifetime> lifetime_{
      std::make_shared<RegistryLifetime>()};
};

class PotRegistry : public IPotRegistry {
public:
  using Clock = IPotRegistry::Clock;
  using TimePoint = IPotRegistry::TimePoint;

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
  PotRegistry() = default;
  ~PotRegistry() override {
    // Members are still intact. Mark the registry dead before they go,
    // so a Potential destroyed from here does not lock a dead mutex.
    if (lifetime_)
      lifetime_->alive.store(false, std::memory_order_release);
  }

  /// Process-lifetime singleton. The instance is allocated on the heap
  /// and never destroyed, so Potential destructors can still record
  /// teardown after C++ static destruction.
  static PotRegistry &get() noexcept;
  void reset();

  // Lifecycle events
  [[nodiscard]] uint64_t on_created(PotType t) noexcept override;
  void on_destroyed(uint64_t id, PotType t, size_t force_calls,
                    TimePoint created_at) override;
  void on_force_call(PotType t) noexcept override;

  // Per-type queries
  [[nodiscard]] size_t type_force_calls(PotType t) const noexcept;
  [[nodiscard]] size_t total_force_calls() const noexcept;
  [[nodiscard]] size_t type_alive(PotType t) const noexcept;

  // JSON output
  void write_summary(const std::string &path = "_potcalls.json") const;
};

} // namespace eonc
