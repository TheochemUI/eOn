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
#include "PotRegistry.h"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace eonc {

namespace {
std::string format_timepoint(PotRegistry::TimePoint tp) {
  auto time_t_val = PotRegistry::Clock::to_time_t(tp);
  auto us = std::chrono::duration_cast<std::chrono::microseconds>(
                tp.time_since_epoch()) %
            std::chrono::seconds{1};
  std::tm tm_buf{};
#ifdef _WIN32
  localtime_s(&tm_buf, &time_t_val);
#else
  localtime_r(&time_t_val, &tm_buf);
#endif
  std::ostringstream oss;
  oss << std::put_time(&tm_buf, "%Y-%m-%dT%H:%M:%S");
  oss << '.' << std::setfill('0') << std::setw(6) << us.count();
  return oss.str();
}

// Escape a string for JSON (handles quotes and backslashes)
std::string json_escape(std::string_view sv) {
  std::string out;
  out.reserve(sv.size());
  for (char c : sv) {
    if (c == '"' || c == '\\')
      out += '\\';
    out += c;
  }
  return out;
}
} // namespace

PotRegistry &PotRegistry::get() noexcept {
  static PotRegistry instance;
  return instance;
}

void PotRegistry::reset() {
  for (auto &ts : m_type_stats) {
    ts.force_calls.store(0, std::memory_order_relaxed);
    ts.created.store(0, std::memory_order_relaxed);
    ts.alive.store(0, std::memory_order_relaxed);
  }
  std::lock_guard<std::mutex> lock(m_records_mutex);
  m_records.clear();
  m_next_id.store(1, std::memory_order_relaxed);
}

uint64_t PotRegistry::on_created(PotType t) noexcept {
  auto idx = static_cast<size_t>(magic_enum::enum_index(t).value_or(0));
  m_type_stats[idx].created.fetch_add(1, std::memory_order_relaxed);
  m_type_stats[idx].alive.fetch_add(1, std::memory_order_relaxed);
  return m_next_id.fetch_add(1, std::memory_order_relaxed);
}

void PotRegistry::on_destroyed(uint64_t id, PotType t, size_t force_calls,
                               TimePoint created_at) {
  auto idx = static_cast<size_t>(magic_enum::enum_index(t).value_or(0));
  m_type_stats[idx].alive.fetch_sub(1, std::memory_order_relaxed);

  InstanceRecord rec{id, t, created_at, Clock::now(), force_calls};
  std::lock_guard<std::mutex> lock(m_records_mutex);
  m_records.push_back(rec);
}

void PotRegistry::on_force_call(PotType t) noexcept {
  auto idx = static_cast<size_t>(magic_enum::enum_index(t).value_or(0));
  m_type_stats[idx].force_calls.fetch_add(1, std::memory_order_relaxed);
}

size_t PotRegistry::type_force_calls(PotType t) const noexcept {
  auto idx = static_cast<size_t>(magic_enum::enum_index(t).value_or(0));
  return m_type_stats[idx].force_calls.load(std::memory_order_relaxed);
}

size_t PotRegistry::total_force_calls() const noexcept {
  size_t total = 0;
  for (const auto &ts : m_type_stats) {
    total += ts.force_calls.load(std::memory_order_relaxed);
  }
  return total;
}

size_t PotRegistry::type_alive(PotType t) const noexcept {
  auto idx = static_cast<size_t>(magic_enum::enum_index(t).value_or(0));
  return m_type_stats[idx].alive.load(std::memory_order_relaxed);
}

void PotRegistry::write_summary(const std::string &path) const {
  std::vector<InstanceRecord> snapshot;
  {
    std::lock_guard<std::mutex> lock(const_cast<std::mutex &>(m_records_mutex));
    snapshot = m_records;
  }

  std::ofstream ofs(path);
  if (!ofs.is_open())
    return;

  ofs << "[\n";
  for (size_t i = 0; i < snapshot.size(); ++i) {
    const auto &r = snapshot[i];
    auto type_name = magic_enum::enum_name(r.type);
    ofs << "  {\n";
    ofs << "    \"id\": " << r.id << ",\n";
    ofs << "    \"type\": \"" << json_escape(type_name) << "\",\n";
    ofs << "    \"created_at\": \"" << format_timepoint(r.created_at)
        << "\",\n";
    ofs << "    \"destroyed_at\": \"" << format_timepoint(r.destroyed_at)
        << "\",\n";
    ofs << "    \"force_calls\": " << r.force_calls << "\n";
    ofs << "  }";
    if (i + 1 < snapshot.size())
      ofs << ',';
    ofs << '\n';
  }
  ofs << "]\n";
}

} // namespace eonc
