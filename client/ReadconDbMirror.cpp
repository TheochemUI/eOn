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
#include "ReadconDbMirror.h"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#if !defined(_WIN32)
#include <dlfcn.h>
#endif

namespace eonc::io {
namespace {

constexpr std::uint64_t kFnvOffset = 14695981039346656037ull;
constexpr std::uint64_t kFnvPrime = 1099511628211ull;
constexpr std::size_t kCorpusGroup = 64;

std::uint64_t fnv1a(std::string_view bytes) {
  std::uint64_t h = kFnvOffset;
  for (unsigned char c : bytes) {
    h ^= static_cast<std::uint64_t>(c);
    h *= kFnvPrime;
  }
  return h == 0 ? 1 : h;
}

std::filesystem::path corpus_dir(const std::filesystem::path &con_path) {
  namespace fs = std::filesystem;
  std::error_code ec;
  fs::path cur = fs::weakly_canonical(con_path, ec);
  if (ec) {
    cur = fs::absolute(con_path);
  }
  cur = cur.parent_path();
  for (int i = 0; i < 8; ++i) {
    if (fs::is_regular_file(cur / "config.ini")) {
      return cur / "readcon.db";
    }
    const fs::path parent = cur.parent_path();
    if (parent == cur) {
      break;
    }
    cur = parent;
  }
  fs::path beside = fs::weakly_canonical(con_path, ec);
  if (ec) {
    beside = fs::absolute(con_path);
  }
  return beside.parent_path() / "readcon.db";
}

std::string canonical_path(const std::filesystem::path &con_path) {
  namespace fs = std::filesystem;
  std::error_code ec;
  fs::path canonical = fs::weakly_canonical(con_path, ec);
  if (ec) {
    canonical = fs::absolute(con_path, ec);
    if (ec) {
      canonical = con_path;
    }
  }
  return canonical.string();
}

#if !defined(_WIN32)
using OpenFn = int (*)(const char *, std::size_t *);
using AppendFn = int (*)(std::size_t, std::uint64_t, const char *, const char *,
                         std::uint32_t *);
using ExtendPathFn = int (*)(std::size_t, std::uint64_t, const char *,
                             std::uint32_t *);
using ExtendStrFn = int (*)(std::size_t, std::uint64_t, const char *,
                            const char *, std::uint32_t *);

struct Api {
  void *lib = nullptr;
  OpenFn open = nullptr;
  AppendFn append_str = nullptr;
  ExtendPathFn extend_path = nullptr;
  ExtendStrFn extend_str = nullptr;
  bool missing = false;
};

Api &api() {
  static Api out;
  if (out.lib != nullptr || out.missing) {
    return out;
  }
  const char *names[] = {"libreadcon_db.so", "libreadcon_db.so.0",
                         "libreadcon_db.dylib"};
  for (const char *name : names) {
    out.lib = dlopen(name, RTLD_LAZY | RTLD_LOCAL);
    if (out.lib != nullptr) {
      break;
    }
  }
  if (out.lib == nullptr) {
    out.missing = true;
    return out;
  }
  out.open = reinterpret_cast<OpenFn>(dlsym(out.lib, "rkrdb_open"));
  out.append_str =
      reinterpret_cast<AppendFn>(dlsym(out.lib, "rkrdb_append_trajectory_str"));
  out.extend_path =
      reinterpret_cast<ExtendPathFn>(dlsym(out.lib, "rkrdb_extend_trajectory"));
  out.extend_str = reinterpret_cast<ExtendStrFn>(
      dlsym(out.lib, "rkrdb_extend_trajectory_str"));
  if (out.open == nullptr) {
    dlclose(out.lib);
    out.lib = nullptr;
    out.append_str = nullptr;
    out.extend_path = nullptr;
    out.extend_str = nullptr;
    out.missing = true;
  }
  return out;
}

std::mutex &handles_mu() {
  static std::mutex mu;
  return mu;
}

std::unordered_map<std::string, std::size_t> &handles() {
  static std::unordered_map<std::string, std::size_t> map;
  return map;
}

std::size_t corpus_handle(const std::filesystem::path &dir) {
  Api &fns = api();
  if (fns.missing || fns.open == nullptr) {
    return static_cast<std::size_t>(-1);
  }
  const std::string key = dir.string();
  std::lock_guard<std::mutex> guard(handles_mu());
  const auto found = handles().find(key);
  if (found != handles().end()) {
    return found->second;
  }
  std::size_t id = 0;
  if (fns.open(key.c_str(), &id) != 0) {
    return static_cast<std::size_t>(-1);
  }
  handles().emplace(key, id);
  return id;
}

struct Batch {
  std::string canonical;
  std::filesystem::path db;
  std::vector<std::string> frames;
};

std::atomic<std::uint64_t> scratch_counter{0};

std::mutex &batch_mu() {
  static std::mutex mu;
  return mu;
}

std::unordered_map<std::string, Batch> &batches() {
  static std::unordered_map<std::string, Batch> map;
  return map;
}

void extend_batch(const std::string &canonical, const std::filesystem::path &db,
                  const std::vector<std::string> &frames) {
  if (frames.empty()) {
    return;
  }
  Api &fns = api();
  if (fns.missing || (fns.extend_str == nullptr && fns.extend_path == nullptr &&
                      fns.append_str == nullptr)) {
    return;
  }
  const std::size_t id = corpus_handle(db);
  if (id == static_cast<std::size_t>(-1)) {
    return;
  }
  std::string text;
  std::size_t bytes = 0;
  for (const auto &frame : frames) {
    bytes += frame.size();
  }
  text.reserve(bytes);
  for (const auto &frame : frames) {
    text += frame;
  }
  const std::uint64_t tid = fnv1a(canonical);
  std::uint32_t nframes = 0;
  if (fns.extend_str != nullptr) {
    fns.extend_str(id, tid, text.c_str(), canonical.c_str(), &nframes);
    return;
  }
  if (fns.extend_path != nullptr) {
    const auto scratch =
        std::filesystem::temp_directory_path() /
        std::format("eon-corpus-{}-{}.con", tid,
                    scratch_counter.fetch_add(1, std::memory_order_relaxed));
    {
      std::ofstream out(scratch, std::ios::binary | std::ios::trunc);
      if (!out) {
        return;
      }
      out.write(text.data(), static_cast<std::streamsize>(text.size()));
      if (!out) {
        std::error_code ec;
        std::filesystem::remove(scratch, ec);
        return;
      }
    }
    fns.extend_path(id, tid, scratch.string().c_str(), &nframes);
    std::error_code ec;
    std::filesystem::remove(scratch, ec);
  }
}

void flush_ready(std::vector<Batch> ready) {
  for (const auto &batch : ready) {
    extend_batch(batch.canonical, batch.db, batch.frames);
  }
}
#endif

} // namespace

std::uint64_t corpus_path_id(std::string_view canonical_path) {
  return fnv1a(canonical_path);
}

void arm_corpus_flush() {
#if !defined(_WIN32)
  // The exit handler must run while these caches are still alive.
  // A handler registered before they are constructed runs after they
  // are destroyed.
  (void)batch_mu();
  (void)batches();
  (void)handles_mu();
  (void)handles();
  (void)api();
  static std::once_flag once;
  std::call_once(once, [] { std::atexit(flush_corpus_appends); });
#endif
}

void mirror_con_corpus(const std::string &path) {
#if !defined(_WIN32)
  arm_corpus_flush();
  namespace fs = std::filesystem;
  std::error_code ec;
  if (!fs::is_regular_file(path, ec) || ec) {
    return;
  }
  std::ifstream in(path);
  if (!in) {
    return;
  }
  std::ostringstream buf;
  buf << in.rdbuf();
  const std::string text = buf.str();
  if (text.find_first_not_of(" \t\r\n") == std::string::npos) {
    return;
  }
  const std::string resolved = canonical_path(path);
  Api &fns = api();
  if (fns.missing || fns.append_str == nullptr) {
    return;
  }
  const std::size_t id = corpus_handle(corpus_dir(resolved));
  if (id == static_cast<std::size_t>(-1)) {
    return;
  }
  const std::uint64_t tid = fnv1a(resolved);
  std::uint32_t nframes = 0;
  fns.append_str(id, tid, text.c_str(), resolved.c_str(), &nframes);
#else
  (void)path;
#endif
}

void note_corpus_append(const std::string &path,
                        const std::string &frame_text) {
#if !defined(_WIN32)
  arm_corpus_flush();
  if (frame_text.find_first_not_of(" \t\r\n") == std::string::npos) {
    return;
  }
  const std::string resolved = canonical_path(path);
  std::vector<Batch> ready;
  {
    std::lock_guard<std::mutex> guard(batch_mu());
    Batch &slot = batches()[resolved];
    slot.canonical = resolved;
    slot.db = corpus_dir(resolved);
    slot.frames.push_back(frame_text);
    if (slot.frames.size() >= kCorpusGroup) {
      ready.push_back(std::move(slot));
      batches().erase(resolved);
    }
  }
  flush_ready(std::move(ready));
#else
  (void)path;
  (void)frame_text;
#endif
}

void flush_corpus_appends() {
#if !defined(_WIN32)
  std::vector<Batch> ready;
  {
    std::lock_guard<std::mutex> guard(batch_mu());
    ready.reserve(batches().size());
    for (auto &entry : batches()) {
      if (!entry.second.frames.empty()) {
        ready.push_back(std::move(entry.second));
      }
    }
    batches().clear();
  }
  flush_ready(std::move(ready));
#endif
}

std::size_t corpus_buffered_frames(const std::string &path) {
#if !defined(_WIN32)
  const std::string resolved = canonical_path(path);
  std::lock_guard<std::mutex> guard(batch_mu());
  const auto found = batches().find(resolved);
  if (found == batches().end()) {
    return 0;
  }
  return found->second.frames.size();
#else
  (void)path;
  return 0;
#endif
}

} // namespace eonc::io
