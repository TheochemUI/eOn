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

#include "eon/EonLogger.h"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <mutex>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef EON_READCON_DB_VERSION
#define EON_READCON_DB_VERSION ""
#endif

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
using LibHandle = HMODULE;
#else
#include <dlfcn.h>
using LibHandle = void *;

extern "C" int rkrdb_open(const char *, std::size_t *);
// The mirror looks the symbol up with dlsym. Naming it here keeps the
// linked libreadcon_db.so on the link line under --as-needed.
[[maybe_unused]] static int (*const kLinkedOpen)(const char *,
                                                 std::size_t *) = rkrdb_open;
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

std::optional<std::filesystem::path>
corpus_dir(const std::filesystem::path &con_path) {
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
  // A stray .con is not a job. Do not mint readcon.db beside it.
  return std::nullopt;
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

using OpenFn = int (*)(const char *, std::size_t *);
using AppendFn = int (*)(std::size_t, std::uint64_t, const char *, const char *,
                         std::uint32_t *);
using ExtendPathFn = int (*)(std::size_t, std::uint64_t, const char *,
                             std::uint32_t *);
using ExtendStrFn = int (*)(std::size_t, std::uint64_t, const char *,
                            const char *, std::uint32_t *);

struct Api {
  LibHandle lib = nullptr;
  OpenFn open = nullptr;
  AppendFn append_str = nullptr;
  ExtendPathFn extend_path = nullptr;
  ExtendStrFn extend_str = nullptr;
  bool missing = false;
  bool from_process = false;
  int epoch = -1;
};

int &mirror_epoch() {
  static int epoch = 0;
  return epoch;
}

void release_lib(LibHandle lib, bool from_process) {
  if (lib == nullptr || from_process) {
    return;
  }
#if defined(_WIN32)
  FreeLibrary(lib);
#else
  dlclose(lib);
#endif
}

void *symbol(LibHandle handle, const char *name) {
#if defined(_WIN32)
  return reinterpret_cast<void *>(GetProcAddress(handle, name));
#else
  return dlsym(handle, name);
#endif
}

LibHandle open_named(const char *name) {
#if defined(_WIN32)
  return LoadLibraryA(name);
#else
  return dlopen(name, RTLD_LAZY | RTLD_LOCAL);
#endif
}

const char *open_error() {
#if defined(_WIN32)
  return "LoadLibrary failed";
#else
  const char *err = dlerror();
  return err != nullptr ? err : "dlopen returned null";
#endif
}

void clear_api(Api &out) {
  release_lib(out.lib, out.from_process);
  out = Api{};
}

void clear_symbols(Api &out) {
  out.open = nullptr;
  out.append_str = nullptr;
  out.extend_path = nullptr;
  out.extend_str = nullptr;
}

// The extend entry points are optional. Older readcon-db builds lack them
// and appended frames then skip the corpus.
bool resolve_symbols(Api &out, LibHandle handle) {
  out.open = reinterpret_cast<OpenFn>(symbol(handle, "rkrdb_open"));
  out.append_str =
      reinterpret_cast<AppendFn>(symbol(handle, "rkrdb_append_trajectory_str"));
  out.extend_path =
      reinterpret_cast<ExtendPathFn>(symbol(handle, "rkrdb_extend_trajectory"));
  out.extend_str = reinterpret_cast<ExtendStrFn>(
      symbol(handle, "rkrdb_extend_trajectory_str"));
  return out.open != nullptr && out.append_str != nullptr;
}

Api &api() {
  static Api out;
  if (out.epoch == mirror_epoch() &&
      (out.lib != nullptr || out.from_process || out.missing)) {
    return out;
  }
  clear_api(out);
  out.epoch = mirror_epoch();
  const char *forced = std::getenv("EON_READCON_DB_LIBRARY");
  if (forced != nullptr) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (forced[0] == '\0' || !fs::is_regular_file(forced, ec) || ec) {
      EONC_LOG_WARNING("[readcon-db] libreadcon_db.so failed to load: {}",
                       forced[0] == '\0' ? "empty path" : forced);
      out.missing = true;
      return out;
    }
    out.lib = open_named(forced);
  } else {
#if !defined(_WIN32)
    if (resolve_symbols(out, RTLD_DEFAULT)) {
      out.from_process = true;
      return out;
    }
    clear_symbols(out);
    const char *names[] = {"libreadcon_db.so", "libreadcon_db.so.0",
                           "libreadcon_db.dylib"};
#else
    // The import library does not pull the DLL in unless some object
    // references rkrdb_*. Load it from PATH (the test runner puts the
    // build directory there) or from an earlier LoadLibrary.
    HMODULE loaded = GetModuleHandleA("readcon_db.dll");
    if (loaded != nullptr && resolve_symbols(out, loaded)) {
      out.lib = loaded;
      out.from_process = true;
      return out;
    }
    clear_symbols(out);
    const char *names[] = {"readcon_db.dll"};
#endif
    for (const char *name : names) {
      out.lib = open_named(name);
      if (out.lib != nullptr) {
        break;
      }
    }
  }
  if (out.lib == nullptr) {
    EONC_LOG_WARNING("[readcon-db] libreadcon_db.so failed to load: {}",
                     open_error());
    out.missing = true;
    return out;
  }
  if (!resolve_symbols(out, out.lib)) {
    release_lib(out.lib, false);
    out.lib = nullptr;
    clear_symbols(out);
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
  if (fns.missing) {
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
  if (fns.missing ||
      (fns.extend_str == nullptr && fns.extend_path == nullptr)) {
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

void flush_ready(std::vector<Batch> ready) {
  for (const auto &batch : ready) {
    extend_batch(batch.canonical, batch.db, batch.frames);
  }
}

void arm_corpus_flush() {
  // The exit handler must run while these caches are still alive.
  // A handler registered before they are constructed runs after they
  // are destroyed.
  static_cast<void>(batch_mu());
  static_cast<void>(batches());
  static_cast<void>(handles_mu());
  static_cast<void>(handles());
  static_cast<void>(api());
  static std::once_flag once;
  std::call_once(once, [] { std::atexit(flush_corpus_appends); });
}

} // namespace

std::uint64_t corpus_path_id(std::string_view canonical_path) {
  return fnv1a(canonical_path);
}

std::uint64_t corpus_blob_id(std::string_view canonical_path,
                             std::string_view text) {
  std::string bytes;
  bytes.reserve(canonical_path.size() + 1 + text.size());
  bytes.append(canonical_path);
  bytes.push_back('\0');
  bytes.append(text);
  return fnv1a(bytes);
}

void mirror_con_corpus(const std::string &path) {
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
  const auto corpus = corpus_dir(resolved);
  if (!corpus) {
    return;
  }
  const std::size_t id = corpus_handle(*corpus);
  if (id == static_cast<std::size_t>(-1)) {
    return;
  }
  const std::uint64_t tid = corpus_blob_id(resolved, text);
  std::uint32_t nframes = 0;
  api().append_str(id, tid, text.c_str(), resolved.c_str(), &nframes);
}

void note_corpus_append(const std::string &path,
                        const std::string &frame_text) {
  arm_corpus_flush();
  if (frame_text.find_first_not_of(" \t\r\n") == std::string::npos) {
    return;
  }
  const std::string resolved = canonical_path(path);
  const auto corpus = corpus_dir(resolved);
  if (!corpus) {
    return;
  }
  std::vector<Batch> ready;
  {
    std::lock_guard<std::mutex> guard(batch_mu());
    Batch &slot = batches()[resolved];
    slot.canonical = resolved;
    slot.db = *corpus;
    slot.frames.push_back(frame_text);
    if (slot.frames.size() >= kCorpusGroup) {
      ready.push_back(std::move(slot));
      batches().erase(resolved);
    }
  }
  flush_ready(std::move(ready));
}

void flush_corpus_appends() {
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
}

std::size_t corpus_buffered_frames(const std::string &path) {
  const std::string resolved = canonical_path(path);
  std::lock_guard<std::mutex> guard(batch_mu());
  const auto found = batches().find(resolved);
  if (found == batches().end()) {
    return 0;
  }
  return found->second.frames.size();
}

bool corpus_extend_ok() {
  const Api &fns = api();
  return !fns.missing && fns.open != nullptr &&
         (fns.extend_str != nullptr || fns.extend_path != nullptr);
}

bool readcon_db_mirror_ok() {
  const Api &fns = api();
  return !fns.missing && fns.open != nullptr && fns.append_str != nullptr;
}

void readcon_db_mirror_reset() { mirror_epoch()++; }

const char *readcon_db_loaded_version() {
  if (!readcon_db_mirror_ok()) {
    return "";
  }
  return EON_READCON_DB_VERSION;
}

} // namespace eonc::io
