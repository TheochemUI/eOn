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

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>

#ifndef EON_READCON_DB_VERSION
#define EON_READCON_DB_VERSION ""
#endif

#if !defined(_WIN32)
#include <dlfcn.h>

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

std::uint64_t fnv1a(const std::string &bytes) {
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

#if !defined(_WIN32)
using OpenFn = int (*)(const char *, std::size_t *);
using AppendFn = int (*)(std::size_t, std::uint64_t, const char *, const char *,
                         std::uint32_t *);

struct Api {
  void *lib = nullptr;
  OpenFn open = nullptr;
  AppendFn append_str = nullptr;
  bool missing = false;
  bool from_process = false;
  int epoch = -1;
};

int &mirror_epoch() {
  static int epoch = 0;
  return epoch;
}

void clear_api(Api &out) {
  if (out.lib != nullptr && !out.from_process) {
    dlclose(out.lib);
  }
  out = Api{};
}

bool resolve_symbols(Api &out, void *handle) {
  out.open = reinterpret_cast<OpenFn>(dlsym(handle, "rkrdb_open"));
  out.append_str = reinterpret_cast<AppendFn>(
      dlsym(handle, "rkrdb_append_trajectory_str"));
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
    out.lib = dlopen(forced, RTLD_LAZY | RTLD_LOCAL);
  } else if (resolve_symbols(out, RTLD_DEFAULT)) {
    out.from_process = true;
    return out;
  } else {
    out.open = nullptr;
    out.append_str = nullptr;
    const char *names[] = {"libreadcon_db.so", "libreadcon_db.so.0",
                           "libreadcon_db.dylib"};
    for (const char *name : names) {
      out.lib = dlopen(name, RTLD_LAZY | RTLD_LOCAL);
      if (out.lib != nullptr) {
        break;
      }
    }
  }
  if (out.lib == nullptr) {
    const char *err = dlerror();
    EONC_LOG_WARNING("[readcon-db] libreadcon_db.so failed to load: {}",
                     err != nullptr ? err : "dlopen returned null");
    out.missing = true;
    return out;
  }
  if (!resolve_symbols(out, out.lib)) {
    dlclose(out.lib);
    out.lib = nullptr;
    out.open = nullptr;
    out.append_str = nullptr;
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
#endif

} // namespace

void mirror_con_corpus(const std::string &path) {
#if !defined(_WIN32)
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
  fs::path canonical = fs::weakly_canonical(path, ec);
  if (ec) {
    canonical = fs::absolute(path);
  }
  const std::string resolved = canonical.string();
  const auto corpus = corpus_dir(canonical);
  if (!corpus) {
    return;
  }
  const std::size_t id = corpus_handle(*corpus);
  if (id == static_cast<std::size_t>(-1)) {
    return;
  }
  const std::uint64_t tid = fnv1a(resolved + '\0' + text);
  std::uint32_t nframes = 0;
  api().append_str(id, tid, text.c_str(), resolved.c_str(), &nframes);
#else
  static_cast<void>(path);
#endif
}

bool readcon_db_mirror_ok() {
#if !defined(_WIN32)
  const Api &fns = api();
  return !fns.missing && fns.open != nullptr && fns.append_str != nullptr;
#else
  return false;
#endif
}

void readcon_db_mirror_reset() {
#if !defined(_WIN32)
  mirror_epoch()++;
#endif
}

const char *readcon_db_loaded_version() {
  if (!readcon_db_mirror_ok()) {
    return "";
  }
  return EON_READCON_DB_VERSION;
}

} // namespace eonc::io
