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

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>

#if !defined(_WIN32)
#include <dlfcn.h>
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

#if !defined(_WIN32)
using OpenFn = int (*)(const char *, std::size_t *);
using AppendFn = int (*)(std::size_t, std::uint64_t, const char *, const char *,
                         std::uint32_t *);

struct Api {
  void *lib = nullptr;
  OpenFn open = nullptr;
  AppendFn append_str = nullptr;
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
  if (out.open == nullptr || out.append_str == nullptr) {
    dlclose(out.lib);
    out.lib = nullptr;
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
  const std::size_t id = corpus_handle(corpus_dir(canonical));
  if (id == static_cast<std::size_t>(-1)) {
    return;
  }
  const std::uint64_t tid = fnv1a(resolved + '\0' + text);
  std::uint32_t nframes = 0;
  api().append_str(id, tid, text.c_str(), resolved.c_str(), &nframes);
#else
  (void)path;
#endif
}

} // namespace eonc::io
