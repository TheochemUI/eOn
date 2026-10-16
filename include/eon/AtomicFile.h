// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cstdio>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace eonc::helpers {

// The destination names only complete, flushed files. The temporary file lives
// beside it so replacement cannot cross a filesystem boundary.
template <typename Writer>
bool atomicWriteFile(const std::string &path, Writer writer) {
  std::string pattern = path + ".tmp.XXXXXX";
  std::vector<char> name(pattern.begin(), pattern.end());
  name.push_back('\0');
#ifdef _WIN32
  if (_mktemp_s(name.data(), name.size()) != 0) {
    return false;
  }
  int fd = _open(name.data(), _O_CREAT | _O_EXCL | _O_WRONLY | _O_BINARY,
                 _S_IREAD | _S_IWRITE);
#else
  int fd = mkstemp(name.data());
#endif
  if (fd < 0) {
    return false;
  }
#ifdef _WIN32
  FILE *file = _fdopen(fd, "wb");
#else
  FILE *file = fdopen(fd, "wb");
#endif
  if (!file) {
#ifdef _WIN32
    _close(fd);
#else
    close(fd);
#endif
    std::remove(name.data());
    return false;
  }
  bool ok = false;
  try {
    ok = writer(file) && !std::ferror(file) && std::fflush(file) == 0;
  } catch (...) {
    std::fclose(file);
    std::remove(name.data());
    throw;
  }
#ifdef _WIN32
  if (ok) {
    ok = _commit(fd) == 0;
  }
#else
  if (ok) {
    ok = fsync(fd) == 0;
  }
#endif
  if (std::fclose(file) != 0) {
    ok = false;
  }
  if (ok) {
#ifdef _WIN32
    ok = MoveFileExA(name.data(), path.c_str(),
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) != 0;
#else
    ok = std::rename(name.data(), path.c_str()) == 0;
#endif
  }
  if (!ok) {
    std::remove(name.data());
  }
  return ok;
}

} // namespace eonc::helpers
