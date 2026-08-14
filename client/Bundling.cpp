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

#include "eon/Bundling.h"
#include "eon/EonLogger.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <system_error>
#include <utility>

namespace fs = std::filesystem;

namespace eonc {

int getBundleSize() {
  int num_bundle = -1;

  for (const auto &entry : fs::directory_iterator(".")) {
    std::string name = entry.path().filename().string();

    if (name[0] == '.') {
      continue;
    }

    // If "config" is not in the filename
    // then skip.
    if (name.find("config") == std::string::npos &&
        name.find("ini") == std::string::npos) {
      continue;
    }

    // Find the last underscore
    auto upos = name.rfind('_');
    if (upos == std::string::npos) {
      continue;
    }

    // Find the last period
    auto dpos = name.rfind('.');
    if (dpos == std::string::npos || dpos <= upos) {
      continue;
    }

    std::string numstr = name.substr(upos + 1, dpos - upos - 1);
    if (!numstr.empty() &&
        std::isdigit(static_cast<unsigned char>(numstr[0]))) {
      int i = std::atoi(numstr.c_str()) + 1;
      if (i > num_bundle) {
        num_bundle = i;
      }
    }
  }

  return num_bundle;
}

int strchrcount(const char *haystack, char needle) {
  int count = 0;
  for (const char *ch = haystack; *ch != '\0'; ch++) {
    if (*ch == needle) {
      count++;
    }
  }
  return count;
}

std::vector<std::string> unbundle(int number) {
  std::vector<std::string> filenames;

  for (const auto &entry : fs::directory_iterator(".")) {
    std::string originalFilename = entry.path().filename().string();

    if (originalFilename[0] == '.') {
      continue;
    }

    int numUnderscores = strchrcount(originalFilename.c_str(), '_');
    if (numUnderscores < 1) {
      continue;
    }

    // Find the last underscore
    auto upos = originalFilename.rfind('_');
    if (upos == std::string::npos) {
      continue;
    }

    // Find the last period
    auto dpos = originalFilename.rfind('.');
    if (dpos == std::string::npos || dpos <= upos) {
      continue;
    }

    std::string numstr = originalFilename.substr(upos + 1, dpos - upos - 1);
    if (!numstr.empty() &&
        std::isdigit(static_cast<unsigned char>(numstr[0]))) {
      int bundleNumber = std::atoi(numstr.c_str());
      if (bundleNumber != number) {
        continue;
      }
    }

    std::string baseName = originalFilename.substr(0, upos);
    // Drop the "_<n>" token and keep everything after it, so a compound
    // extension (".con.gz", ".con.zst") survives whichever side of the
    // underscore it sits on.
    size_t extStart =
        originalFilename.find_first_not_of("0123456789", upos + 1);
    if (extStart == std::string::npos) {
      extStart = originalFilename.size();
    }
    std::string newFilename = baseName + originalFilename.substr(extStart);

    try {
      fs::copy_file(originalFilename, newFilename,
                    fs::copy_options::overwrite_existing);
      filenames.push_back(newFilename);
    } catch (const fs::filesystem_error &e) {
      EONC_LOG_ERROR("unbundle: problem copying {} to {}: {}", originalFilename,
                     newFilename, e.what());
    }
  }

  return filenames;
}

void deleteUnbundledFiles(const std::vector<std::string> &unbundledFilenames) {
  for (const auto &filename : unbundledFilenames) {
    std::error_code ec;
    fs::remove(filename, ec);
    if (ec) {
      EONC_LOG_ERROR("deleteUnbundledFiles: cannot remove {}: {}", filename,
                     ec.message());
    }
  }
}

namespace {

/// Split off the extension the bundle number has to be inserted before.
/// Compression suffixes are carried with the extension they wrap, so
/// "results.con.gz" splits as {"results", ".con.gz"} and the bundled name
/// stays a readable .con.gz rather than "results.con_3.gz".
std::pair<std::string, std::string>
splitBundleExtension(const std::string &name) {
  const size_t pos = name.find_last_of('.');
  if (pos == std::string::npos) {
    return {name, ""};
  }
  const std::string tail = name.substr(pos);
  if (pos > 0 && (tail == ".gz" || tail == ".zst")) {
    const size_t inner = name.find_last_of('.', pos - 1);
    if (inner != std::string::npos) {
      return {name.substr(0, inner), name.substr(inner)};
    }
  }
  return {name.substr(0, pos), tail};
}

} // namespace

void bundle(int number, const std::vector<std::string> &filenames,
            std::vector<std::string> *bundledFilenames) {
  for (const auto &filename : filenames) {
    const auto [baseName, ext] = splitBundleExtension(filename);
    const std::string newFilename =
        baseName + "_" + std::to_string(number) + ext;

    try {
      fs::rename(filename, newFilename);
      bundledFilenames->push_back(newFilename);
    } catch (const fs::filesystem_error &e) {
      EONC_LOG_ERROR("bundle: cannot rename {} to {}: {}", filename,
                     newFilename, e.what());
    }
  }
}

} // namespace eonc
