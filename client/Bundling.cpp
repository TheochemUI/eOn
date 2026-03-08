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

#include "Bundling.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
using namespace std;

namespace fs = std::filesystem;

int getBundleSize(void) {
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
    std::string ext = originalFilename.substr(dpos + 1);
    std::string newFilename = baseName + "." + ext;

    try {
      fs::copy_file(originalFilename, newFilename,
                    fs::copy_options::overwrite_existing);
    } catch (const fs::filesystem_error &e) {
      std::cerr << "error: unbundle: problem copying " << originalFilename
                << " to " << newFilename << ": " << e.what() << '\n';
    }

    filenames.push_back(newFilename);
  }

  return filenames;
}

void deleteUnbundledFiles(const std::vector<std::string> &unbundledFilenames) {
  for (const auto &filename : unbundledFilenames) {
    fs::remove(filename);
  }
}

void bundle(int number, const std::vector<std::string> &filenames,
            std::vector<std::string> *bundledFilenames) {
  for (const auto &filename : filenames) {
    std::string newFilename = filename;

    size_t pos = newFilename.find_last_of('.');
    if (pos != std::string::npos) {
      std::string fileEnding = newFilename.substr(pos + 1);
      newFilename = newFilename.substr(0, pos) + "_" + std::to_string(number) +
                    "." + fileEnding;
    } else {
      newFilename = newFilename + "_" + std::to_string(number);
    }

    try {
      fs::rename(filename, newFilename);
    } catch (const fs::filesystem_error &e) {
      std::cerr << "error: bundle: cannot rename " << filename << " to "
                << newFilename << ": " << e.what() << '\n';
    }

    bundledFilenames->push_back(newFilename);
  }
}
