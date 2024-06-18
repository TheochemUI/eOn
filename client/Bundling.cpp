
#include "Bundling.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <filesystem>
#include <iostream>
#include <unistd.h>

namespace fs = std::filesystem;

int getBundleSize(void) {
  DIR *dir;
  struct dirent *dp;
  int num_bundle = -1;

  dir = opendir(".");

  if (!dir) {
    perror("opendir");
    return num_bundle;
  }

  // Find the highest numbered bundle file
  // and return that number.
  while ((dp = readdir(dir))) {
    if (dp->d_name[0] == '.') {
      continue;
    }

    // If "config" is not in the filename
    // then skip.
    if (strstr(dp->d_name, "config") == NULL &&
        strstr(dp->d_name, "ini") == NULL) {
      continue;
    }

    // Find the last underscore
    char *ch = strrchr(dp->d_name, '_');
    if (ch == NULL) {
      continue;
    } else {
      ch += 1;
    }
    // Find the last period
    char *cch = strrchr(dp->d_name, '.');
    if (cch == NULL)
      continue;
    *cch = '\0';
    if (isdigit(*ch)) {
      int i = atoi(ch) + 1;
      if (i > num_bundle) {
        num_bundle = i;
      }
    }
  }

  closedir(dir);

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
  DIR *dir;
  struct dirent *dp;

  dir = opendir(".");

  if (!dir) {
    perror("opendir");
    return filenames;
  }

  while ((dp = readdir(dir))) {
    int bundleNumber;
    std::string originalFilename(dp->d_name);

    if (dp->d_name[0] == '.') {
      continue;
    }

    int numUnderscores = strchrcount(dp->d_name, '_');
    if (numUnderscores < 1) {
      continue;
    }

    // Find the last underscore
    char *ch = strrchr(dp->d_name, '_');
    if (ch == NULL) {
      continue;
    } else {
      ch += 1;
    }
    // Find the last period
    char *cch = strrchr(dp->d_name, '.');
    if (cch == NULL)
      continue;
    *cch = '\0';
    if (isdigit(*ch)) {
      bundleNumber = atoi(ch);
      if (bundleNumber != number) {
        continue;
      }
    }

    *(ch - 1) = '\0';

    std::string newFilename =
        std::string(dp->d_name) + "." + std::string(cch + 1);

    try {
      fs::copy_file(originalFilename, newFilename,
                    fs::copy_options::overwrite_existing);
    } catch (const fs::filesystem_error &e) {
      std::cerr << "error: unbundle: problem copying " << originalFilename
                << " to " << newFilename << ": " << e.what() << '\n';
    }

    filenames.push_back(newFilename);
  }

  closedir(dir);
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
