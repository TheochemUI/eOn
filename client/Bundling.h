#ifndef BUNDLING_H
#define BUNDLING_H
#include <string>
#include <vector>

int getBundleSize(void);
std::vector<std::string> unbundle(int number);
void bundle(int number, const std::vector<std::string> &filenames,
            std::vector<std::string> *bundledFilenames);
void deleteUnbundledFiles(const std::vector<std::string> &unbundledFilenames);
#endif
