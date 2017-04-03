
#ifndef BUNDLING_H
#define BUNDLING_H
#include <vector>
#include <string>

int getBundleSize(void);
std::vector<std::string> unbundle(int number);
void bundle(int number, std::vector<std::string> filenames,
std::vector<std::string> *bundledFilenames);
void deleteUnbundledFiles(std::vector<std::string> unbundledFilenames);
#endif
