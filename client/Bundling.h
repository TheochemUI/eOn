//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

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
