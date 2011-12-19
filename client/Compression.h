//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef COMPRESSION_H
#define COMPRESSIONH

#include <string>
#include <vector>

int create_archive(char *outname, char *path, 
                   const std::vector<std::string> &filenames);
int extract_archive(char *filename);

#endif
