#ifndef COMPRESSION_H
#define COMPRESSIONH

#include <string>
#include <vector>

int create_archive(char *outname, char *path, 
                   const std::vector<std::string> &filenames);
int extract_archive(char *filename);

#endif
