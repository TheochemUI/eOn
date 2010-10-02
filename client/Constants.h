/*
 *===============================================
 *  EON Constants.h
 *===============================================
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

/* Collection of constants */

#define STRING_SIZE 512
//const long potentialNewPotential = 0;

namespace constants
{
    const std::string READ("rb"); // set the file to read
    const std::string APPEND("ab"); // set the file to append
    const std::string WRITE("wb"); // set the file to write
}
#endif
