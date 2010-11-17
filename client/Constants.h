//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
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
