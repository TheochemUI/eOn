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

#include "client/io/writers/con/ConWriter.hpp"

namespace eonc::io {

/**
 * @brief Writes the Matter object to a .con file.
 *
 * This function writes the details of the Matter object, such as atomic
 * positions, cell dimensions, and atomic types, to a .con file. The file format
 * is consistent with eOn-generated files.
 *
 * @param mat The Matter object containing atomic data.
 * @param fout The output file stream where the data will be written.
 * @return True if the file is written successfully, false otherwise.
 *
 * @note The function assumes that atoms are grouped by their atomic numbers.
 *       If the Matter object is empty (i.e., contains no atoms), the function
 *       returns false.
 *       The function uses unique masses and atomic numbers while preserving the
 *       order of their first occurrences.
 */
bool ConWriter::writeImpl(const Matter &mat, std::ofstream &fout) {
  return writeBase(mat, fout, false);
}

} // namespace eonc::io
