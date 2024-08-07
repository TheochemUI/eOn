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
#include "client/io/IOHelpers.hpp"
#include <fmt/core.h>
#include <fstream>
#include <map>

namespace eonc::io {

/**
 * @brief Extracts unique elements from the input VectorType while
 * preserving their insertion order.
 *
 * This function takes an VectorType as input and returns a new vector
 * containing only the unique elements from the input vector, while maintaining
 * their original insertion order. It leverages a std::map to ensure uniqueness
 * and uses the insertion order as a value to preserve the sequence.
 *
 * @param vec The input VectorType containing elements that may have
 * duplicates.
 * @return VectorType An VectorType containing unique elements in their
 * original insertion order.
 *
 * @note The function uses a std::map to store each unique element and its
 * insertion order. It then sorts the elements based on their order of
 * appearance in the input vector to ensure that the resultant vector preserves
 * the original sequence of first occurrences.
 *
 * @todo(rg) The map may be a bit slow [1], and ideally this might be something
 * for Matter to have setup instead of traversing for every write out [1]
 * https://stackoverflow.com/a/30153545/1895378
 *
 * @example
 *   VectorType vec(10);
 *   vec << 1, 2, 3, 2, 1, 4, 5, 4, 6, 5;
 *   VectorType uniqueVec = getUniqueElements(vec);
 *   // uniqueVec now contains: 1, 2, 3, 4, 5, 6
 */
VectorType getUniqueValues(const VectorType &vec) {
  // Map to store unique elements and their insertion order
  std::map<double, size_t> elementOrder;
  size_t order = 0;

  // Iterate through the input vector and store each unique element with its
  // insertion order
  for (auto i = 0; i < vec.size(); ++i) {
    if (elementOrder.find(vec(i)) == elementOrder.end()) {
      elementOrder[vec(i)] = order++;
    }
  }

  // Convert the map to a vector of pairs to sort by insertion order
  std::vector<std::pair<double, size_t>> orderedElements(elementOrder.begin(),
                                                         elementOrder.end());
  std::sort(orderedElements.begin(), orderedElements.end(),
            [](const auto &a, const auto &b) { return a.second < b.second; });

  // Create a new VectorType to hold the unique elements in their original
  // order
  VectorType uniqueVec(orderedElements.size());
  for (size_t i = 0; i < orderedElements.size(); ++i) {
    uniqueVec(i) = orderedElements[i].first;
  }

  return uniqueVec;
}

Vector<size_t> getUniqueCounts(const Vector<size_t> &vec) {
  // Map to store unique elements and their counts
  std::map<double, size_t> elementCounts;

  // Iterate through the input vector and count each unique element
  for (auto i = 0; i < vec.size(); ++i) {
    elementCounts[vec(i)]++;
  }

  // Create a new variable to hold the counts in the original order of their
  // first appearance
  Vector<size_t> uniqueCounts(elementCounts.size());
  size_t index = 0;

  for (auto i = 0; i < vec.size(); ++i) {
    if (elementCounts.find(vec(i)) != elementCounts.end()) {
      uniqueCounts(index++) = elementCounts[vec(i)];
      // Remove the element to ensure it's counted only once
      elementCounts.erase(vec(i));
    }
  }

  return uniqueCounts;
}

/**
 * @brief Ensures that the output file stream is open.
 *
 * This function attempts to open the specified output file stream.
 * If the `append` flag is true, it opens the file in append mode.
 * If the file does not exist, it throws an exception.
 * If the `append` flag is false, it opens the file in truncate mode.
 *
 * @param file The output file stream to be opened.
 * @param filename The path of the file to be opened.
 * @param append Flag indicating whether to open the file in append mode.
 *
 * @throws std::runtime_error If the file does not exist when appending.
 */
void ensureFileOpen(std::ofstream &file, const std::filesystem::path &filename,
                    bool append) {
  if (append) {
    file.open(filename, std::ios::out | std::ios::app);
    if (!file.is_open()) {
      throw std::runtime_error(
          fmt::format("Cannot open file {} for appending", filename.c_str()));
    }
  } else {
    file.open(filename, std::ios::out | std::ios::trunc);
    if (!file.is_open()) {
      throw std::runtime_error(
          fmt::format("Cannot open file {} for writing", filename.c_str()));
    }
  }
}

/**
 * @brief Ensures that the input file stream is open.
 *
 * This function attempts to open the specified input file stream.
 * If the file does not exist or cannot be opened, it throws an exception.
 *
 * @param file The input file stream to be opened.
 * @param filename The path of the file to be opened.
 *
 * @throws std::runtime_error If the file does not exist or cannot be opened.
 */
void ensureFileOpen(std::ifstream &file,
                    const std::filesystem::path &filename) {
  file.open(filename, std::ios::in);
  if (!file.is_open()) {
    throw std::runtime_error(
        fmt::format("Cannot open file {}", filename.c_str()));
  }
}

} // namespace eonc::io