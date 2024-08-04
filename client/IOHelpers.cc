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
#include "IOHelpers.hpp"
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

} // namespace eonc::io
