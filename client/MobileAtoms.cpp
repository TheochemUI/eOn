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
#include "eon/MobileAtoms.h"

#include <cctype>
#include <string>
#include <unordered_set>
#include <vector>

namespace eonc {

bool atomListMeansAll(const std::string &atomList) {
  if (atomList.empty()) {
    return true;
  }
  // Trim and case-fold "all"
  size_t b = 0;
  while (b < atomList.size() &&
         std::isspace(static_cast<unsigned char>(atomList[b]))) {
    ++b;
  }
  size_t e = atomList.size();
  while (e > b && std::isspace(static_cast<unsigned char>(atomList[e - 1]))) {
    --e;
  }
  if (e <= b) {
    return true;
  }
  if (e - b != 3) {
    return false;
  }
  auto lower = [](char c) {
    return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  };
  return lower(atomList[b]) == 'a' && lower(atomList[b + 1]) == 'l' &&
         lower(atomList[b + 2]) == 'l';
}

VectorXi freeAtomIndices(const Matter *matter) {
  const long n = matter->numberOfAtoms();
  std::vector<int> free;
  free.reserve(static_cast<size_t>(matter->numberOfFreeAtoms()));
  for (long i = 0; i < n; ++i) {
    if (!matter->getFixed(i)) {
      free.push_back(static_cast<int>(i));
    }
  }
  VectorXi out(static_cast<Eigen::Index>(free.size()));
  for (Eigen::Index k = 0; k < out.size(); ++k) {
    out(k) = free[static_cast<size_t>(k)];
  }
  return out;
}

VectorXi resolveMobileAtoms(const Matter *matter, const std::string &atomList) {
  if (atomListMeansAll(atomList)) {
    return freeAtomIndices(matter);
  }
  const long n = matter->numberOfAtoms();
  std::vector<int> mobile;
  std::unordered_set<int> seen;
  std::string token;
  for (size_t p = 0; p <= atomList.size(); ++p) {
    const char c = (p < atomList.size()) ? atomList[p] : ',';
    if (c == ',' || c == ' ' || c == '\t' || p == atomList.size()) {
      if (!token.empty()) {
        try {
          const long idx = std::stol(token);
          if (idx >= 0 && idx < n && !matter->getFixed(idx) &&
              seen.insert(static_cast<int>(idx)).second) {
            mobile.push_back(static_cast<int>(idx));
          }
        } catch (const std::exception &) {
          // skip non-integer tokens
        }
        token.clear();
      }
    } else {
      token.push_back(c);
    }
  }
  VectorXi out(static_cast<Eigen::Index>(mobile.size()));
  for (Eigen::Index k = 0; k < out.size(); ++k) {
    out(k) = mobile[static_cast<size_t>(k)];
  }
  return out;
}

VectorXi resolveMobileAtoms(const Matter *matter, const VectorXi &candidates) {
  const long n = matter->numberOfAtoms();
  std::vector<int> mobile;
  std::unordered_set<int> seen;
  mobile.reserve(static_cast<size_t>(candidates.size()));
  for (Eigen::Index k = 0; k < candidates.size(); ++k) {
    const long idx = candidates(k);
    if (idx >= 0 && idx < n && !matter->getFixed(idx) &&
        seen.insert(static_cast<int>(idx)).second) {
      mobile.push_back(static_cast<int>(idx));
    }
  }
  VectorXi out(static_cast<Eigen::Index>(mobile.size()));
  for (Eigen::Index k = 0; k < out.size(); ++k) {
    out(k) = mobile[static_cast<size_t>(k)];
  }
  return out;
}

VectorXd packMobileRows(const AtomMatrix &full, const VectorXi &mobile) {
  VectorXd packed(3 * mobile.size());
  for (Eigen::Index a = 0; a < mobile.size(); ++a) {
    packed.segment<3>(3 * a) = full.row(mobile(a));
  }
  return packed;
}

void unpackMobileRows(const VectorXd &packed, const VectorXi &mobile,
                      AtomMatrix &full) {
  for (Eigen::Index a = 0; a < mobile.size(); ++a) {
    full.row(mobile(a)) = packed.segment<3>(3 * a);
  }
}

VectorXd mobileForces(Matter *matter, const VectorXi &mobile) {
  const AtomMatrix forces = matter->getForces();
  return packMobileRows(forces, mobile);
}

} // namespace eonc
