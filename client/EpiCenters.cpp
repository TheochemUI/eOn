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
#include "eon/EpiCenters.h"
#include "eon/HelperFunctions.h"

#include <climits>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace eonc::helpers;

namespace {
long pickFromHits(const std::vector<long> &hits, const char *what) {
  if (hits.empty()) {
    throw std::runtime_error(what);
  }
  const long pick =
      static_cast<long>(randomDouble(static_cast<long>(hits.size())));
  return hits[static_cast<size_t>(pick)];
}
} // namespace

long eonc::EpiCenters::cnaEpiCenter(const Matter *matter,
                                    double neighborCutoff) {
  if (!matter) {
    throw std::invalid_argument("EpiCenters: null Matter");
  }
  const long nAtoms = matter->numberOfAtoms();
  std::vector<long> cnaList(static_cast<size_t>(nAtoms));
  cna(cnaList.data(), matter, neighborCutoff);

  std::vector<long> hits;
  hits.reserve(static_cast<size_t>(nAtoms));
  for (long i = 0; i < nAtoms; ++i) {
    if (cnaList[static_cast<size_t>(i)] == 2 && !matter->getFixed(i)) {
      hits.push_back(i);
    }
  }
  return pickFromHits(hits, "EpiCenters: no free non-FCC/HCP atom");
}

long eonc::EpiCenters::minCoordinatedEpiCenter(const Matter *matter,
                                               double neighborCutoff) {
  if (!matter) {
    throw std::invalid_argument("EpiCenters: null Matter");
  }
  const long nAtoms = matter->numberOfAtoms();
  auto minCoordinatedList = std::make_unique<bool[]>(static_cast<size_t>(nAtoms));
  const long minCoordinationVal = minCoordination(matter, neighborCutoff);
  coordinationLessOrEqual(minCoordinatedList.get(), minCoordinationVal, matter,
                          neighborCutoff);

  std::vector<long> hits;
  hits.reserve(static_cast<size_t>(nAtoms));
  for (long i = 0; i < nAtoms; ++i) {
    if (minCoordinatedList[static_cast<size_t>(i)] && !matter->getFixed(i)) {
      hits.push_back(i);
    }
  }
  return pickFromHits(hits, "EpiCenters: no free minimally coordinated atom");
}

long eonc::EpiCenters::lastAtom(const Matter *matter) {
  if (!matter) {
    throw std::invalid_argument("EpiCenters: null Matter");
  }
  const long nAtoms = matter->numberOfAtoms();
  if (nAtoms <= 0) {
    throw std::runtime_error("EpiCenters: lastAtom on empty Matter");
  }
  return nAtoms - 1;
}

long eonc::EpiCenters::randomFreeAtomEpiCenter(const Matter *matter) {
  if (!matter) {
    throw std::invalid_argument("EpiCenters: null Matter");
  }
  const long nAtoms = matter->numberOfAtoms();
  std::vector<long> hits;
  hits.reserve(static_cast<size_t>(nAtoms));
  for (long i = 0; i < nAtoms; ++i) {
    if (!matter->getFixed(i)) {
      hits.push_back(i);
    }
  }
  return pickFromHits(hits, "EpiCenters: no free atom");
}

void eonc::EpiCenters::cna(long *cna, const Matter *matter,
                           double neighborCutoff) {
  long nAtoms = matter->numberOfAtoms();
  std::vector<int> nFCC(nAtoms, 0);
  std::vector<int> nHCP(nAtoms, 0);
  std::vector<std::vector<int>> neighborLists(nAtoms);

  for (long i = 0; i < nAtoms - 1; i++) {
    for (long j = i + 1; j < nAtoms; j++) {
      double diffR = matter->distance(i, j);
      if (diffR < neighborCutoff) {
        neighborLists[i].push_back(static_cast<int>(j));
        neighborLists[j].push_back(static_cast<int>(i));
      }
    }
  }

  for (long a2 = 0; a2 < nAtoms; a2++) {
    const auto &nbs2 = neighborLists[a2];
    for (size_t n2 = 0; n2 < nbs2.size(); n2++) {
      int a1 = nbs2[n2];
      if (a1 < a2) {
        std::vector<int> common;
        const auto &nbs1 = neighborLists[a1];
        for (size_t n1 = 0; n1 < nbs1.size(); n1++) {
          int a3 = nbs1[n1];
          for (size_t m2 = 0; m2 < nbs2.size(); m2++) {
            if (a3 == nbs2[m2])
              common.push_back(a3);
          }
        }
        if (common.size() == 4) {
          int nBonds = 0;
          int bondsSum = 0;
          for (int j2 = 1; j2 < 4; j2++) {
            const auto &nbs = neighborLists[common[j2]];
            for (int j1 = 0; j1 < j2; j1++) {
              for (size_t n = 0; n < nbs.size(); n++) {
                if (common[j1] == nbs[n]) {
                  nBonds++;
                  bondsSum += j1 + j2;
                  break;
                }
              }
            }
          }
          if (nBonds == 2) {
            if (bondsSum == 6) {
              nFCC[a1]++;
              nFCC[a2]++;
            } else {
              nHCP[a1]++;
              nHCP[a2]++;
            }
          }
        }
      }
    }
  }
  // 0: fcc (421), 1: hcp (422), 2: other
  for (long i = 0; i < nAtoms; i++) {
    if (neighborLists[i].size() == 12) {
      if (nFCC[i] == 12)
        cna[i] = 0;
      else if (nFCC[i] == 6 && nHCP[i] == 6)
        cna[i] = 1;
      else
        cna[i] = 2;
    } else {
      cna[i] = 2;
    }
  }
}

void eonc::EpiCenters::coordination(long *coordinationVal, const Matter *matter,
                                    double neighborCutoff) {
  long nAtoms = matter->numberOfAtoms();
  for (long i = 0; i < nAtoms; i++)
    coordinationVal[i] = 0;

  for (long i = 0; i < nAtoms - 1; i++) {
    for (long j = i + 1; j < nAtoms; j++) {
      double diffR = matter->distance(i, j);
      if (diffR < neighborCutoff) {
        coordinationVal[i]++;
        coordinationVal[j]++;
      }
    }
  }
}

void eonc::EpiCenters::coordinationLessOrEqual(bool *result,
                                               long coordinationMaxVal,
                                               const Matter *matter,
                                               double neighborCutoff) {
  long nAtoms = matter->numberOfAtoms();
  std::vector<long> coordinationVal(nAtoms);

  coordination(coordinationVal.data(), matter, neighborCutoff);

  for (long i = 0; i < nAtoms; i++) {
    result[i] =
        (coordinationVal[i] <= coordinationMaxVal) && !matter->getFixed(i);
  }
}

long eonc::EpiCenters::listedAtomEpiCenter(const Matter *matter,
                                           const std::vector<long> &atomList) {
  if (!matter) {
    throw std::invalid_argument("EpiCenters: null Matter");
  }
  long nAtoms = matter->numberOfAtoms();
  std::vector<long> freeAtoms;
  // Lone -1 is every free atom (akmc-al / ListedAtoms). A mixed list
  // still treats negatives as out of range.
  if (atomList.size() == 1 && atomList[0] == -1) {
    for (long i = 0; i < nAtoms; ++i) {
      if (!matter->getFixed(i)) {
        freeAtoms.push_back(i);
      }
    }
  } else {
    for (long idx : atomList) {
      if (idx < 0 || idx >= nAtoms) {
        continue;
      }
      const long row = matter->mapFileRow(idx);
      if (row >= 0 && row < nAtoms && !matter->getFixed(row)) {
        freeAtoms.push_back(row);
      }
    }
  }
  if (freeAtoms.empty()) {
    throw std::runtime_error("Listed atoms are all frozen");
  }
  // randomDouble(N) is [0, N); size-1 dropped the last listed / last free atom.
  long pick =
      static_cast<long>(randomDouble(static_cast<long>(freeAtoms.size())));
  return freeAtoms[pick];
}

long eonc::EpiCenters::minCoordination(const Matter *matter,
                                       double neighborCutoff) {
  long nAtoms = matter->numberOfAtoms();
  std::vector<long> coordinationVal(nAtoms);

  coordination(coordinationVal.data(), matter, neighborCutoff);

  long minVal = LONG_MAX;
  for (long i = 0; i < nAtoms; i++) {
    if (coordinationVal[i] < minVal && !matter->getFixed(i)) {
      minVal = coordinationVal[i];
    }
  }
  return minVal;
}
