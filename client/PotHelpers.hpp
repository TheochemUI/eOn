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
#pragma once

#include "C_Structs.h"
#include <algorithm>
#include <numeric>
#include <vector>

namespace eonc::pot {
// Typically this is done by the caller, however it is here as a sanity check
void zeroForceOut(const size_t &nAtoms, ForceOut *efvd);

template <typename T> class registry {
public:
  std::vector<size_t> instanceForceCalls;
  std::vector<bool> instanceActive;

  registry() { addInstance(); }

  void incrementForceCalls() {
    for (size_t i = 0; i < instanceActive.size(); ++i) {
      if (instanceActive[i]) {
        instanceForceCalls[i]++;
        break;
      }
    }
  }

  size_t getTotalForceCalls() const {
    return std::accumulate(instanceForceCalls.begin(), instanceForceCalls.end(),
                           0);
  }

  size_t getInstances() const {
    return std::count(instanceActive.begin(), instanceActive.end(), true);
  }

  void addInstance() {
    instanceForceCalls.push_back(0);
    instanceActive.push_back(true);
  }

  void removeInstance() {
    for (size_t i = 0; i < instanceActive.size(); ++i) {
      if (instanceActive[i]) {
        instanceActive[i] = false;
        break;
      }
    }
  }
};

} // namespace eonc::pot
