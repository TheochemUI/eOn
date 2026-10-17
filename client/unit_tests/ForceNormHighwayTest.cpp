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

#include "eon/ForceNorm.h"
#include "eon/PbcSimd.h"

#include "catch2/catch_amalgamated.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

#ifdef WITH_HIGHWAY
#include "hwy/targets.h"
#endif

namespace {

uint64_t doubleBits(double value) {
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof bits);
  return bits;
}

void requireSameBits(double got, double ref) {
  const uint64_t gotBits = doubleBits(got);
  const uint64_t refBits = doubleBits(ref);
  INFO("got " << got << " bits " << std::hex << gotBits);
  INFO("ref " << ref << " bits " << std::hex << refBits);
  REQUIRE(gotBits == refBits);
}

#ifdef WITH_HIGHWAY
struct HwyTargetReset {
  HwyTargetReset() = default;
  ~HwyTargetReset() { hwy::SetSupportedTargetsForTest(0); }
};
#endif

// Runs fn once per compiled target the CPU can execute, then clears the
// Highway test mask. Without Highway there is only the scalar path.
template <typename Fn> void forEachTarget(Fn &&fn) {
#ifdef WITH_HIGHWAY
  [[maybe_unused]] const HwyTargetReset reset;
  const std::vector<int64_t> targets = hwy::SupportedAndGeneratedTargets();
  std::string ran;
  for (int64_t target : targets) {
    if (!ran.empty()) {
      ran.push_back(' ');
    }
    ran += hwy::TargetName(target);
  }
  std::cout << "Highway targets compared: " << ran << '\n';
  REQUIRE_FALSE(targets.empty());
  for (int64_t target : targets) {
    hwy::SetSupportedTargetsForTest(target);
    fn();
  }
#else
  std::cout << "Highway targets compared: scalar\n";
  fn();
#endif
}

double drawMantissa(std::mt19937_64 &rng) {
  const uint64_t mantissa = rng() & 0x000fffffffffffffULL;
  const uint64_t exponent = 1000 + (rng() % 40);
  const uint64_t sign = (rng() & 1ULL) << 63;
  const uint64_t bits = sign | (exponent << 52) | mantissa;
  double value = 0.0;
  std::memcpy(&value, &bits, sizeof value);
  return value;
}

struct ForceCase {
  std::string name;
  std::vector<double> forces;
  std::vector<double> fixed;
  long nAtoms = 0;
};

void checkForceCase(const ForceCase &item) {
  const double *forces = item.forces.empty() ? nullptr : item.forces.data();
  const double *fixed = item.fixed.empty() ? nullptr : item.fixed.data();
  const double ref =
      eonc::detail::maxFreeAtomForceNormScalar(forces, fixed, 0, item.nAtoms);
  const double got = eonc::maxFreeAtomForceNorm(forces, fixed, item.nAtoms);
  INFO(item.name);
  requireSameBits(got, ref);
}

// One atom whose unfused sum of squares differs in the last bit from a
// contracted fma under GCC 13 -ffp-contract=fast -mavx2.
ForceCase fusedDiffers(long nAtoms, long at) {
  ForceCase item;
  item.name =
      "fma_last_bit n=" + std::to_string(nAtoms) + " at=" + std::to_string(at);
  item.nAtoms = nAtoms;
  item.forces.assign(static_cast<size_t>(3 * nAtoms), 0.0);
  const size_t base = static_cast<size_t>(3 * at);
  item.forces[base] = -0x1.6f6c93e34b87ep+12;
  item.forces[base + 1] = -0x1.0661fa48720cfp-10;
  item.forces[base + 2] = 0x1.8d13a5e25d153p-8;
  return item;
}

} // namespace

TEST_CASE("MaxFreeAtomForceNorm matches the scalar bits on every target",
          "[highway][fp-contract]") {
  std::vector<ForceCase> cases;
  const long lengths[] = {1, 3, 7, 8, 9, 63, 64, 65, 1000};
  std::mt19937_64 rng(0xC0FFEEULL);
  for (long n : lengths) {
    ForceCase item;
    item.name = "random_mantissa";
    item.nAtoms = n;
    item.forces.resize(static_cast<size_t>(3 * n));
    item.fixed.assign(item.forces.size(), 0.0);
    for (double &component : item.forces) {
      component = drawMantissa(rng);
    }
    for (long i = 0; i < n; ++i) {
      if (i % 5 == 0) {
        item.fixed[static_cast<size_t>(3 * i)] = 1.0;
        item.fixed[static_cast<size_t>(3 * i + 1)] = 1.0;
        item.fixed[static_cast<size_t>(3 * i + 2)] = 1.0;
      } else if (i % 5 == 1) {
        item.fixed[static_cast<size_t>(3 * i)] = 1.0;
      }
    }
    cases.push_back(std::move(item));
  }

  cases.push_back(fusedDiffers(1, 0));
  cases.push_back(fusedDiffers(8, 0));
  cases.push_back(fusedDiffers(9, 8));

  {
    ForceCase item;
    item.name = "minus_zero";
    item.nAtoms = 4;
    item.forces = {-0.0, 0.0,  0.0, 3.0, 4.0, 0.0,
                   0.0,  -0.0, 1.0, 1.0, 1.0, 1.0};
    cases.push_back(std::move(item));
  }
  {
    ForceCase item;
    item.name = "subnormal";
    item.nAtoms = 3;
    const double tiny = std::numeric_limits<double>::denorm_min();
    item.forces = {tiny, 0.0, -0.0, 0.5, tiny, 0.25, 0.0, 0.0, tiny};
    cases.push_back(std::move(item));
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  {
    ForceCase item;
    item.name = "nan_only";
    item.nAtoms = 1;
    item.forces = {0.25, nan, -1.0};
    cases.push_back(std::move(item));
  }
  {
    ForceCase item;
    item.name = "nan_in_vector";
    item.nAtoms = 8;
    item.forces.assign(24, 0.25);
    item.forces[1] = nan;
    cases.push_back(std::move(item));
  }
  {
    ForceCase item;
    item.name = "nan_in_tail";
    item.nAtoms = 9;
    item.forces.assign(27, 0.25);
    item.forces[25] = nan;
    cases.push_back(std::move(item));
  }
  {
    ForceCase item;
    item.name = "nan_on_fixed_atom";
    item.nAtoms = 8;
    item.forces.assign(24, 0.25);
    item.forces[1] = nan;
    item.fixed.assign(24, 0.0);
    item.fixed[0] = 1.0;
    item.fixed[1] = 1.0;
    item.fixed[2] = 1.0;
    cases.push_back(std::move(item));
  }

  forEachTarget([&]() {
    for (const ForceCase &item : cases) {
      checkForceCase(item);
    }
  });
}

TEST_CASE("PBC wraps match the scalar formulas on every target",
          "[highway][fp-contract]") {
  // The scalar helpers in PbcSimd.cpp are file-local. These are the same
  // formulas the header documents: x - floor(x + 0.5), and
  // (x + 1) - trunc(x + 1).
  std::vector<std::vector<double>> inputs;
  const size_t lengths[] = {0, 1, 3, 7, 8, 9, 63, 64, 65};
  std::mt19937_64 rng(0xA11CEULL);
  for (size_t n : lengths) {
    std::vector<double> values(n);
    for (double &value : values) {
      value = drawMantissa(rng);
    }
    inputs.push_back(std::move(values));
  }
  inputs.push_back({-0.0, 0.0, 0.5, -0.5, 1.0, -1.0, 1.5, -1.25});
  inputs.push_back(
      {std::numeric_limits<double>::denorm_min(), -0x1.8p-1, 0x1.0p+2});

  forEachTarget([&]() {
    for (const std::vector<double> &input : inputs) {
      std::vector<double> gotMin = input;
      std::vector<double> refMin = input;
      if (!gotMin.empty()) {
        eonc::pbc::wrapMinimumImage(gotMin.data(), gotMin.size());
      }
      for (size_t i = 0; i < refMin.size(); ++i) {
        refMin[i] = refMin[i] - std::floor(refMin[i] + 0.5);
      }
      REQUIRE(gotMin.size() == refMin.size());
      for (size_t i = 0; i < gotMin.size(); ++i) {
        requireSameBits(gotMin[i], refMin[i]);
      }

      std::vector<double> gotLegacy = input;
      std::vector<double> refLegacy = input;
      if (!gotLegacy.empty()) {
        eonc::pbc::wrapLegacyUnit(gotLegacy.data(), gotLegacy.size());
      }
      for (size_t i = 0; i < refLegacy.size(); ++i) {
        const double shifted = refLegacy[i] + 1.0;
        refLegacy[i] = shifted - std::trunc(shifted);
      }
      for (size_t i = 0; i < gotLegacy.size(); ++i) {
        requireSameBits(gotLegacy[i], refLegacy[i]);
      }
    }
  });
}
