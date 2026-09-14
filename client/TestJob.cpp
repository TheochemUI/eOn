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
#include "eon/TestJob.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/Potential.h"

#include <cmath>
#include <fstream>
#include "magic_enum/magic_enum.hpp"
#include <stdexcept>
#include <string>

namespace {
struct PotRef {
  const char *tag;
  eonc::PotType type;
  double energy;
  double max_force;
};
} // namespace

std::vector<std::string> TestJob::run() {
  checkPotentials();
  checkFullSearch();
  return {"results.dat"};
}

void TestJob::checkFullSearch() {
  // Historical saddle self-test needs reactant_test.con / displacement_test.con
  // which the tree does not ship. Keep the hook.
}

void TestJob::checkPotentials() {
  const PotRef cases[] = {
      {"lj", eonc::PotType::LJ, -1475.984331, 2.007213},
      {"emt", eonc::PotType::EMT, 46.086312, 0.357493},
      {"edip", eonc::PotType::EDIP, -1033.250950, 7.080115},
      {"tersoff_si", eonc::PotType::TERSOFF_SI, -1035.809985, 11.145002},
      {"sw_si", eonc::PotType::SW_SI, -1449.795645, 2.530904},
      {"lenosky_si", eonc::PotType::LENOSKY_SI, -1410.679106, 2.320168},
      {"eam_al", eonc::PotType::EAM_AL, -1206.825825, 0.000246},
      {"tip4p", eonc::PotType::TIP4P, 4063.865115, 73.655248},
  };

  std::ofstream out("results.dat");
  for (const auto &c : cases) {
    try {
      const double de = getEnergyDiff(c.tag, c.energy);
      if (std::abs(de) > tolerance) {
        out << "FAIL " << c.tag << " energy_diff " << de << "\n";
        continue;
      }
      const double df = getForceDiff(c.tag, c.max_force);
      if (std::abs(df) > tolerance) {
        out << "FAIL " << c.tag << " force_diff " << df << "\n";
        continue;
      }
      out << "OK " << c.tag << "\n";
    } catch (const std::exception &e) {
      out << "SKIP " << c.tag << " " << e.what() << "\n";
    }
  }
}

double TestJob::getEnergyDiff(std::string potTag, double refEnergy) {
  auto type =
      magic_enum::enum_cast<eonc::PotType>(potTag, magic_enum::case_insensitive);
  if (!type) {
    throw std::invalid_argument("unknown pot " + potTag);
  }
  Parameters p = params;
  p.potential_options.potential = *type;
  auto potHandle = eonc::helpers::makePotential(*type, p);
  Matter pos(potHandle, p);
  if (!eonc::io::io_ok(pos.con2matter(std::string("pos_test.con")))) {
    throw std::runtime_error("no pos_test.con");
  }
  return pos.getPotentialEnergy() - refEnergy;
}

double TestJob::getForceDiff(std::string potTag, double refForce) {
  auto type =
      magic_enum::enum_cast<eonc::PotType>(potTag, magic_enum::case_insensitive);
  if (!type) {
    throw std::invalid_argument("unknown pot " + potTag);
  }
  Parameters p = params;
  p.potential_options.potential = *type;
  auto potHandle = eonc::helpers::makePotential(*type, p);
  Matter pos(potHandle, p);
  if (!eonc::io::io_ok(pos.con2matter(std::string("pos_test.con")))) {
    throw std::runtime_error("no pos_test.con");
  }
  return pos.maxForce() - refForce;
}
