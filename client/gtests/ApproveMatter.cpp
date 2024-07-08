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
#include "ApprovalTests.hpp"
#include "Matter.h"
#include "catch2/catch_amalgamated.hpp"
#include <iomanip>
#include <iostream>
#include <string>

std::string toSpaceSeparatedString(const Eigen::VectorXd &vec) {
  std::ostringstream oss;
  for (int i = 0; i < vec.size(); ++i) {
    oss << vec[i];
    if (i < vec.size() - 1) {
      oss << " ";
    }
  }
  return oss.str();
}

double perAtomNormConst(const Matter &matter) {
  return const_cast<Matter &>(matter).perAtomNorm(matter);
}

double distanceToConst(const Matter &matter) {
  return const_cast<Matter &>(matter).distanceTo(matter);
}

double getPotentialEnergyConst(const Matter &matter) {
  return const_cast<Matter &>(matter).getPotentialEnergy();
}

double maxForceConst(const Matter &matter) {
  return const_cast<Matter &>(matter).maxForce();
}

auto getForcesConst(const Matter &matter) {
  return const_cast<Matter &>(matter).getForces();
}

auto getForcesVConst(const Matter &matter) {
  return const_cast<Matter &>(matter).getForcesV();
}

auto getEnergyVarianceConst(const Matter &matter) {
  return const_cast<Matter &>(matter).getEnergyVariance();
}

auto getMechanicalEnergyConst(const Matter &matter) {
  return const_cast<Matter &>(matter).getMechanicalEnergy();
}

auto getForcesFreeConst(const Matter &matter) {
  return const_cast<Matter &>(matter).getForcesFree();
}

auto getForcesFreeVConst(const Matter &matter) {
  return const_cast<Matter &>(matter).getForcesFreeV();
}

auto getBiasForcesConst(const Matter &matter) {
  return const_cast<Matter &>(matter).getBiasForces();
}

std::ostream &operator<<(std::ostream &os, const Matter &matter) {
  os << std::setprecision(8);
  os << "Matter has: " << std::endl;

  // Methods without parameters
  os << "Number of atoms: \n" << matter.numberOfAtoms() << std::endl;
  os << "Cell: \n" << matter.getCell() << std::endl;
  os << "Potential calls: \n" << matter.getPotentialCalls() << std::endl;
  os << "Positions: \n" << matter.getPositions() << std::endl;
  os << "Positions (Vector): \n"
     << toSpaceSeparatedString(matter.getPositionsV()) << std::endl;
  os << "Positions (Free): \n" << matter.getPositionsFree() << std::endl;
  os << "Positions (Free Vector): \n"
     << toSpaceSeparatedString(matter.getPositionsFreeV()) << std::endl;
  os << "Velocities: \n" << matter.getVelocities() << std::endl;
  // TODO(rg): These have no initializers
  // os << "Potential: \n" << (matter.getPotential() ? "Exists" : \n"None") <<
  // std::endl; os << "Bias Forces: \n" << getBiasForcesConst(matter) <<
  // std::endl;
  os << "Forces: \n" << getForcesConst(matter) << std::endl;
  os << "Forces (Vector): \n"
     << toSpaceSeparatedString(getForcesVConst(matter)) << std::endl;
  os << "Forces (Free): \n" << getForcesFreeConst(matter) << std::endl;
  os << "Forces (Free Vector): \n"
     << toSpaceSeparatedString(getForcesFreeVConst(matter)) << std::endl;
  os << "Masses: \n" << matter.getMasses() << std::endl;
  os << "Atomic Numbers: \n" << matter.getAtomicNrs() << std::endl;
  os << "Fixed Atoms: \n" << matter.getFixed(3) << std::endl;
  os << "Energy Variance: \n" << getEnergyVarianceConst(matter) << std::endl;
  os << "Potential Energy: \n" << getPotentialEnergyConst(matter) << std::endl;
  os << "Kinetic Energy: \n" << matter.getKineticEnergy() << std::endl;
  os << "Mechanical Energy: \n"
     << getMechanicalEnergyConst(matter) << std::endl;
  os << "Number of Free Atoms: \n" << matter.numberOfFreeAtoms() << std::endl;
  os << "Number of Fixed Atoms: \n" << matter.numberOfFixedAtoms() << std::endl;
  // TODO(rg): This is weirdly broken
  // os << "Force Calls: \n" << matter.getForceCalls() << std::endl;
  os << "Maximum Force: \n" << maxForceConst(matter) << std::endl;
  os << "Free (Matrix): \n" << matter.getFree() << std::endl;
  os << "Free (Vector): \n"
     << toSpaceSeparatedString(matter.getFreeV()) << std::endl;

  // Methods with parameters
  os << "Distance to self: \n" << distanceToConst(matter) << std::endl;
  os << "Per atom norm with self: \n" << perAtomNormConst(matter) << std::endl;
  os << "Position of atom 0, axis 0: \n"
     << matter.getPosition(0, 0) << std::endl;
  os << "Mass of atom 0: \n" << matter.getMass(0) << std::endl;
  os << "Atomic number of atom 0: \n" << matter.getAtomicNr(0) << std::endl;
  os << "Is atom 0 fixed: \n" << matter.getFixed(0) << std::endl;
  os << "Distance between atoms 0 and 1: \n"
     << matter.distance(0, 1) << std::endl;
  os << "P-distance between atoms 0 and 1, axis 0: \n"
     << matter.pdistance(0, 1, 0) << std::endl;
  os << "Distance between same atom in two configurations: \n"
     << matter.distance(matter, 0) << std::endl;

  AtomMatrix placeholderMatrix = Eigen::MatrixXd::Zero(0, 3);
  Eigen::VectorXd placeholderVector = Eigen::VectorXd::Zero(0);
  os << "PBC Matrix: \n" << matter.pbc(placeholderMatrix) << std::endl;
  os << "PBC Vector: \n" << matter.pbcV(placeholderVector) << std::endl;

  return os;
}

std::vector<Matter> getMatter() {
  // Return test data for Matter
  // TODO(rg): Add more objects
  auto params = std::make_shared<Parameters>();
  auto pot_default = helper_functions::makePotential(PotType::LJ, params);
  auto m1 = Matter(pot_default, params);
  std::string confile("pos.con"); // Sulfolene
  m1.con2matter(confile);
  for (size_t idx{0}; idx < 3; idx++) {
    m1.setFixed(idx, true);
  }
  return {m1};
}

TEST_CASE("VerifyMatter") {
  ApprovalTests::Approvals::verifyAll("matter", getMatter());
}
