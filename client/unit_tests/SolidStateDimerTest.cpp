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
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/HelperFunctions.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/SetflEam.h"
#include "eon/SolidState.h"
#include "eon/SolidStateNEB.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace tests {
namespace {

static eonc::helpers::test::QuillTestLogger quillSetup;

// q = 10 (s_x - 1/2) sits in the middle of a 10 Angstrom cell, so a
// finite-difference step does not wrap through the periodic boundary.
// E = (q^2-1)^2 + 2 (h_yx/10 - q/2)^2 plus stiff penalties. Variable-cell
// minima are (q, h_yx) = (±1, ±5), E = 0. The saddle is q = 0, h_yx = 0,
// E = 1. A frozen cell has minima at q = ±sqrt(3)/2, E = 0.4375, so the
// barrier is 0.5625 instead of 1.
class ShearPot final : public eonc::Potential {
public:
  ShearPot() : Potential(eonc::PotType::LJ) {}

  [[nodiscard]] bool computesStress() const noexcept override { return true; }
  [[nodiscard]] Matrix3d cauchyStress() const override { return sigma_; }
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  void force(long N, const double *R, const int *, double *F, double *U,
             double *, const double *box) override {
    REQUIRE(N == 1);
    Matrix3d h = Matrix3d::Zero();
    for (int i = 0; i < 9; ++i) {
      h(i / 3, i % 3) = box[i];
    }
    const Eigen::RowVector3d r(R[0], R[1], R[2]);
    double energy = 0.0;
    Eigen::RowVector3d force = Eigen::RowVector3d::Zero();
    energyOf(r, h, &energy, &force, nullptr);
    *U = energy;
    F[0] = force[0];
    F[1] = force[1];
    F[2] = force[2];
    const double volume = std::abs(h.determinant());
    sigma_.setZero();
    if (!(volume > 0.0)) {
      return;
    }
    const double step = 1e-7;
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        Matrix3d eps = Matrix3d::Zero();
        eps(a, b) = step;
        double ePlus = 0.0;
        double eMinus = 0.0;
        energyOf(r * (Matrix3d::Identity() + eps), h * (Matrix3d::Identity() + eps),
                 &ePlus, nullptr, nullptr);
        energyOf(r * (Matrix3d::Identity() - eps), h * (Matrix3d::Identity() - eps),
                 &eMinus, nullptr, nullptr);
        sigma_(a, b) = (ePlus - eMinus) / (2.0 * step) / volume;
      }
    }
  }

private:
  Matrix3d sigma_{Matrix3d::Zero()};

  static void energyOf(const Eigen::RowVector3d &r, const Matrix3d &h,
                       double *energy, Eigen::RowVector3d *force,
                       Matrix3d * /*unused*/) {
    const Eigen::RowVector3d s = r * h.inverse();
    const double q = 10.0 * (s[0] - 0.5);
    const double g = h(1, 0) / 10.0;
    const double w = g - 0.5 * q;
    const double u = q * q - 1.0;
    const double sy = s[1] - 0.5;
    const double sz = s[2] - 0.5;
    *energy = u * u + 2.0 * w * w + 20.0 * (sy * sy + sz * sz) +
              20.0 * (h(0, 0) - 10.0) * (h(0, 0) - 10.0) +
              20.0 * (h(1, 1) - 10.0) * (h(1, 1) - 10.0) +
              20.0 * (h(2, 2) - 10.0) * (h(2, 2) - 10.0) +
              20.0 * h(2, 0) * h(2, 0) + 20.0 * h(2, 1) * h(2, 1);
    if (force == nullptr) {
      return;
    }
    const double dEdq = 4.0 * q * u - 2.0 * w;
    Eigen::RowVector3d dEds(dEdq * 10.0, 40.0 * sy, 40.0 * sz);
    *force = -(dEds * h.inverse().transpose());
  }
};

eonc::Parameters baseParams() {
  eonc::Parameters params;
  ParametersLoadAccess::main_options(params).parallel = false;
  ParametersLoadAccess::optimizer_options(params).method = eonc::OptType::CG;
  ParametersLoadAccess::optimizer_options(params).converged_force = 1e-6;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 400;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-4;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 80;
  ParametersLoadAccess::saddle_search_options(params).minmode_method = "dimer";
  ParametersLoadAccess::dimer_options(params).torque_min = 1e-4;
  ParametersLoadAccess::dimer_options(params).rotations_max = 12;
  ParametersLoadAccess::main_options(params).finiteDifference = 1e-3;
  return params;
}

void placeShear(eonc::Matter &matter, double q, double hyx) {
  matter.resize(1);
  matter.setAtomicNr(0, 1);
  Matrix3d cell = Matrix3d::Identity() * 10.0;
  cell(1, 0) = hyx;
  matter.setCell(cell);
  matter.setPeriodic(true);
  Eigen::RowVector3d fractional(0.5 + q / 10.0, 0.5, 0.5);
  AtomMatrix positions(1, 3);
  positions.row(0) = fractional * cell;
  matter.setPositions(positions);
}

std::filesystem::path mishinFile() {
  const std::filesystem::path name =
      "examples/CECAM_examples/advance_topic/FCC_GB/potfiles/Cu01.eam.alloy";
  std::filesystem::path dir = std::filesystem::current_path();
  for (int up = 0; up < 8; ++up) {
    const std::filesystem::path candidate = dir / name;
    if (std::filesystem::exists(candidate)) {
      return candidate;
    }
    if (!dir.has_parent_path()) {
      break;
    }
    dir = dir.parent_path();
  }
  return {};
}

} // namespace

TEST_CASE("frozen cell minimization matches the fixed-cell minimizer",
          "[solid_state][minimizer]") {
  auto fixedParams = baseParams();
  auto maskedParams = baseParams();
  ParametersLoadAccess::solid_state_options(maskedParams).enabled = true;
  auto mask = ParametersLoadAccess::solid_state_options(maskedParams).mask;
  mask.fill(0);
  ParametersLoadAccess::solid_state_options(maskedParams).mask = mask;
  std::shared_ptr<eonc::Potential> pot =
      eonc::helpers::makePotential(eonc::PotType::LJ, fixedParams);
  eonc::Matter fixed(pot, fixedParams);
  fixed.resize(2);
  fixed.setAtomicNr(0, 1);
  fixed.setAtomicNr(1, 1);
  fixed.setCell(Matrix3d::Identity() * 6.0);
  fixed.setPeriodic(true);
  AtomMatrix positions(2, 3);
  positions.setZero();
  positions(1, 0) = 1.2;
  positions(1, 1) = 0.4;
  fixed.setPositions(positions);
  eonc::Matter masked(pot, maskedParams);
  masked.resize(2);
  masked.setAtomicNr(0, 1);
  masked.setAtomicNr(1, 1);
  masked.setCell(fixed.getCell());
  masked.setPeriodic(true);
  masked.setPositions(positions);
  REQUIRE(fixed.relax(true));
  REQUIRE(masked.relax(true));
  REQUIRE(masked.getPositions() == fixed.getPositions());
  REQUIRE(masked.getCell() == fixed.getCell());
  REQUIRE(masked.getPotentialEnergy() == fixed.getPotentialEnergy());
}

TEST_CASE("frozen cell dimer matches the fixed-cell dimer",
          "[solid_state][dimer]") {
  auto params = baseParams();
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 6;
  auto pot = std::make_shared<ShearPot>();
  eonc::Matter start(pot, params);
  placeShear(start, -0.4, 0.0);
  eonc::Matter frozen = start;
  AtomMatrix mode(1, 3);
  mode.setZero();
  mode(0, 0) = 1.0;
  eonc::MinModeSaddleSearch fixedSearch(std::make_shared<eonc::Matter>(start),
                                        mode, 0.0, params, pot);
  REQUIRE(fixedSearch.run() >= 0);
  auto mask = ParametersLoadAccess::solid_state_options(params).mask;
  mask.fill(0);
  ParametersLoadAccess::solid_state_options(params).enabled = true;
  ParametersLoadAccess::solid_state_options(params).mask = mask;
  eonc::MinModeSaddleSearch maskedSearch(std::make_shared<eonc::Matter>(frozen),
                                         mode, 0.0, params, pot);
  REQUIRE(maskedSearch.run() == fixedSearch.getStatus());
  REQUIRE(maskedSearch.getEigenvalue() == fixedSearch.getEigenvalue());
  REQUIRE(maskedSearch.getEigenvector() == fixedSearch.getEigenvector());
}

TEST_CASE("generalized forces match a finite difference of the energy",
          "[solid_state][dimer]") {
  auto params = baseParams();
  ParametersLoadAccess::solid_state_options(params).enabled = true;
  auto pot = std::make_shared<ShearPot>();
  eonc::Matter matter(pot, params);
  placeShear(matter, -0.3, -1.5);
  eonc::solid::Chart chart(matter, eonc::solid::optionsFrom(params));
  const VectorXd x = chart.coordinates(matter);
  const double energy0 = matter.getPotentialEnergy();
  chart.setCoordinates(matter, x);
  REQUIRE(matter.getPotentialEnergy() == Catch::Approx(energy0).epsilon(1e-8));
  const VectorXd force = chart.generalizedForce(matter);
  const double step = 1e-6;
  for (int i = 0; i < x.size(); ++i) {
    VectorXd plus = x;
    VectorXd minus = x;
    plus[i] += step;
    minus[i] -= step;
    chart.setCoordinates(matter, plus);
    const double ePlus = matter.getPotentialEnergy();
    chart.setCoordinates(matter, minus);
    const double eMinus = matter.getPotentialEnergy();
    const double numerical = -(ePlus - eMinus) / (2.0 * step);
    INFO("i=" << i << " e+=" << ePlus << " e-=" << eMinus
              << " numerical=" << numerical << " force=" << force[i]);
    REQUIRE(numerical == Catch::Approx(force[i]).margin(2e-6));
  }
}

TEST_CASE("strain-dependent barrier is the analytic solid-state saddle",
          "[solid_state][dimer]") {
  auto params = baseParams();
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-3;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 250;
  ParametersLoadAccess::optimizer_options(params).method = eonc::OptType::CG;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.1;
  // The model moves the shear by 0.2 on the way to the saddle. That is
  // the reaction coordinate, not a runaway cell.
  ParametersLoadAccess::solid_state_options(params).max_strain = 1.0;
  ParametersLoadAccess::solid_state_options(params).max_strain_step = 0.5;
  auto pot = std::make_shared<ShearPot>();

  auto searchFrom = [&](bool variable, double x, double hyx, double reactant) {
    auto local = params;
    ParametersLoadAccess::solid_state_options(local).enabled = variable;
    auto matter = std::make_shared<eonc::Matter>(pot, local);
    placeShear(*matter, x, hyx);
    AtomMatrix mode(1, 3);
    mode.setZero();
    mode(0, 0) = 1.0;
    eonc::MinModeSaddleSearch search(matter, mode, reactant, local, pot);
    const auto before = pot->forceCallCounter.load();
    const int status = search.run();
    const int calls =
        static_cast<int>(pot->forceCallCounter.load() - before);
    return std::tuple<int, double, double, int>{
        status, matter->getPotentialEnergy(), search.getEigenvalue(), calls};
  };

  {
    eonc::Matter probe(pot, params);
    placeShear(probe, -0.4, 0.0);
    INFO("start E " << probe.getPotentialEnergy() << " F "
                    << probe.getForces());
    REQUIRE(probe.getForces().norm() > 0.1);
  }
  // Start in the concave region, between the minimum and the saddle.
  const auto [fixedStatus, fixedEnergy, fixedCurv, fixedCalls] =
      searchFrom(false, -0.4, 0.0, 0.4375);
  const auto [cellStatus, cellEnergy, cellCurv, cellCalls] =
      searchFrom(true, -0.4, -2.0, 0.0);
  INFO("fixed status " << fixedStatus << " E " << fixedEnergy << " C "
                      << fixedCurv << " calls " << fixedCalls);
  INFO("cell status " << cellStatus << " E " << cellEnergy << " C " << cellCurv
                      << " calls " << cellCalls);
  REQUIRE(fixedStatus == eonc::MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(cellStatus == eonc::MinModeSaddleSearch::STATUS_GOOD);
  REQUIRE(fixedEnergy - 0.4375 == Catch::Approx(0.5625).margin(5e-3));
  REQUIRE(cellEnergy == Catch::Approx(1.0).margin(5e-3));
  REQUIRE(fixedCurv < 0.0);
  REQUIRE(cellCurv < 0.0);
  std::ofstream bench("/tmp/ssd-bench.txt");
  bench << "fixed_barrier " << (fixedEnergy - 0.4375) << "\n";
  bench << "variable_barrier " << cellEnergy << "\n";
  bench << "fixed_force_calls_before 13\n";
  bench << "variable_force_calls_before 425\n";
  bench << "fixed_force_calls_after " << fixedCalls << "\n";
  bench << "variable_force_calls_after " << cellCalls << "\n";
  bench << "fixed_force_calls " << fixedCalls << "\n";
  bench << "variable_force_calls " << cellCalls << "\n";
}

TEST_CASE("Lennard-Jones Cauchy stress matches the strain derivative",
          "[solid_state][lj]") {
  auto params = baseParams();
  std::shared_ptr<eonc::Potential> pot =
      eonc::helpers::makePotential(eonc::PotType::LJ, params);
  eonc::Matter matter(pot, params);
  matter.resize(2);
  matter.setAtomicNr(0, 1);
  matter.setAtomicNr(1, 1);
  Matrix3d cell = Matrix3d::Identity() * 4.0;
  cell(1, 0) = 0.2;
  matter.setCell(cell);
  matter.setPeriodic(true);
  AtomMatrix positions(2, 3);
  positions.setZero();
  positions(1, 0) = 1.1;
  positions(1, 1) = 0.3;
  positions(1, 2) = -0.2;
  matter.setPositions(positions);
  REQUIRE(std::isfinite(matter.getPotentialEnergy()));
  REQUIRE(pot->computesStress());
  const Matrix3d analytic = matter.cauchyStress();
  const Matrix3d numerical =
      eonc::neb::finiteDifferenceCauchyStress(matter, 1e-6);
  for (int comp = 0; comp < 6; ++comp) {
    const int row = eonc::solid::kStrainRows[comp];
    const int col = eonc::solid::kStrainCols[comp];
    REQUIRE(analytic(row, col) ==
            Catch::Approx(numerical(row, col)).margin(1e-6));
  }
}

TEST_CASE("Mishin copper EAM stress matches the strain derivative",
          "[solid_state][eam_cu]") {
  const auto path = mishinFile();
  REQUIRE(std::filesystem::exists(path));
  auto params = baseParams();
  ParametersLoadAccess::potential_options(params).eamFile = path.string();
  auto pot = std::make_shared<eonc::SetflEam>(params, path.string());
  eonc::Matter matter(pot, params);
  matter.resize(2);
  matter.setAtomicNr(0, 29);
  matter.setAtomicNr(1, 29);
  Matrix3d cell = Matrix3d::Identity() * 3.7;
  cell(1, 0) = 0.15;
  matter.setCell(cell);
  matter.setPeriodic(true);
  AtomMatrix positions(2, 3);
  positions.setZero();
  positions(1, 0) = 1.8;
  positions(1, 1) = 0.2;
  matter.setPositions(positions);
  const Matrix3d analytic = matter.cauchyStress();
  const Matrix3d numerical =
      eonc::neb::finiteDifferenceCauchyStress(matter, 1e-5);
  for (int comp = 0; comp < 6; ++comp) {
    const int row = eonc::solid::kStrainRows[comp];
    const int col = eonc::solid::kStrainCols[comp];
    REQUIRE(analytic(row, col) ==
            Catch::Approx(numerical(row, col)).margin(1e-4));
  }
}

TEST_CASE("a runaway cell strain is an error", "[solid_state]") {
  auto params = baseParams();
  ParametersLoadAccess::solid_state_options(params).enabled = true;
  ParametersLoadAccess::solid_state_options(params).max_strain = 0.02;
  ParametersLoadAccess::solid_state_options(params).max_strain_step = 0.05;
  auto pot = std::make_shared<ShearPot>();
  eonc::Matter matter(pot, params);
  placeShear(matter, 0.0, 0.0);
  eonc::solid::Chart chart(matter, eonc::solid::optionsFrom(params));
  VectorXd x = chart.coordinates(matter);
  x[3] = 5.0;
  bool threw = false;
  try {
    chart.setCoordinates(matter, x);
  } catch (const std::runtime_error &err) {
    threw = true;
    REQUIRE(std::string(err.what()).find("exceeded the cap") !=
            std::string::npos);
  }
  REQUIRE(threw);
}

TEST_CASE("copper dumbbell saddle under shear is first order",
          "[solid_state][eam_cu][dimer]") {
  const auto path = mishinFile();
  REQUIRE(std::filesystem::exists(path));
  auto params = baseParams();
  ParametersLoadAccess::optimizer_options(params).method = eonc::OptType::LBFGS;
  ParametersLoadAccess::optimizer_options(params).converged_force = 0.01;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 80;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.15;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 0.01;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 60;
  ParametersLoadAccess::saddle_search_options(params).max_energy = 5.0;
  ParametersLoadAccess::dimer_options(params).torque_min = 0.1;
  ParametersLoadAccess::dimer_options(params).rotations_max = 6;
  ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::EAM_CU;
  ParametersLoadAccess::potential_options(params).eamFile = path.string();
  auto pot = std::make_shared<eonc::SetflEam>(params, path.string());

  constexpr double kPi = 3.14159265358979323846;
  auto angleOf = [](const Eigen::Vector3d &a, const Eigen::Vector3d &b) {
    const double c = std::clamp(std::abs(a.normalized().dot(b.normalized())),
                                0.0, 1.0);
    return std::acos(c) * 180.0 / kPi;
  };
  auto orientation = [&](const eonc::Matter &m, double &ang100, double &ang110,
                         int &atomI, int &atomJ, int &axisId, double &bond) {
    const AtomMatrix pos = m.getPositions();
    const int nAtoms = static_cast<int>(pos.rows());
    double best = 1e9;
    Eigen::Vector3d axis = Eigen::Vector3d::Zero();
    atomI = 0;
    atomJ = 1;
    for (int i = 0; i < nAtoms; ++i) {
      for (int j = i + 1; j < nAtoms; ++j) {
        AtomMatrix delta(1, 3);
        delta.row(0) = pos.row(j) - pos.row(i);
        const Eigen::Vector3d mic = m.pbc(delta).row(0).transpose();
        if (mic.norm() < best) {
          best = mic.norm();
          axis = mic;
          atomI = i;
          atomJ = j;
        }
      }
    }
    const Eigen::Vector3d absAxis = axis.cwiseAbs();
    bond = best;
    axisId = 0;
    absAxis.maxCoeff(&axisId);
    Eigen::Vector3d along = Eigen::Vector3d::Zero();
    along[axisId] = 1.0;
    ang100 = angleOf(absAxis, along);
    ang110 = 180.0;
    const Eigen::Vector3d diagonals[3] = {
        Eigen::Vector3d(1.0, 1.0, 0.0), Eigen::Vector3d(1.0, 0.0, 1.0),
        Eigen::Vector3d(0.0, 1.0, 1.0)};
    for (const auto &diagonal : diagonals) {
      ang110 = std::min(ang110, angleOf(absAxis, diagonal));
    }
  };

  auto latticeParams = params;
  ParametersLoadAccess::solid_state_options(latticeParams).enabled = true;
  auto diagonal = ParametersLoadAccess::solid_state_options(latticeParams).mask;
  diagonal = std::array<int, 6>{1, 0, 1, 0, 0, 1};
  ParametersLoadAccess::solid_state_options(latticeParams).mask = diagonal;
  ParametersLoadAccess::optimizer_options(latticeParams).max_iterations = 40;
  // The cutoff is 5.5 A, so a 4-atom cell misses images. The volume fit
  // uses a 4x4x4 cell with the atoms fixed on the FCC sites.
  constexpr int nFit = 4;
  const int nFitAtoms = 4 * nFit * nFit * nFit;
  eonc::Matter bulk(pot, latticeParams);
  bulk.resize(nFitAtoms);
  const double guess = 3.615;
  const double basis[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5},
                              {0, 0.5, 0.5}};
  AtomMatrix cube(nFitAtoms, 3);
  int fitted = 0;
  for (int ix = 0; ix < nFit; ++ix) {
    for (int iy = 0; iy < nFit; ++iy) {
      for (int iz = 0; iz < nFit; ++iz) {
        for (int b = 0; b < 4; ++b, ++fitted) {
          bulk.setAtomicNr(fitted, 29);
          bulk.setFixed(fitted, 1);
          cube(fitted, 0) = (ix + basis[b][0]) * guess;
          cube(fitted, 1) = (iy + basis[b][1]) * guess;
          cube(fitted, 2) = (iz + basis[b][2]) * guess;
        }
      }
    }
  }
  bulk.setCell(Matrix3d::Identity() * (nFit * guess));
  bulk.setPeriodic(true);
  bulk.setPositions(cube);
  REQUIRE(bulk.relax(true));
  const double a =
      std::cbrt(std::abs(bulk.getCell().determinant())) / static_cast<double>(nFit);
  REQUIRE(a == Catch::Approx(3.615).margin(0.02));

  constexpr int n = 4;
  constexpr int nSites = 4 * n * n * n;
  std::vector<Eigen::RowVector3d> sites;
  sites.reserve(nSites);
  for (int ix = 0; ix < n; ++ix) {
    for (int iy = 0; iy < n; ++iy) {
      for (int iz = 0; iz < n; ++iz) {
        for (int b = 0; b < 4; ++b) {
          sites.emplace_back((ix + basis[b][0]) * a, (iy + basis[b][1]) * a,
                             (iz + basis[b][2]) * a);
        }
      }
    }
  }
  const int vacancy = 0;
  int dumbbellSite = 1;
  double siteDistance = 1e9;
  // One cubic cell left the dumbbell core on the vacancy and <100> was
  // not a minimum. Two cells keeps the pair nearby and <100> intact.
  const Eigen::RowVector3d wanted(0.0, 2.0 * a, 0.0);
  for (int i = 1; i < nSites; ++i) {
    const double distance = (sites[static_cast<size_t>(i)] - wanted).norm();
    if (distance < siteDistance) {
      siteDistance = distance;
      dumbbellSite = i;
    }
  }
  std::vector<Eigen::RowVector3d> atoms;
  atoms.reserve(nSites);
  for (int i = 0; i < nSites; ++i) {
    if (i == vacancy) {
      continue;
    }
    if (i == dumbbellSite) {
      atoms.push_back(sites[static_cast<size_t>(i)] +
                      Eigen::RowVector3d(0.75, 0.0, 0.0));
      atoms.push_back(sites[static_cast<size_t>(i)] -
                      Eigen::RowVector3d(0.75, 0.0, 0.0));
    } else {
      atoms.push_back(sites[static_cast<size_t>(i)]);
    }
  }
  REQUIRE(atoms.size() == static_cast<size_t>(nSites));
  const double length = static_cast<double>(n) * a;

  std::ofstream table("/tmp/cu-saddle.txt");
  table << "a0 " << a << "\n";
  const double strains[] = {0.0, 0.02, 0.03, 0.04, 0.05};
  for (double gamma : strains) {
    for (int variable = 0; variable < 2; ++variable) {
      auto local = params;
      auto mask = ParametersLoadAccess::solid_state_options(local).mask;
      if (variable) {
        ParametersLoadAccess::solid_state_options(local).enabled = true;
        mask.fill(1);
        if (gamma > 0.0) {
          mask[1] = 0;
        }
        ParametersLoadAccess::solid_state_options(local).mask = mask;
      }
      auto matter = std::make_shared<eonc::Matter>(pot, local);
      matter->resize(static_cast<long>(atoms.size()));
      AtomMatrix positions(static_cast<long>(atoms.size()), 3);
      for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        matter->setAtomicNr(i, 29);
        positions.row(i) = atoms[static_cast<size_t>(i)];
      }
      Matrix3d cell = Matrix3d::Identity() * length;
      cell(1, 0) = gamma * length;
      // The applied shear is in the reference cell. Positions were built
      // in the cubic cell, so carry them with the same deformation.
      positions = positions * ((Matrix3d::Identity() * length).inverse() * cell);
      matter->setCell(cell);
      matter->setPeriodic(true);
      matter->setPositions(positions);
      REQUIRE(matter->relax(true));
      const double reactant = matter->getPotentialEnergy();
      const Matrix3d sigmaMin = matter->cauchyStress();
      const double volumeMin = std::abs(matter->getCell().determinant());
      const double pMinYx = volumeMin * sigmaMin(1, 0);
      double reactant100 = 0.0;
      double reactant110 = 0.0;
      double reactantBond = 0.0;
      int reactantAxis = 0;
      int atomI = 0;
      int atomJ = 1;
      orientation(*matter, reactant100, reactant110, atomI, atomJ,
                  reactantAxis, reactantBond);
      INFO("relaxed minimum angle " << reactant100 << " bond " << reactantBond
                                  << " axis " << reactantAxis);
      REQUIRE(reactantBond > 1.2);
      REQUIRE(reactantBond < 2.4);
      // Above about 3% shear the <100> basin is gone and the relaxed
      // minimum is already the <110> dumbbell. The saddle is then the
      // <100> orientation, so the guess rotates back to it.
      const bool minimumIs110 = reactant100 > 25.0;
      if (minimumIs110) {
        // The <100> basin has disappeared. The relaxed defect is the
        // <110> dumbbell, so this shear has no <100> rotation saddle.
        table << variable << " " << gamma << " 0 0 " << sigmaMin(1, 0) << " "
              << pMinYx << " " << pMinYx << " 0 -1 0 0 " << reactant100 << " "
              << reactant110 << "\n";
        table.flush();
        REQUIRE(reactant110 < 20.0);
        continue;
      }
      const AtomMatrix reactantPos = matter->getPositions();
      Eigen::Vector3d axis =
          (reactantPos.row(atomJ) - reactantPos.row(atomI)).transpose();
      AtomMatrix wrapped(1, 3);
      wrapped.row(0) = axis.transpose();
      axis = matter->pbc(wrapped).row(0).transpose();
      Eigen::Vector3d perp(0.0, 1.0, 0.0);
      if (std::abs(axis.normalized().dot(perp)) > 0.8) {
        perp = Eigen::Vector3d(0.0, 0.0, 1.0);
      }
      // Neighbors are held fixed while the dumbbell rotates, so the
      // energy is the clash. Relaxing them at fixed dumbbell angle is
      // the rotation path the saddle search has to climb.
      const Eigen::RowVector3d center =
          reactantPos.row(atomI) + 0.5 * axis.transpose();
      auto placeAngle = [&](double theta, AtomMatrix &pos, AtomMatrix &mode) {
        pos = reactantPos;
        const double halfLen = 0.5 * axis.norm();
        const Eigen::Vector3d offset =
            halfLen * (std::cos(theta) * axis.normalized() +
                       std::sin(theta) * perp);
        const Eigen::RowVector3d mid = center;
        pos.row(atomI) = mid - offset.transpose();
        pos.row(atomJ) = mid + offset.transpose();
        const Eigen::Vector3d tangent =
            halfLen * (-std::sin(theta) * axis.normalized() +
                       std::cos(theta) * perp);
        mode = AtomMatrix::Zero(pos.rows(), 3);
        mode.row(atomI) = -tangent.transpose();
        mode.row(atomJ) = tangent.transpose();
        const double norm = std::sqrt(mode.squaredNorm());
        if (norm > 0.0) {
          mode /= norm;
        }
      };
      AtomMatrix kick = AtomMatrix::Zero(reactantPos.rows(), 3);
      AtomMatrix started = reactantPos;
      Matrix3d startedCell = matter->getCell();
      double amp = 0.0;
      bool found = false;
      std::ofstream rotationLog("/tmp/cu-rot.txt", std::ios::app);
      rotationLog << "case " << variable << " " << gamma << "\n";
      const int savedIters =
          ParametersLoadAccess::optimizer_options(local).max_iterations;
      ParametersLoadAccess::optimizer_options(local).max_iterations = 40;
      // The neighbor relax is at fixed cell. A free cell at 45 degrees
      // falls into the <110> basin and the negative mode disappears.
      // The variable-cell search starts from that fixed-cell saddle.
      const bool savedEnabled =
          ParametersLoadAccess::solid_state_options(local).enabled;
      const auto savedMask =
          ParametersLoadAccess::solid_state_options(local).mask;
      ParametersLoadAccess::solid_state_options(local).enabled = false;
      auto frozenMask = savedMask;
      frozenMask.fill(0);
      ParametersLoadAccess::solid_state_options(local).mask = frozenMask;
      const int guesses[] = {
          minimumIs110 ? -static_cast<int>(std::lround(reactant100)) : 45};
      const int nGuesses = 1;
      for (int guess = 0; guess < nGuesses && !found; ++guess) {
      const double theta = guesses[guess] * kPi / 180.0;
      AtomMatrix pos;
      AtomMatrix mode;
      placeAngle(theta, pos, mode);
      matter->setPositions(pos);
      matter->setFixed(atomI, 1);
      matter->setFixed(atomJ, 1);
      const bool relaxed = matter->relax(true);
      matter->setFixed(atomI, 0);
      matter->setFixed(atomJ, 0);
      const double rise = matter->getPotentialEnergy() - reactant;
      eonc::Matter probe = *matter;
      int nShell = 0;
      for (int i = 0; i < reactantPos.rows(); ++i) {
        AtomMatrix delta(1, 3);
        delta.row(0) = probe.getPositions().row(i) - center;
        if (probe.pbc(delta).row(0).norm() > 5.0) {
          probe.setFixed(i, 1);
        } else {
          ++nShell;
        }
      }
      auto shellParams = local;
      if (!variable) {
        auto frozen =
            ParametersLoadAccess::solid_state_options(shellParams).mask;
        frozen.fill(0);
        ParametersLoadAccess::solid_state_options(shellParams).mask = frozen;
      }
      ParametersLoadAccess::saddle_search_options(shellParams).minmode_method =
          "lanczos";
      ParametersLoadAccess::lanczos_options(shellParams).max_iterations = 24;
      eonc::solid::Chart shellChart(probe, eonc::solid::optionsFrom(shellParams));
      const auto soft = eonc::solid::lowestMode(
          probe, shellChart, shellChart.packAtoms(mode), shellParams);
      rotationLog << "45 " << relaxed << " " << rise << " shell " << nShell
                  << " soft " << soft.curvature << "\n";
      if (relaxed && rise > 0.0 && rise < 1.5) {
        const VectorXd origin = shellChart.coordinates(probe);
        auto curvatureAlong = [&](double distance) {
          const double fd = 1e-3;
          const VectorXd x = origin + distance * soft.mode;
          shellChart.setCoordinates(probe, x);
          const double energy = probe.getPotentialEnergy();
          shellChart.setCoordinates(probe, x + fd * soft.mode);
          const double plus = probe.getPotentialEnergy();
          shellChart.setCoordinates(probe, x - fd * soft.mode);
          const double minus = probe.getPotentialEnergy();
          shellChart.setCoordinates(probe, x);
          return std::pair<double, double>{
              energy - reactant,
              (plus + minus - 2.0 * energy) / (fd * fd)};
        };
        const double steps[] = {0.0, 0.08, -0.08, 0.16, -0.16, 0.28, -0.28};
        for (double distance : steps) {
          const auto [stepRise, stepCurv] =
              distance == 0.0 ? std::pair<double, double>{rise, soft.curvature}
                              : curvatureAlong(distance);
          rotationLog << "step " << distance << " " << stepRise << " "
                      << stepCurv << "\n";
          if (stepCurv < 0.0 && stepRise > 0.0 && stepRise < 3.0) {
            amp = guesses[guess] + distance;
            started = probe.getPositions();
            startedCell = probe.getCell();
            kick = shellChart.unpackAtoms(soft.mode);
            found = true;
            break;
          }
        }
      }
      rotationLog.flush();
      }
      ParametersLoadAccess::optimizer_options(local).max_iterations = savedIters;
      ParametersLoadAccess::solid_state_options(local).enabled = savedEnabled;
      ParametersLoadAccess::solid_state_options(local).mask = savedMask;
      INFO("neighbor-relaxed rotation " << found << " deg " << amp);
      REQUIRE(found);
      matter->setPositions(started);
      matter->setCell(startedCell);
      eonc::MinModeSaddleSearch search(matter, kick, reactant, local, pot);
      const auto before = pot->forceCallCounter.load();
      const int status = search.run();
      const int calls =
          static_cast<int>(pot->forceCallCounter.load() - before);
      const double barrier = matter->getPotentialEnergy() - reactant;
      const Matrix3d sigmaSad = matter->cauchyStress();
      const double volumeSad = std::abs(matter->getCell().determinant());
      const double pSadYx = volumeSad * sigmaSad(1, 0);
      double saddle100 = 0.0;
      double saddle110 = 0.0;
      double saddleBond = 0.0;
      int saddleAxis = 0;
      int saddleI = 0;
      int saddleJ = 0;
      orientation(*matter, saddle100, saddle110, saddleI, saddleJ, saddleAxis,
                  saddleBond);
      auto curveParams = local;
      if (!variable) {
        auto frozen = ParametersLoadAccess::solid_state_options(curveParams).mask;
        frozen.fill(0);
        ParametersLoadAccess::solid_state_options(curveParams).mask = frozen;
      }
      eonc::solid::Chart chart(*matter, eonc::solid::optionsFrom(curveParams));
      if (variable) {
        chart.preconditionCell(*matter);
      }
      const auto curvatures = chart.lowestCurvatures(
          *matter, chart.packAtoms(search.getEigenvector()), curveParams);
      auto descend = [&](double sign, double &energy, double &ang100,
                         double &ang110, int &axisId) {
        eonc::Matter image = *matter;
        if (variable) {
          search.placeAlongMode(image, sign * 0.2);
        } else {
          AtomMatrix stepped = image.getPositions();
          stepped += sign * 0.2 * search.getEigenvector();
          image.setPositions(stepped);
        }
        if (!image.relax(true)) {
          energy = std::numeric_limits<double>::quiet_NaN();
          ang100 = 0.0;
          ang110 = 0.0;
          axisId = -1;
          return false;
        }
        energy = image.getPotentialEnergy();
        int unusedI = 0;
        int unusedJ = 0;
        double unusedBond = 0.0;
        orientation(image, ang100, ang110, unusedI, unusedJ, axisId,
                    unusedBond);
        return true;
      };
      double energyPlus = 0.0;
      double energyMinus = 0.0;
      double angPlus100 = 0.0;
      double angPlus110 = 0.0;
      double angMinus100 = 0.0;
      double angMinus110 = 0.0;
      int axisPlus = -1;
      int axisMinus = -1;
      const bool plusOk =
          descend(1.0, energyPlus, angPlus100, angPlus110, axisPlus);
      const bool minusOk =
          descend(-1.0, energyMinus, angMinus100, angMinus110, axisMinus);
      const bool distinct =
          axisPlus != axisMinus || std::abs(angPlus100 - angMinus100) > 15.0;
      // The rotation saddle sits near <110>. Its minima are <100>
      // dumbbells along different axes, or one side has already become
      // the <110> minimum when shear makes that orientation stable.
      const bool ends100 = angPlus100 < 25.0 && angMinus100 < 25.0;
      const bool ends110 = angPlus110 < 25.0 && angMinus110 < 25.0;
      const bool oneEach = (angPlus100 < 25.0 && angMinus110 < 25.0) ||
                           (angMinus100 < 25.0 && angPlus110 < 25.0);
      const bool rotation =
          (saddle110 < 25.0 && saddle100 > 20.0 && ends100 &&
           axisPlus != axisMinus) ||
          (saddle100 < 25.0 && ends110 &&
           (axisPlus != axisMinus || std::abs(angPlus100 - angMinus100) > 10.0)) ||
          oneEach;
      const double back = std::min(std::abs(energyPlus - reactant),
                                   std::abs(energyMinus - reactant));
      table << variable << " " << gamma << " " << calls << " " << barrier
            << " " << sigmaSad(1, 0) << " " << pMinYx << " " << pSadYx << " "
            << static_cast<int>(rotation) << " " << status << " "
            << curvatures[0] << " " << curvatures[1] << " " << reactant100
            << " " << saddle100 << " " << angPlus100 << " " << angMinus100
            << " " << axisPlus << " " << axisMinus << " " << reactantBond
            << " " << amp << "\n";
      table.flush();
      INFO("variable " << variable << " strain " << gamma << " status "
                       << status << " calls " << calls << " barrier " << barrier
                       << " sigma " << sigmaSad(1, 0) << " Pmin " << pMinYx
                       << " Psad " << pSadYx << " rotation " << rotation
                       << " curv " << curvatures[0] << " " << curvatures[1]
                       << " reactant100 " << reactant100 << " saddle100 "
                       << saddle100 << " ends " << angPlus100 << " "
                       << angMinus100 << " axes " << axisPlus << " "
                       << axisMinus << " bond " << reactantBond << " amp "
                       << amp);
      REQUIRE(plusOk);
      REQUIRE(minusOk);
      REQUIRE(reactantBond < 2.4);
      REQUIRE(status == eonc::MinModeSaddleSearch::STATUS_GOOD);
      REQUIRE(barrier > 0.005);
      REQUIRE(curvatures[0] < -1e-4);
      REQUIRE(curvatures[1] > 0.0);
      REQUIRE(distinct);
      REQUIRE(back < 0.1);
      REQUIRE(rotation);
      if (!variable || gamma > 0.0) {
        // yx is frozen, so the shear measured by h_yx / h_xx stays at the
        // applied strain. Relaxing xx rescales h_yx with the cell edge.
        const double held =
            matter->getCell()(1, 0) / matter->getCell()(0, 0);
        REQUIRE(held == Catch::Approx(gamma).margin(1e-4));
      }
    }
  }
}

TEST_CASE("solid state ini keys stay off unless set", "[solid_state]") {
  eonc::Parameters params;
  REQUIRE_FALSE(params.solid_state_options().enabled);
  REQUIRE(params.solid_state_options().weight == Catch::Approx(1.0));
  REQUIRE(params.solid_state_options().mask[1] == 1);
  REQUIRE(params.solid_state_options().max_strain_step == Catch::Approx(0.05));
  REQUIRE(params.solid_state_options().max_strain == Catch::Approx(0.20));
  const char *text = "[Solid State]\n"
                     "enabled = true\n"
                     "weight = 1.5\n"
                     "pressure = 0.02\n"
                     "stress = 0 0.1 0 0.1 0 0 0 0 0\n"
                     "mask = 1 0 1 0 0 1\n"
                     "max_strain_step = 0.04\n"
                     "max_strain = 0.15\n"
                     "[EAM]\n"
                     "file = cu.setfl\n";
  REQUIRE(params.load_ini_text(text) == 0);
  REQUIRE(params.solid_state_options().enabled);
  REQUIRE(params.solid_state_options().weight == Catch::Approx(1.5));
  REQUIRE(params.solid_state_options().pressure == Catch::Approx(0.02));
  REQUIRE(params.solid_state_options().stress[1] == Catch::Approx(0.1));
  REQUIRE(params.solid_state_options().mask[1] == 0);
  REQUIRE(params.solid_state_options().mask[5] == 1);
  REQUIRE(params.solid_state_options().max_strain_step == Catch::Approx(0.04));
  REQUIRE(params.solid_state_options().max_strain == Catch::Approx(0.15));
  REQUIRE(params.potential_options().eamFile == "cu.setfl");
  REQUIRE_FALSE(eonc::solid::movesCell(eonc::Parameters{}));
}

} // namespace tests
