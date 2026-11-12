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

#include <filesystem>
#include <fstream>
#include <sstream>

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

TEST_CASE("copper dumbbell saddle under shear moves the cell",
          "[solid_state][eam_cu][dimer]") {
  const auto path = mishinFile();
  REQUIRE(std::filesystem::exists(path));
  auto params = baseParams();
  ParametersLoadAccess::optimizer_options(params).method = eonc::OptType::CG;
  ParametersLoadAccess::optimizer_options(params).converged_force = 0.02;
  ParametersLoadAccess::optimizer_options(params).max_iterations = 60;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 0.05;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 40;
  ParametersLoadAccess::dimer_options(params).torque_min = 0.05;
  ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::EAM_CU;
  ParametersLoadAccess::potential_options(params).eamFile = path.string();
  auto pot = std::make_shared<eonc::SetflEam>(params, path.string());

  const double a = 3.615;
  const int n = 2;
  const int nLattice = 4 * n * n * n;
  AtomMatrix lattice(nLattice, 3);
  const double basis[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5},
                              {0, 0.5, 0.5}};
  int placed = 0;
  for (int ix = 0; ix < n; ++ix) {
    for (int iy = 0; iy < n; ++iy) {
      for (int iz = 0; iz < n; ++iz) {
        for (int b = 0; b < 4; ++b, ++placed) {
          lattice(placed, 0) = (ix + basis[b][0]) * a;
          lattice(placed, 1) = (iy + basis[b][1]) * a;
          lattice(placed, 2) = (iz + basis[b][2]) * a;
        }
      }
    }
  }
  // Vacancy at the last site, <100> dumbbell on the first site.
  AtomMatrix dumbbell(nLattice, 3);
  dumbbell.topRows(nLattice - 1) = lattice.topRows(nLattice - 1);
  const Eigen::RowVector3d site = lattice.row(0);
  dumbbell.row(0) = site + Eigen::RowVector3d(0.72, 0.0, 0.0);
  dumbbell.row(nLattice - 1) = site - Eigen::RowVector3d(0.72, 0.0, 0.0);

  auto run = [&](bool variable) {
    auto local = params;
    ParametersLoadAccess::solid_state_options(local).enabled = variable;
    if (variable) {
      // Shear the cell. The diagonal stays free as well.
      auto mask = ParametersLoadAccess::solid_state_options(local).mask;
      mask.fill(1);
      ParametersLoadAccess::solid_state_options(local).mask = mask;
      ParametersLoadAccess::solid_state_options(local).stress[1] = 0.001;
      ParametersLoadAccess::solid_state_options(local).stress[3] = 0.001;
    }
    auto matter = std::make_shared<eonc::Matter>(pot, local);
    matter->resize(nLattice);
    for (int i = 0; i < nLattice; ++i) {
      matter->setAtomicNr(i, 29);
    }
    matter->setCell(Matrix3d::Identity() * (n * a));
    matter->setPeriodic(true);
    matter->setPositions(dumbbell);
    matter->relax(true);
    const double reactant = matter->getPotentialEnergy();
    AtomMatrix kicked = matter->getPositions();
    kicked.row(0) += Eigen::RowVector3d(-0.15, 0.25, 0.0);
    kicked.row(nLattice - 1) += Eigen::RowVector3d(0.15, -0.25, 0.0);
    matter->setPositions(kicked);
    AtomMatrix mode = AtomMatrix::Zero(nLattice, 3);
    mode.row(0) = Eigen::RowVector3d(-1.0, 1.0, 0.0);
    mode.row(nLattice - 1) = Eigen::RowVector3d(1.0, -1.0, 0.0);
    eonc::MinModeSaddleSearch search(matter, mode, reactant, local, pot);
    const int status = search.run();
    Eigen::Vector3d axis = (matter->getPositions().row(0) -
                            matter->getPositions().row(nLattice - 1))
                               .transpose();
    AtomMatrix wrapped(1, 3);
    wrapped.row(0) = axis.transpose();
    axis = matter->pbc(wrapped).row(0).transpose();
    const double angle = std::atan2(std::abs(axis[1]), std::abs(axis[0])) *
                         180.0 / 3.141592653589793;
    return std::tuple<int, double, double, double>{
        status, matter->getPotentialEnergy() - reactant,
        std::abs(matter->getCell()(1, 0)), angle};
  };

  const auto [fixedStatus, fixedBarrier, fixedShear, fixedAngle] = run(false);
  const auto [cellStatus, cellBarrier, cellShear, cellAngle] = run(true);
  INFO("fixed status " << fixedStatus << " barrier " << fixedBarrier
                      << " shear " << fixedShear << " angle " << fixedAngle);
  INFO("cell status " << cellStatus << " barrier " << cellBarrier << " shear "
                      << cellShear << " angle " << cellAngle);
  // A fixed cell cannot follow the shear that turns <100> into <110>.
  // The same kick under an applied shear moves the cell and lowers the
  // energy the search reaches.
  REQUIRE(fixedShear == Catch::Approx(0.0).margin(1e-8));
  REQUIRE(cellShear > 0.2);
  REQUIRE(cellBarrier < fixedBarrier - 0.5);
  REQUIRE(cellAngle > 20.0);
  std::ofstream bench("/tmp/ssd-bench.txt", std::ios::app);
  bench << "cu_fixed_barrier " << fixedBarrier << "\n";
  bench << "cu_variable_barrier " << cellBarrier << "\n";
  bench << "cu_fixed_shear " << fixedShear << "\n";
  bench << "cu_variable_shear " << cellShear << "\n";
  bench << "cu_fixed_angle " << fixedAngle << "\n";
  bench << "cu_variable_angle " << cellAngle << "\n";
}

TEST_CASE("solid state ini keys stay off unless set", "[solid_state]") {
  eonc::Parameters params;
  REQUIRE_FALSE(params.solid_state_options().enabled);
  REQUIRE(params.solid_state_options().weight == Catch::Approx(1.0));
  REQUIRE(params.solid_state_options().mask[1] == 1);
  const char *text = "[Solid State]\n"
                     "enabled = true\n"
                     "weight = 1.5\n"
                     "pressure = 0.02\n"
                     "stress = 0 0.1 0 0.1 0 0 0 0 0\n"
                     "mask = 1 0 1 0 0 1\n"
                     "[EAM]\n"
                     "file = cu.setfl\n";
  REQUIRE(params.load_ini_text(text) == 0);
  REQUIRE(params.solid_state_options().enabled);
  REQUIRE(params.solid_state_options().weight == Catch::Approx(1.5));
  REQUIRE(params.solid_state_options().pressure == Catch::Approx(0.02));
  REQUIRE(params.solid_state_options().stress[1] == Catch::Approx(0.1));
  REQUIRE(params.solid_state_options().mask[1] == 0);
  REQUIRE(params.solid_state_options().mask[5] == 1);
  REQUIRE(params.potential_options().eamFile == "cu.setfl");
  REQUIRE_FALSE(eonc::solid::movesCell(eonc::Parameters{}));
}

} // namespace tests
