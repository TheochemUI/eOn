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
#include <cmath>
#include <limits>
#include <stdexcept>

#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Matter.h"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// Guards the RgpotAdapter + parameter mapping for every migrated arm: the
// same pinned numbers as the pre-migration in-tree kernels, evaluated
// through helpers::makePotential exactly as jobs do. Kernel numerics
// themselves are pinned upstream in rgpot's own test suite.

TEST_CASE("Adapter LJ matches the pinned reference", "[pot][rgpot-adapter]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE(matter->getPeriodic());
  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-39.965351).epsilon(1e-4));
  REQUIRE(matter->getForces().rowwise().norm().maxCoeff() ==
          Catch::Approx(0.004704).epsilon(1e-2));
}

TEST_CASE("Adapter clones preserve an independent rgpot kernel",
          "[pot][rgpot-adapter]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  const std::shared_ptr<Potential> clone = pot->clonePotential();
  REQUIRE(clone);
  REQUIRE(clone.get() != pot.get());

  Matter original(pot, params);
  Matter copied(clone, params);
  REQUIRE(eonc::io::io_ok(original.con2matter(std::string("reactant.con"))));
  REQUIRE(eonc::io::io_ok(copied.con2matter(std::string("reactant.con"))));
  REQUIRE(copied.getPotentialEnergy() ==
          Catch::Approx(-39.965351).epsilon(1e-4));
  REQUIRE(copied.getForces().rowwise().norm().maxCoeff() ==
          Catch::Approx(0.004704).epsilon(1e-2));
}

TEST_CASE("Adapter LJCluster matches the pinned reference",
          "[pot][rgpot-adapter]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential =
      PotType::LJCLUSTER;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE_FALSE(matter->getPeriodic());
  REQUIRE(matter->getPotentialEnergy() ==
          Catch::Approx(-39.965379).epsilon(1e-4));
  AtomMatrix pos = matter->getPositionsCopy();
  pos(0, 2) = -1.0e-6;
  matter->setPositions(pos);
  REQUIRE(matter->getPositions()(0, 2) == Catch::Approx(-1.0e-6).margin(1e-12));
}

TEST_CASE("Adapter LJ reports coincident atoms as a domain error",
          "[pot][rgpot-adapter]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  Matter matter(pot, params);
  REQUIRE(eonc::io::io_ok(matter.con2matter(std::string("reactant.con"))));
  AtomMatrix positions = matter.getPositionsCopy();
  positions.row(1) = positions.row(0);
  matter.setPositions(positions);
  REQUIRE_THROWS_AS(matter.getPotentialEnergy(), std::domain_error);
}

// rgpot's Morse rewrites a neighbour scratch on every call, so it asks for
// an instance per thread.
TEST_CASE("Adapter Morse evaluates and reports per-instance",
          "[pot][rgpot-adapter]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::MORSE_PT;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));

  REQUIRE_FALSE(pot->isThreadSafe());
  REQUIRE_FALSE(pot->isSharedInstanceThreadSafe());
  REQUIRE(pot->needsPerImageInstance());

  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
  REQUIRE(matter->getForces().allFinite());
}

TEST_CASE("Adapter ZBL maps config.ini cutoffs through",
          "[pot][rgpot-adapter]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::ZBL;
  ParametersLoadAccess::zbl_options(params).cut_inner = 2.0;
  ParametersLoadAccess::zbl_options(params).cut_global = 2.5;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  REQUIRE(std::isfinite(matter->getPotentialEnergy()));
}

// A batch of one repeated system passes even when an implementation
// evaluates the first system and broadcasts the answer, so the two systems
// here differ. LJ takes the default per-system loop rather than a native
// batch, which is the arm every migrated kernel inherits.
TEST_CASE("Adapter forceBatch agrees with per-system force calls",
          "[pot][rgpot-adapter]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));

  // Checked, because the displacement below indexes the first atom: an
  // unread file would leave an empty matrix and turn this into a segfault
  // instead of a failure naming the missing input.
  Matter m0(pot, params);
  REQUIRE(eonc::io::io_ok(m0.con2matter(std::string("reactant.con"))));
  Matter m1(pot, params);
  REQUIRE(eonc::io::io_ok(m1.con2matter(std::string("reactant.con"))));
  AtomMatrix displaced = m1.getPositionsCopy();
  displaced(0, 0) += 0.15;
  m1.setPositions(displaced);

  const long nAtoms = m0.numberOfAtoms();
  const auto nrs0 = m0.getAtomicNrs();
  const auto nrs1 = m1.getAtomicNrs();
  const auto cell0 = m0.getCell();
  const auto cell1 = m1.getCell();
  const AtomMatrix pos0 = m0.getPositionsCopy();
  const AtomMatrix pos1 = m1.getPositionsCopy();

  AtomMatrix singleF0 = AtomMatrix::Zero(nAtoms, 3);
  AtomMatrix singleF1 = AtomMatrix::Zero(nAtoms, 3);
  double singleU0{0}, singleU1{0}, var{0};
  pot->force(nAtoms, pos0.data(), nrs0.data(), singleF0.data(), &singleU0, &var,
             cell0.data());
  pot->force(nAtoms, pos1.data(), nrs1.data(), singleF1.data(), &singleU1, &var,
             cell1.data());
  // The displacement has to matter, or the comparison below is vacuous.
  REQUIRE(singleU0 != Catch::Approx(singleU1));

  AtomMatrix batchF0 = AtomMatrix::Zero(nAtoms, 3);
  AtomMatrix batchF1 = AtomMatrix::Zero(nAtoms, 3);
  const double *positions[2] = {pos0.data(), pos1.data()};
  const int *atomicNrs[2] = {nrs0.data(), nrs1.data()};
  double *forces[2] = {batchF0.data(), batchF1.data()};
  const double *boxes[2] = {cell0.data(), cell1.data()};
  double energies[2] = {0.0, 0.0};
  double variances[2] = {0.0, 0.0};
  pot->forceBatch(2, nAtoms, positions, atomicNrs, forces, energies, variances,
                  boxes);

  REQUIRE(energies[0] == Catch::Approx(singleU0));
  REQUIRE(energies[1] == Catch::Approx(singleU1));
  REQUIRE((batchF0 - singleF0).cwiseAbs().maxCoeff() < 1e-10);
  REQUIRE((batchF1 - singleF1).cwiseAbs().maxCoeff() < 1e-10);
}

// Central differences of the mass-weighted gradient at steps h and h/2.
// The h/2 product carries an O(h^2) truncation of a third of their
// difference, so the difference itself bounds it threefold. Positions on
// a dyadic grid, masses that are powers of four, dyadic direction
// components and a power-of-two step make every displaced coordinate
// exact, so no step rounding enters the differences. LJ is pairwise, so
// a force component sums at most `atoms` pair terms and rounds by at most
// atoms * eps times the sum of their magnitudes, taken from dimers.
TEST_CASE("Adapter LJ mass-weighted Hessian-vector product matches central "
          "differences",
          "[pot][rgpot-adapter]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  REQUIRE(pot->supportsMassWeightedHessianVectorProduct());

  Matter matter(pot, params);
  REQUIRE(eonc::io::io_ok(matter.con2matter(std::string("reactant.con"))));
  const long atoms = matter.numberOfAtoms();
  AtomMatrix x0 = matter.getPositionsCopy();
  for (long i = 0; i < atoms; ++i) {
    matter.setMass(i, std::ldexp(1.0, 2 * static_cast<int>(i % 4)));
    for (long axis = 0; axis < 3; ++axis) {
      x0(i, axis) = std::ldexp(std::round(std::ldexp(x0(i, axis), 24)), -24);
    }
  }
  matter.setPositions(x0);
  matter.setFixed(0, 1);
  matter.setFixed(5, 1);
  VectorXi freeAtoms(atoms - 2);
  for (long i = 0, k = 0; i < atoms; ++i) {
    if (i != 0 && i != 5) {
      freeAtoms(k++) = static_cast<int>(i);
    }
  }
  const long dim = 3 * freeAtoms.size();
  VectorXd v(dim);
  for (long k = 0; k < dim; ++k) {
    v(k) = static_cast<double>(k % 7 - 3) / 8.0;
  }

  const VectorXd analytic =
      pot->massWeightedHessianVectorProduct(matter, v, freeAtoms);
  REQUIRE(analytic.size() == dim);

  const auto numbers = matter.getAtomicNrs();
  const Matrix3d cell = matter.getCell();
  auto gradient = [&](double step) {
    AtomMatrix x = x0;
    for (long k = 0; k < freeAtoms.size(); ++k) {
      const double scale = step / std::sqrt(matter.getMass(freeAtoms(k)));
      for (long axis = 0; axis < 3; ++axis) {
        x(freeAtoms(k), axis) += scale * v(3 * k + axis);
      }
    }
    AtomMatrix f = AtomMatrix::Zero(atoms, 3);
    double u{0}, variance{0};
    pot->force(atoms, x.data(), numbers.data(), f.data(), &u, &variance,
               cell.data());
    VectorXd g(dim);
    for (long k = 0; k < freeAtoms.size(); ++k) {
      const double scale = 1.0 / std::sqrt(matter.getMass(freeAtoms(k)));
      for (long axis = 0; axis < 3; ++axis) {
        g(3 * k + axis) = -f(freeAtoms(k), axis) * scale;
      }
    }
    return g;
  };
  const double h = std::ldexp(1.0, -17);
  const VectorXd coarse = (gradient(h) - gradient(-h)) / (2.0 * h);
  const VectorXd fine = (gradient(h / 2) - gradient(-h / 2)) / h;
  AtomMatrix pairSum = AtomMatrix::Zero(atoms, 3);
  for (long i = 0; i < atoms; ++i) {
    for (long j = i + 1; j < atoms; ++j) {
      AtomMatrix pair(2, 3);
      pair.row(0) = x0.row(i);
      pair.row(1) = x0.row(j);
      const int pairNumbers[2] = {numbers(i), numbers(j)};
      AtomMatrix f = AtomMatrix::Zero(2, 3);
      double u{0}, variance{0};
      pot->force(2, pair.data(), pairNumbers, f.data(), &u, &variance,
                 cell.data());
      pairSum.row(i) += f.row(0).cwiseAbs();
      pairSum.row(j) += f.row(1).cwiseAbs();
    }
  }
  const double sumRounding =
      static_cast<double>(atoms) * std::numeric_limits<double>::epsilon();
  for (long k = 0; k < dim; ++k) {
    const long atom = freeAtoms(k / 3);
    const double rounding = 2.0 * sumRounding * pairSum(atom, k % 3) /
                            (std::sqrt(matter.getMass(atom)) * h);
    INFO("component " << k);
    REQUIRE(std::abs(analytic(k) - fine(k)) <=
            std::abs(coarse(k) - fine(k)) + rounding);
  }
}

TEST_CASE("Adapter kernels without an analytic product report none",
          "[pot][rgpot-adapter]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential =
      PotType::LJCLUSTER;
  auto pot =
      eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  REQUIRE_FALSE(pot->supportsMassWeightedHessianVectorProduct());
  Matter matter(pot, params);
  REQUIRE(eonc::io::io_ok(matter.con2matter(std::string("reactant.con"))));
  const VectorXi freeAtoms = VectorXi::LinSpaced(matter.numberOfAtoms(), 0,
                                                 matter.numberOfAtoms() - 1);
  REQUIRE_THROWS_AS(pot->massWeightedHessianVectorProduct(
                        matter, VectorXd::Ones(3 * freeAtoms.size()), freeAtoms),
                    std::logic_error);
}

} /* namespace tests */
