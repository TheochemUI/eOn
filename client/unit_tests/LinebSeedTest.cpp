/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#include "eon/LinebSeed.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"
#include "eon/Potential.h"

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {

std::shared_ptr<Matter> loadReactant(Parameters &params) {
  ParametersLoadAccess::main_options(params).parallel = false;
  ParametersLoadAccess::potential_options(params).potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::MORSE_PT, params));
  auto m = std::make_shared<Matter>(pot, params);
  REQUIRE(eonc::io::io_ok(m->con2matter(std::string("reactant.con"))));
  return m;
}

std::vector<std::shared_ptr<Matter>>
threeImageSeed(const Matter &reactant, const AtomMatrix &midShift) {
  std::vector<std::shared_ptr<Matter>> band;
  band.push_back(std::make_shared<Matter>(reactant));
  auto mid = std::make_shared<Matter>(reactant);
  mid->setPositions(reactant.getPositions() + midShift);
  band.push_back(mid);
  auto end = std::make_shared<Matter>(reactant);
  end->setPositions(reactant.getPositions() + 2.0 * midShift);
  band.push_back(end);
  return band;
}

} // namespace

TEST_CASE("lineb seed refuses species identity mismatch", "[lineb][seed]") {
  Parameters params;
  auto ref = loadReactant(params);
  auto image = std::make_shared<Matter>(*ref);
  REQUIRE(image->numberOfAtoms() >= 1);
  image->setAtomicNr(0, ref->getAtomicNr(0) + 1);
  REQUIRE_FALSE(eonc::lineb::matchingSeedImage(*ref, *image));
  REQUIRE_THROWS_WITH(eonc::lineb::requireMatchingSeedImage(*ref, *image),
                      Catch::Matchers::ContainsSubstring("species"));
}

TEST_CASE("lineb seed refuses mass mismatch", "[lineb][seed]") {
  Parameters params;
  auto ref = loadReactant(params);
  auto image = std::make_shared<Matter>(*ref);
  image->setMass(0, ref->getMass(0) + 1.0);
  REQUIRE_THROWS_WITH(eonc::lineb::requireMatchingSeedMasses(*ref, *image),
                      Catch::Matchers::ContainsSubstring("masses"));
}

TEST_CASE("lineb seed refuses moved frozen substrate", "[lineb][seed]") {
  Parameters params;
  auto ref = loadReactant(params);
  REQUIRE(ref->numberOfAtoms() >= 2);
  ref->setFixed(0, 1);
  auto image = std::make_shared<Matter>(*ref);
  AtomMatrix pos = image->getPositions();
  pos.row(0) += Eigen::RowVector3d(0.2, 0.0, 0.0);
  image->setPositions(pos);
  REQUIRE_FALSE(eonc::lineb::matchingSeedImage(*ref, *image));
  REQUIRE_THROWS_WITH(eonc::lineb::requireMatchingSeedImage(*ref, *image),
                      Catch::Matchers::ContainsSubstring("frozen-substrate"));
}

TEST_CASE("lineb seed refuses non-finite positions when packing",
          "[lineb][seed]") {
  Parameters params;
  auto ref = loadReactant(params);
  AtomMatrix shift = AtomMatrix::Zero(ref->numberOfAtoms(), 3);
  shift.col(0).setConstant(0.05);
  auto band = threeImageSeed(*ref, shift);
  AtomMatrix bad = band[1]->getPositions();
  bad(0, 0) = std::numeric_limits<double>::quiet_NaN();
  // Bypass Matter::setPositions validation path by writing through a copy
  // that still reports matching masks: setPositions accepts NaN then packs.
  band[1]->setPositions(bad);
  std::vector<long> freeAtoms;
  std::vector<double> masses;
  eonc::lineb::collectFreeAtoms(*ref, freeAtoms, masses);
  std::vector<double> packed;
  REQUIRE_THROWS_WITH(
      eonc::lineb::packSeedCartesian(band, freeAtoms, false, packed),
      Catch::Matchers::ContainsSubstring("non-finite"));
}

TEST_CASE("lineb pack recovers continuous free coords across a PBC crossing",
          "[lineb][seed]") {
  Parameters params;
  auto ref = loadReactant(params);
  REQUIRE(ref->getPeriodic());
  const double lx = ref->getCell()(0, 0);
  REQUIRE(lx > 1.0);

  // Continuous free path: 0, +0.4 Lx, +0.85 Lx. Matter::setPositions rewraps
  // the last image into the cell; packing must still emit the continuous
  // Cartesian end without mutating Matter.
  AtomMatrix midShift = AtomMatrix::Zero(ref->numberOfAtoms(), 3);
  midShift.col(0).setConstant(0.4 * lx);
  auto band = threeImageSeed(*ref, midShift);
  AtomMatrix endShift = AtomMatrix::Zero(ref->numberOfAtoms(), 3);
  endShift.col(0).setConstant(0.85 * lx);
  const AtomMatrix continuousEnd = ref->getPositions() + endShift;
  band[2]->setPositions(continuousEnd);
  const AtomMatrix wrappedEnd = band[2]->getPositions();
  REQUIRE((wrappedEnd - continuousEnd).cwiseAbs().maxCoeff() > 0.5 * lx - 1e-6);

  std::vector<long> freeAtoms;
  std::vector<double> masses;
  eonc::lineb::collectFreeAtoms(*ref, freeAtoms, masses);
  std::vector<double> packed;
  eonc::lineb::packSeedCartesian(band, freeAtoms, false, packed);

  REQUIRE((band[2]->getPositions() - wrappedEnd).cwiseAbs().maxCoeff() < 1e-12);

  const size_t nFree = freeAtoms.size();
  REQUIRE(packed.size() == 3 * nFree * 3);
  for (size_t a = 0; a < nFree; ++a) {
    const long atom = freeAtoms[a];
    const size_t base = (2 * nFree + a) * 3;
    REQUIRE(packed[base + 0] ==
            Catch::Approx(continuousEnd(atom, 0)).margin(1e-8));
    REQUIRE(packed[base + 1] ==
            Catch::Approx(continuousEnd(atom, 1)).margin(1e-8));
    REQUIRE(packed[base + 2] ==
            Catch::Approx(continuousEnd(atom, 2)).margin(1e-8));
  }
}

TEST_CASE("lineb pack admits endpoints without potential gating",
          "[lineb][seed]") {
  Parameters params;
  auto ref = loadReactant(params);
  AtomMatrix shift = AtomMatrix::Zero(ref->numberOfAtoms(), 3);
  shift.col(0).setConstant(0.05);
  auto band = threeImageSeed(*ref, shift);
  // Endpoints may sit at or below E; packing must not evaluate or gate on V.
  std::vector<long> freeAtoms;
  std::vector<double> masses;
  eonc::lineb::collectFreeAtoms(*ref, freeAtoms, masses);
  std::vector<double> packed;
  REQUIRE_NOTHROW(
      eonc::lineb::packSeedCartesian(band, freeAtoms, false, packed));
  REQUIRE(packed.size() == band.size() * freeAtoms.size() * 3);
  for (size_t a = 0; a < freeAtoms.size(); ++a) {
    const long atom = freeAtoms[a];
    REQUIRE(packed[a * 3 + 0] ==
            Catch::Approx(ref->getPositions()(atom, 0)).margin(1e-12));
  }
}

TEST_CASE("lineb pack refuses periodic molecular seed", "[lineb][seed]") {
  Parameters params;
  auto ref = loadReactant(params);
  REQUIRE(ref->getPeriodic());
  for (long i = 0; i < ref->numberOfAtoms(); ++i) {
    ref->setFixed(i, 0);
  }
  AtomMatrix shift = AtomMatrix::Zero(ref->numberOfAtoms(), 3);
  shift.col(0).setConstant(0.05);
  auto band = threeImageSeed(*ref, shift);
  std::vector<long> freeAtoms;
  std::vector<double> masses;
  eonc::lineb::collectFreeAtoms(*ref, freeAtoms, masses);
  std::vector<double> packed;
  REQUIRE_THROWS_WITH(
      eonc::lineb::packSeedCartesian(band, freeAtoms, true, packed),
      Catch::Matchers::ContainsSubstring("periodic"));
}

TEST_CASE("lineb pack uses raw Cartesian when non-periodic", "[lineb][seed]") {
  Parameters params;
  auto ref = loadReactant(params);
  ref->setPeriodic(false);
  AtomMatrix midShift = AtomMatrix::Zero(ref->numberOfAtoms(), 3);
  midShift.col(0).setConstant(1.25);
  auto band = threeImageSeed(*ref, midShift);
  std::vector<long> freeAtoms;
  std::vector<double> masses;
  eonc::lineb::collectFreeAtoms(*ref, freeAtoms, masses);
  std::vector<double> packed;
  eonc::lineb::packSeedCartesian(band, freeAtoms, false, packed);
  const size_t nFree = freeAtoms.size();
  for (size_t a = 0; a < nFree; ++a) {
    const long atom = freeAtoms[a];
    const size_t base = (2 * nFree + a) * 3;
    REQUIRE(packed[base + 0] ==
            Catch::Approx(ref->getPositions()(atom, 0) + 2.5).margin(1e-10));
  }
}

} // namespace tests
