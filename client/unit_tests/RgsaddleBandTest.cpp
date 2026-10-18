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

#include "../ForEachImage.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/NudgedElasticBand.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "eon/XtsciBand.h"
#include <memory>
#include <rgsaddle.h>
#include <tuple>
#include <vector>

#include <cmath>

namespace eonc {
int xtsciBandSurface(void *user, void *request);
} // namespace eonc

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("default LBFGS NEB matches one xtsci band step", "[neb][rgsaddle]") {
  auto stepped = [](OptType method) {
    Parameters params;
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    ParametersLoadAccess::neb_options(params).opt_method = method;
    ParametersLoadAccess::neb_options(params).image_count = 3;
    ParametersLoadAccess::neb_options(params).max_iterations = 2;
    ParametersLoadAccess::neb_options(params).force_tolerance = 1e-12;
    ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
    ParametersLoadAccess::neb_options(params).climbing_image.ocineb.use_mmf =
        false;
    ParametersLoadAccess::neb_options(params).initialization.method =
        NEBInit::LINEAR;
    ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
    ParametersLoadAccess::optimizer_options(params).xtsci.method = "lbfgs";
    auto pot = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(PotType::LJ, params));
    auto reactant = std::make_shared<Matter>(pot, params);
    auto product = std::make_shared<Matter>(pot, params);
    reactant->con2matter(std::string("reactant.con"));
    product->con2matter(std::string("reactant.con"));
    auto shifted = product->getPositions();
    shifted(0, 0) += 0.5;
    product->setPositions(shifted);
    NudgedElasticBand neb(reactant, product, params, pot);
    static_cast<void>(neb.compute());
    return neb.path[1]->getPositions();
  };
  const AtomMatrix lbfgs = stepped(OptType::LBFGS);
  const AtomMatrix xtsci = stepped(OptType::XTSCI);
  REQUIRE(lbfgs.isApprox(xtsci, 0.0));
  REQUIRE(lbfgs.cwiseAbs().maxCoeff() > 0.0);
}

TEST_CASE("rgsaddle band steps a short LJ path", "[neb][rgsaddle]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::neb_options(params).opt_method = OptType::XTSCI;
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).max_iterations = 3;
  ParametersLoadAccess::neb_options(params).force_tolerance = 0.01;
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
  ParametersLoadAccess::neb_options(params).climbing_image.ocineb.use_mmf =
      false;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::LINEAR;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  ParametersLoadAccess::optimizer_options(params).xtsci.method = "fire";

  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  reactant->con2matter(std::string("reactant.con"));
  product->con2matter(std::string("reactant.con"));
  auto pos = product->getPositions();
  pos(0, 0) += 0.5;
  product->setPositions(pos);

  NudgedElasticBand neb(reactant, product, params, pot);
  const auto status = neb.compute();
  REQUIRE(status != NudgedElasticBand::NEBStatus::INIT);
  for (const auto &image : neb.path) {
    REQUIRE(std::isfinite(image->getPotentialEnergy()));
  }
}

TEST_CASE("rgsaddle band evaluates each moved image once per step",
          "[neb][rgsaddle][force_calls]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::neb_options(params).opt_method = OptType::XTSCI;
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).max_iterations = 4;
  ParametersLoadAccess::neb_options(params).force_tolerance = 1e-8;
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
  ParametersLoadAccess::neb_options(params).climbing_image.ocineb.use_mmf =
      false;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::LINEAR;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.05;
  ParametersLoadAccess::optimizer_options(params).xtsci.method = "fire";
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  reactant->con2matter(std::string("reactant.con"));
  product->con2matter(std::string("reactant.con"));
  auto shifted = product->getPositions();
  shifted(0, 0) += 0.5;
  product->setPositions(shifted);
  const size_t before = pot->forceCallCounter;
  NudgedElasticBand neb(reactant, product, params, pot);
  static_cast<void>(neb.compute());
  const size_t calls = pot->forceCallCounter - before;
  // Two endpoints and the first band update, then at most two
  // evaluations of the interior per FIRE step (rgmin's start point and its
  // trial). The fixed endpoints and the band update after each step come
  // from the images' own caches; resending the band every step cost about
  // 3.2 evaluations of the whole band per step.
  const size_t steps = 4;
  const size_t interior = static_cast<size_t>(neb.numImages);
  CAPTURE(calls);
  REQUIRE(calls <= 2 + interior + 2 * steps * interior);
#if RGSADDLE_ABI_MINOR >= 5
  // From ABI minor 5 the session keeps its last evaluation, so a warm
  // step carries exactly the interior images once.
  REQUIRE(calls <= 2 + interior + steps * interior);
#endif
}

TEST_CASE("rgsaddle band surface serves whole, interior and one-image "
          "requests",
          "[neb][rgsaddle]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::neb_options(params).opt_method = OptType::XTSCI;
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::LINEAR;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  reactant->con2matter(std::string("reactant.con"));
  product->con2matter(std::string("reactant.con"));
  auto shifted = product->getPositions();
  shifted(0, 0) += 0.5;
  product->setPositions(shifted);
  NudgedElasticBand neb(reactant, product, params, pot);
  eonc::XtsciBand band(neb, params);

  const long band_n = neb.numImages + 2;
  const long atoms = neb.atoms;
  const long dof = 3 * atoms;
  // Reference energies of every image, and the positions of every image
  // displaced a little so each request is a new evaluation.
  std::vector<double> all(static_cast<size_t>(band_n * dof));
  std::vector<double> reference(static_cast<size_t>(band_n));
  for (long i = 0; i < band_n; ++i) {
    AtomMatrix p = neb.path[static_cast<size_t>(i)]->getPositions();
    p(1, 1) += 0.01 * static_cast<double>(i + 1);
    Matter probe(*neb.path[static_cast<size_t>(i)]);
    probe.setPositions(p);
    reference[static_cast<size_t>(i)] = probe.getPotentialEnergy();
    std::copy(p.data(), p.data() + dof, all.begin() + i * dof);
  }

  auto request = [&](long rows, long first, uint64_t flags, long image) {
    std::vector<double> e(static_cast<size_t>(rows));
    std::vector<double> g(static_cast<size_t>(rows * dof));
    rgsaddle_surface_request_t req{};
    req.version = RGSADDLE_VERSION_INIT;
    req.flags = flags;
    req.n_images = (flags & RGSADDLE_REQ_ONE_IMAGE) != 0 ? band_n : rows;
    req.n_atoms = atoms;
    req.positions = all.data() + first * dof;
    req.energies = e.data();
    req.gradients = g.data();
    req.image = image;
    const int rc = eonc::xtsciBandSurface(&band, &req);
    return std::make_pair(rc, e);
  };

  // The whole band: row r is image r.
  auto [rcAll, eAll] = request(band_n, 0, 0, -1);
  REQUIRE(rcAll == RGSADDLE_OK);
  for (long r = 0; r < band_n; ++r) {
    REQUIRE(eAll[static_cast<size_t>(r)] ==
            Catch::Approx(reference[static_cast<size_t>(r)]));
  }
  // Interior only: row r is image r + 1.
  auto [rcIn, eIn] = request(band_n - 2, 1, 0, -1);
  REQUIRE(rcIn == RGSADDLE_OK);
  for (long r = 0; r < band_n - 2; ++r) {
    REQUIRE(eIn[static_cast<size_t>(r)] ==
            Catch::Approx(reference[static_cast<size_t>(r + 1)]));
  }
  // One image with its band index.
  auto [rcOne, eOne] = request(1, 2, RGSADDLE_REQ_ONE_IMAGE, 2);
  REQUIRE(rcOne == RGSADDLE_OK);
  REQUIRE(eOne[0] == Catch::Approx(reference[2]));
  // Any other row count is a shape error, not a silent mismatch.
  auto [rcBad, eBad] = request(band_n - 1, 0, 0, -1);
  REQUIRE(rcBad == RGSADDLE_SHAPE);
}

TEST_CASE("rgsaddle band surface evaluates its images on the image pool",
          "[neb][rgsaddle][parallel]") {
  using eonc::helpers::test::PerImageLJ;
  auto run = [](bool parallel, std::shared_ptr<PerImageLJ::Threads> threads) {
    Parameters params;
    ParametersLoadAccess::main_options(params).parallel = parallel;
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    ParametersLoadAccess::neb_options(params).opt_method = OptType::FIRE;
    ParametersLoadAccess::neb_options(params).image_count = 7;
    ParametersLoadAccess::neb_options(params).max_iterations = 20;
    ParametersLoadAccess::neb_options(params).force_tolerance = 1e-12;
    ParametersLoadAccess::neb_options(params).climbing_image.enabled = true;
    ParametersLoadAccess::neb_options(params).initialization.method =
        NEBInit::LINEAR;
    ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
    ParametersLoadAccess::optimizer_options(params).max_move = 0.05;
    PotRegistry local;
    auto pot = std::make_shared<PerImageLJ>(threads, local);
    auto reactant = std::make_shared<Matter>(pot, params);
    auto product = std::make_shared<Matter>(pot, params);
    reactant->con2matter(std::string("reactant.con"));
    product->con2matter(std::string("reactant.con"));
    auto shifted = product->getPositions();
    shifted(0, 0) += 0.5;
    shifted(0, 1) -= 0.3;
    product->setPositions(shifted);
    NudgedElasticBand neb(reactant, product, params, pot);
    static_cast<void>(neb.compute());
    std::vector<AtomMatrix> positions;
    std::vector<double> energies;
    for (long i = 0; i <= neb.numImages + 1; ++i) {
      positions.push_back(neb.path[static_cast<size_t>(i)]->getPositions());
      energies.push_back(
          neb.path[static_cast<size_t>(i)]->getPotentialEnergy());
    }
    return std::make_tuple(positions, energies, local.total_force_calls(),
                           neb.lastIteration());
  };
  auto serialThreads = std::make_shared<PerImageLJ::Threads>();
  auto poolThreads = std::make_shared<PerImageLJ::Threads>();
  const auto [sPos, sE, sCalls, sIter] = run(false, serialThreads);
  const auto [pPos, pE, pCalls, pIter] = run(true, poolThreads);
  REQUIRE(sIter == pIter);
  REQUIRE(sIter > 0);
  // Two endpoints, then at least the interior once per step: the
  // registry of the job's potential sees every per-image instance.
  REQUIRE(sCalls == pCalls);
  REQUIRE(sCalls >= 2 + 7 * static_cast<size_t>(sIter));
  for (size_t i = 0; i < sPos.size(); ++i) {
    REQUIRE(sPos[i] == pPos[i]);
    REQUIRE(sE[i] == pE[i]);
  }
  REQUIRE(serialThreads->ids.size() == 1);
  if (eonc::detail::ImagePool::instance().threads() > 0) {
    REQUIRE(poolThreads->ids.size() > 1);
  }
}

} // namespace tests
