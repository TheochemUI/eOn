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
#include <algorithm>
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

TEST_CASE("rgsaddle band keeps a zero force on fixed atoms",
          "[neb][rgsaddle]") {
  // eOn drops the mean force of an image only when every atom is free.
  // Removing it with fixed atoms present pushes those atoms; the host pins
  // them back after every step, and the L-BFGS band then records that
  // correction as a curvature pair and cycles at max_move.
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::neb_options(params).opt_method = OptType::LBFGS;
  ParametersLoadAccess::neb_options(params).image_count = 3;
  ParametersLoadAccess::neb_options(params).max_iterations = 3;
  ParametersLoadAccess::neb_options(params).force_tolerance = 1e-12;
  ParametersLoadAccess::neb_options(params).climbing_image.enabled = false;
  ParametersLoadAccess::neb_options(params).climbing_image.ocineb.use_mmf =
      false;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::LINEAR;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  reactant->con2matter(std::string("reactant.con"));
  product->con2matter(std::string("reactant.con"));
  const long nAtoms = reactant->numberOfAtoms();
  for (long j = nAtoms / 2; j < nAtoms; ++j) {
    reactant->setFixed(j, 1);
    product->setFixed(j, 1);
  }
  auto shifted = product->getPositions();
  shifted(0, 0) += 0.5;
  product->setPositions(shifted);
  NudgedElasticBand neb(reactant, product, params, pot);
  static_cast<void>(neb.compute());
  double freeForce = 0.0;
  for (long i = 1; i <= neb.numImages; ++i) {
    const AtomMatrix &force = *neb.projectedForce[static_cast<size_t>(i)];
    for (long j = 0; j < nAtoms; ++j) {
      CAPTURE(i, j);
      if (neb.path[static_cast<size_t>(i)]->getFixed(j)) {
        REQUIRE(force.row(j).cwiseAbs().maxCoeff() == 0.0);
      } else {
        freeForce = std::max(freeForce, force.row(j).cwiseAbs().maxCoeff());
      }
    }
  }
  REQUIRE(freeForce > 0.0);
}

TEST_CASE("rgsaddle band hands its highest image and tangents to the host",
          "[neb][rgsaddle][climbing_image]") {
  // The session projects the force, so the host never runs its own
  // projection while the band steps; sp.con (the highest image), the
  // climbing image and neb.dat (the tangents) must still describe the
  // accepted band, not the first one.
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::neb_options(params).opt_method = OptType::LBFGS;
  ParametersLoadAccess::neb_options(params).image_count = 5;
  ParametersLoadAccess::neb_options(params).max_iterations = 15;
  ParametersLoadAccess::neb_options(params).force_tolerance = 1e-12;
  ParametersLoadAccess::neb_options(params).climbing_image.ocineb.use_mmf =
      false;
  ParametersLoadAccess::neb_options(params).initialization.method =
      NEBInit::LINEAR;
  ParametersLoadAccess::neb_options(params).endpoints.minimize = false;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
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
  REQUIRE(neb.lastIteration() > 1);

  long highest = 1;
  for (long i = 2; i <= neb.numImages; ++i) {
    if (neb.path[static_cast<size_t>(i)]->getPotentialEnergy() >
        neb.path[static_cast<size_t>(highest)]->getPotentialEnergy()) {
      highest = i;
    }
  }
  REQUIRE(static_cast<long>(neb.maxEnergyImage) == highest);
  REQUIRE(neb.climbingImage >= 0);
  REQUIRE(neb.climbingImage <= neb.numImages);

  std::vector<AtomMatrix> held;
  for (long i = 1; i <= neb.numImages; ++i) {
    held.push_back(*neb.tangent[static_cast<size_t>(i)]);
  }
  // Every image holds its evaluation, so this projection costs no force
  // call and rebuilds the tangents from the final band.
  const size_t before = pot->forceCallCounter;
  neb.updateForces();
  REQUIRE(pot->forceCallCounter == before);
  for (long i = 1; i <= neb.numImages; ++i) {
    CAPTURE(i);
    REQUIRE(held[static_cast<size_t>(i - 1)].isApprox(
        *neb.tangent[static_cast<size_t>(i)], 1e-12));
  }
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

// Three interior LJ images on the neb_morse cluster, one atom moved 0.5.
struct LJBand {
  Parameters params;
  std::shared_ptr<Potential> pot;

  LJBand(OptType method, std::string precon = "none",
         std::string stepper = "auto") {
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    auto &neb = ParametersLoadAccess::neb_options(params);
    neb.opt_method = method;
    neb.image_count = 3;
    neb.max_iterations = 3;
    neb.force_tolerance = 1e-12;
    neb.climbing_image.enabled = false;
    neb.climbing_image.ocineb.use_mmf = false;
    neb.initialization.method = NEBInit::LINEAR;
    neb.endpoints.minimize = false;
    neb.band.precon = std::move(precon);
    neb.band.stepper = std::move(stepper);
    ParametersLoadAccess::optimizer_options(params).max_move = 0.2;
  }

  std::unique_ptr<NudgedElasticBand> neb(long atoms = 0) {
    pot = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(PotType::LJ, params));
    auto reactant = std::make_shared<Matter>(pot, params);
    auto product = std::make_shared<Matter>(pot, params);
    reactant->con2matter(std::string("reactant.con"));
    product->con2matter(std::string("reactant.con"));
    if (atoms > 0) {
      const AtomMatrix kept = reactant->getPositions().topRows(atoms);
      for (auto *m : {reactant.get(), product.get()}) {
        m->resize(atoms);
        m->setAtomicNr(0, 1);
        m->setPositions(kept);
      }
    }
    auto shifted = product->getPositions();
    shifted(0, 0) += 0.5;
    product->setPositions(shifted);
    return std::make_unique<NudgedElasticBand>(reactant, product, params, pot);
  }

  // The middle interior image after three steps.
  AtomMatrix run() {
    auto band = neb();
    static_cast<void>(band->compute());
    REQUIRE(std::isfinite(band->path[2]->getPotentialEnergy()));
    return band->path[2]->getPositions();
  }
};

TEST_CASE("band_precon and band_stepper default off and refuse other values",
          "[neb][rgsaddle][band_precon]") {
  using Catch::Matchers::ContainsSubstring;
  LJBand off(OptType::LBFGS);
  REQUIRE(off.params.neb_options().band.precon == "none");
  REQUIRE(off.params.neb_options().band.stepper == "auto");
  auto neb = off.neb();
  ParametersLoadAccess::neb_options(off.params).band.precon = "pair";
  REQUIRE_THROWS_WITH(eonc::XtsciBand(*neb, off.params),
                      ContainsSubstring("band_precon is none or exp"));
  ParametersLoadAccess::neb_options(off.params).band.precon = "none";
  ParametersLoadAccess::neb_options(off.params).band.stepper = "rk4";
  REQUIRE_THROWS_WITH(eonc::XtsciBand(*neb, off.params),
                      ContainsSubstring("band_stepper is auto or ode12r"));
  // An eOn optimizer does not step the rgsaddle band.
  REQUIRE_THROWS_WITH(LJBand(OptType::CG, "exp").neb()->compute(),
                      ContainsSubstring("act on the rgsaddle band"));
}

#ifndef WITH_RGSADDLE_BAND_PRECON
TEST_CASE("an rgsaddle older than 0.3.0 refuses band_precon and ode12r",
          "[neb][rgsaddle][band_precon]") {
  using Catch::Matchers::ContainsSubstring;
  REQUIRE_THROWS_WITH(LJBand(OptType::LBFGS, "exp").neb()->compute(),
                      ContainsSubstring("rgsaddle 0.3.0 or newer"));
  REQUIRE_THROWS_WITH(LJBand(OptType::FIRE, "none", "ode12r").neb()->compute(),
                      ContainsSubstring("rgsaddle 0.3.0 or newer"));
}
#else
TEST_CASE("the exp metric and ode12r each change the rgsaddle band step",
          "[neb][rgsaddle][band_precon]") {
  const AtomMatrix lbfgs = LJBand(OptType::LBFGS).run();
  const AtomMatrix fire = LJBand(OptType::FIRE).run();
  const AtomMatrix ode = LJBand(OptType::FIRE, "none", "ode12r").run();
  REQUIRE_FALSE(lbfgs.isApprox(LJBand(OptType::LBFGS, "exp").run(), 1e-10));
  REQUIRE_FALSE(fire.isApprox(LJBand(OptType::FIRE, "exp").run(), 1e-10));
  REQUIRE_FALSE(fire.isApprox(ode, 1e-10));
  REQUIRE_FALSE(
      ode.isApprox(LJBand(OptType::FIRE, "exp", "ode12r").run(), 1e-10));
}

TEST_CASE("the exp metric refuses a non-uniform spring and one atom",
          "[neb][rgsaddle][band_precon]") {
  using Catch::Matchers::ContainsSubstring;
  LJBand weighted(OptType::LBFGS, "exp");
  ParametersLoadAccess::neb_options(weighted.params).spring.weighting.enabled =
      true;
  REQUIRE_THROWS_WITH(weighted.neb()->compute(),
                      ContainsSubstring("uniform spring"));
  LJBand om(OptType::LBFGS, "exp");
  ParametersLoadAccess::neb_options(om.params).spring.om.enabled = true;
  REQUIRE_THROWS_WITH(om.neb()->compute(), ContainsSubstring("uniform spring"));
  LJBand single(OptType::LBFGS, "exp");
  auto neb = single.neb(1);
  REQUIRE_THROWS_WITH(eonc::XtsciBand(*neb, single.params),
                      ContainsSubstring("two atoms"));
}

TEST_CASE("the exp metric refuses fixed atoms",
          "[neb][rgsaddle][band_precon]") {
  using Catch::Matchers::ContainsSubstring;
  LJBand fixture(OptType::LBFGS, "exp");
  auto neb = fixture.neb();
  for (auto &image : neb->path) {
    image->setFixed(1, 1);
  }
  REQUIRE_THROWS_WITH(eonc::XtsciBand(*neb, fixture.params),
                      ContainsSubstring("refused with fixed atoms"));
  auto &band = ParametersLoadAccess::neb_options(fixture.params).band;
  band.precon = "none";
  band.stepper = "ode12r";
  REQUIRE_NOTHROW(eonc::XtsciBand(*neb, fixture.params));
}
#endif

TEST_CASE("a one-row band request evaluates off the band",
          "[neb][rgsaddle][band_precon]") {
  LJBand fixture(OptType::LBFGS);
  auto neb = fixture.neb();
  eonc::XtsciBand band(*neb, fixture.params);
  const AtomMatrix held = neb->path[1]->getPositions();
  AtomMatrix p = held;
  p(1, 1) += 0.01;
  Matter probe(*neb->path[1]);
  probe.setPositions(p);
  std::vector<double> e(1);
  std::vector<double> g(static_cast<size_t>(3 * neb->atoms));
  rgsaddle_surface_request_t req{};
  req.version = RGSADDLE_VERSION_INIT;
  req.n_images = 1;
  req.n_atoms = neb->atoms;
  req.positions = p.data();
  req.energies = e.data();
  req.gradients = g.data();
  req.image = -1;
  REQUIRE(eonc::xtsciBandSurface(&band, &req) == RGSADDLE_OK);
  REQUIRE(e[0] == Catch::Approx(probe.getPotentialEnergy()));
  REQUIRE(g[0] == Catch::Approx(-probe.getForces().data()[0]));
  REQUIRE(neb->path[1]->getPositions() == held);
}

TEST_CASE("a one-row band request with an image or flags is a shape error",
          "[neb][rgsaddle][band_precon]") {
  LJBand fixture(OptType::LBFGS);
  auto neb = fixture.neb();
  eonc::XtsciBand band(*neb, fixture.params);
  const AtomMatrix held = neb->path[1]->getPositions();
  AtomMatrix p = held;
  p(1, 1) += 0.01;
  auto request = [&](uint64_t flags, int64_t image) {
    std::vector<double> e(1, 0.0);
    std::vector<double> g(static_cast<size_t>(3 * neb->atoms), 0.0);
    rgsaddle_surface_request_t req{};
    req.version = RGSADDLE_VERSION_INIT;
    req.flags = flags;
    req.n_images = 1;
    req.n_atoms = neb->atoms;
    req.positions = p.data();
    req.energies = e.data();
    req.gradients = g.data();
    req.image = image;
    const int rc = eonc::xtsciBandSurface(&band, &req);
    return std::make_pair(rc, e[0]);
  };
  const std::pair<uint64_t, int64_t> refused[] = {
      {0, 0}, {0, 1}, {RGSADDLE_REQ_ONE_IMAGE << 1, -1}};
  for (const auto &[flags, image] : refused) {
    const auto [rc, energy] = request(flags, image);
    REQUIRE(rc == RGSADDLE_SHAPE);
    REQUIRE(energy == 0.0);
  }
  REQUIRE(request(0, -1).first == RGSADDLE_OK);
  REQUIRE(neb->path[1]->getPositions() == held);
}

} // namespace tests
