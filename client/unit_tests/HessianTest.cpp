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

#include "Hessian.h"
#include "Matter.h"
#include "Parameters.h"
#include "SafeMath.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

#include <cstdio>
#include <fstream>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("Hessian on LJ cluster is symmetric", "[hessian]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  REQUIRE(matter->numberOfAtoms() == 13);

  // Use a small subset (3 atoms) to keep the test fast
  VectorXi subAtoms(3);
  subAtoms << 0, 1, 2;

  Hessian hess(params, matter.get());
  MatrixXd H = hess.getHessian(matter.get(), subAtoms);

  REQUIRE(H.rows() == 9);
  REQUIRE(H.cols() == 9);

  for (long i = 0; i < H.rows(); i++) {
    for (long j = 0; j < i; j++) {
      REQUIRE_THAT(H(i, j), Catch::Matchers::WithinAbs(H(j, i), 1e-6));
    }
  }
}

TEST_CASE("Hessian getFreqs returns finite eigenvalues", "[hessian]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  VectorXi subAtoms(3);
  subAtoms << 0, 1, 2;

  Hessian hess(params, matter.get());
  VectorXd freqs = hess.getFreqs(matter.get(), subAtoms);

  REQUIRE(freqs.size() == 9);
  for (long i = 0; i < freqs.size(); i++) {
    REQUIRE(std::isfinite(freqs(i)));
  }
}

TEST_CASE("Hessian getFreqs rejects out-of-range atom indices", "[hessian]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  VectorXi bad(2);
  bad << 0, 99; // 99 is OOB for 13-atom LJ cluster

  Hessian hess(params, matter.get());
  VectorXd freqs = hess.getFreqs(matter.get(), bad);
  REQUIRE(freqs.size() == 0);
}

TEST_CASE("Hessian mobile atom_list yields 3*n_mobile square matrix",
          "[hessian]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  // Mobile/displaced set = hybrid/PHVA-class active list
  VectorXi subAtoms(2);
  subAtoms << 0, 1;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  Hessian hess(params, matter.get());
  MatrixXd H = hess.getHessian(matter.get(), subAtoms);
  REQUIRE(H.rows() == 6);
  REQUIRE(H.cols() == 6);
}

TEST_CASE("Hessian column checkpoint resume matches full FD", "[hessian]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  params.hessian_options.fd_scheme = "one_sided";
  const std::string ckpt = "hessian_resume_test.ckpt";
  std::remove(ckpt.c_str());

  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  VectorXi subAtoms(2);
  subAtoms << 0, 1;

  // Full run writes columns + symmetrizes; also writes ckpt each column then
  // deletes it on success.
  params.hessian_options.resume = false;
  params.hessian_options.checkpoint_path = ckpt;
  Hessian hessFull(params, matter.get());
  MatrixXd Hfull = hessFull.getHessian(matter.get(), subAtoms);
  REQUIRE(Hfull.rows() == 6);
  REQUIRE(!std::ifstream(ckpt).good()); // cleared on success

  // Mid-run interrupt simulation: run with resume+ckpt, but pre-seed a
  // checkpoint that has *unsymmetrized* FD rows 0..1 only. Build those rows
  // by a throwaway full run that keeps the ckpt by using resume=false and
  // manually constructing the file from a controlled partial via second
  // Hessian that we interrupt by writing next_col=2 ourselves after one
  // successful full run's logic: use forces FD for cols 0-1 only (same as
  // Hessian.cpp one_sided) so resume continues cols 2..5 then symmetrizes.
  {
    Matter matterTemp(*matter);
    const double dr = params.main_options.finiteDifference;
    const int nAtoms = matter->numberOfAtoms();
    const int size = 6;
    AtomMatrix pos = matter->getPositions();
    AtomMatrix posDisplace = AtomMatrix::Zero(nAtoms, 3);
    AtomMatrix force0 = matterTemp.getForces();
    MatrixXd Hpart = MatrixXd::Zero(size, size);
    for (int i = 0; i < 2; ++i) {
      posDisplace.setZero();
      posDisplace(subAtoms(i / 3), i % 3) = dr;
      matterTemp.setPositions(pos + posDisplace);
      AtomMatrix forcePlus = matterTemp.getForces();
      for (int j = 0; j < size; ++j) {
        const double dF = forcePlus(subAtoms(j / 3), j % 3) -
                          force0(subAtoms(j / 3), j % 3);
        Hpart(i, j) = -dF / dr;
        const double effMass =
            std::sqrt(matter->getMass(subAtoms(j / 3)) *
                      matter->getMass(subAtoms(i / 3)));
        Hpart(i, j) = eonc::safemath::safe_div(Hpart(i, j), effMass, 0.0);
      }
    }
    matterTemp.setPositions(pos);
    std::ofstream out(ckpt);
    out << "eon_hess_ckpt " << size << " 2\n";
    out.precision(17);
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        out << Hpart(i, j) << (j + 1 == size ? '\n' : ' ');
      }
    }
  }

  params.hessian_options.resume = true;
  params.hessian_options.checkpoint_path = ckpt;
  Hessian hessRes(params, matter.get());
  MatrixXd Hres = hessRes.getHessian(matter.get(), subAtoms);

  REQUIRE(Hres.rows() == Hfull.rows());
  for (long i = 0; i < Hfull.rows(); ++i) {
    for (long j = 0; j < Hfull.cols(); ++j) {
      REQUIRE_THAT(Hres(i, j), Catch::Matchers::WithinAbs(Hfull(i, j), 1e-6));
    }
  }
  REQUIRE(!std::ifstream(ckpt).good());
}

TEST_CASE("Hessian central fd_scheme produces finite symmetric H",
          "[hessian]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  params.hessian_options.fd_scheme = "central";
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  VectorXi subAtoms(3);
  subAtoms << 0, 1, 2;

  Hessian hess(params, matter.get());
  MatrixXd H = hess.getHessian(matter.get(), subAtoms);
  REQUIRE(H.rows() == 9);
  for (long i = 0; i < H.rows(); ++i) {
    for (long j = 0; j < i; ++j) {
      REQUIRE_THAT(H(i, j), Catch::Matchers::WithinAbs(H(j, i), 1e-6));
    }
    REQUIRE(std::isfinite(H(i, i)));
  }
}

TEST_CASE("Hessian on Pt frozen layers system handles mixed fixed/free",
          "[hessian][morse_pt]") {
  Parameters params;
  params.potential_options.potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("../Pt_Heptamer_FrozenLayers/pos.con"));

  // Build moved atom list (free atoms only)
  long nAtoms = matter->numberOfAtoms();
  std::vector<long> freeList;
  for (long i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i))
      freeList.push_back(i);
  }
  REQUIRE(freeList.size() > 0);
  REQUIRE(freeList.size() < static_cast<size_t>(nAtoms));

  VectorXi moved(freeList.size());
  for (size_t i = 0; i < freeList.size(); i++) {
    moved[i] = freeList[i];
  }

  Hessian hess(params, matter.get());
  VectorXd freqs = hess.getFreqs(matter.get(), moved);

  // SVN crashes here (VectorXi size mismatch). Our fix handles it.
  long expectedSize = static_cast<long>(freeList.size()) * 3;
  REQUIRE(freqs.size() == expectedSize);
  for (long i = 0; i < freqs.size(); i++) {
    REQUIRE(std::isfinite(freqs(i)));
  }
}

} /* namespace tests */
