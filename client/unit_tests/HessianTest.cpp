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

#include "eon/Hessian.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/AtomicFile.h"
#include "eon/ConFileIO.h"
#include "eon/CutoffNeighbors.h"
#include "eon/Davidson.h"
#include "eon/FiniteDifference.h"
#include "eon/Lanczos.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/SafeMath.h"
#include "eon/Tunneling.h"
#include "eon/VesinNeighbors.h"
#include "eon/potentials/RgpotAdapter/RgpotAdapter.h"
#include "rgpot/LennardJones/LJPot.hpp"
#include <array>
#include <map>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

// Meson runs every neb_morse suite in this directory at once. Hessian
// truncates hessian.dat in the current directory, and that path is also a
// tracked fixture other tests copy. A private directory keeps the open off
// that shared file (Windows denies the write while a copy holds it).
class HessianScratch {
protected:
  std::filesystem::path originalDir;
  std::filesystem::path workdir;
  std::filesystem::path ptPos;

  HessianScratch()
      : originalDir(std::filesystem::current_path()) {
    static int counter = 0;
    workdir = std::filesystem::temp_directory_path() /
              ("eon_test_hess_" + std::to_string(counter++));
    std::filesystem::remove_all(workdir);
    std::filesystem::create_directories(workdir);
    std::filesystem::copy_file(originalDir / "reactant.con",
                               workdir / "reactant.con");
    ptPos = originalDir / ".." / "Pt_Heptamer_FrozenLayers" / "pos.con";
    std::filesystem::current_path(workdir);
  }

  ~HessianScratch() {
    std::error_code ec;
    std::filesystem::current_path(originalDir, ec);
    std::filesystem::remove_all(workdir, ec);
  }
};

TEST_CASE_METHOD(HessianScratch, "Hessian on LJ cluster is symmetric",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
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

  std::ifstream saved("hessian.dat");
  MatrixXd restored(H.rows(), H.cols());
  for (long i = 0; i < H.rows(); ++i) {
    for (long j = 0; j < H.cols(); ++j) {
      REQUIRE(static_cast<bool>(saved >> restored(i, j)));
      REQUIRE(restored(i, j) == H(i, j));
    }
  }

  for (long i = 0; i < H.rows(); i++) {
    for (long j = 0; j < i; j++) {
      REQUIRE_THAT(H(i, j), Catch::Matchers::WithinAbs(H(j, i), 1e-6));
    }
  }
}

TEST_CASE_METHOD(HessianScratch, "Hessian getFreqs returns finite eigenvalues",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
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

TEST_CASE_METHOD(HessianScratch,
                 "Hessian getFreqs rejects out-of-range atom indices",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  VectorXi bad(2);
  bad << 0, 99; // 99 is OOB for 13-atom LJ cluster

  Hessian hess(params, matter.get());
  VectorXd freqs = hess.getFreqs(matter.get(), bad);
  REQUIRE(freqs.size() == 0);
}

TEST_CASE_METHOD(HessianScratch,
                 "Hessian mobile phva_atoms yields 3*n_mobile square matrix",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  // Mobile/displaced set = hybrid/PHVA-class active list
  VectorXi subAtoms(2);
  subAtoms << 0, 1;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));

  Hessian hess(params, matter.get());
  MatrixXd H = hess.getHessian(matter.get(), subAtoms);
  REQUIRE(H.rows() == 6);
  REQUIRE(H.cols() == 6);
}

TEST_CASE_METHOD(HessianScratch,
                 "Hessian column checkpoint resume matches full FD",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::hessian_options(params).fd_scheme = "one_sided";
  const std::string ckpt = "hessian_resume_test.ckpt";
  std::remove(ckpt.c_str());

  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  VectorXi subAtoms(2);
  subAtoms << 0, 1;

  // Full run writes columns + symmetrizes; also writes ckpt each column then
  // deletes it on success.
  ParametersLoadAccess::hessian_options(params).resume = false;
  ParametersLoadAccess::hessian_options(params).checkpoint_path = ckpt;
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
    const double dr = params.main_options().finiteDifference;
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
        const double dF =
            forcePlus(subAtoms(j / 3), j % 3) - force0(subAtoms(j / 3), j % 3);
        Hpart(i, j) = -dF / dr;
        const double effMass = std::sqrt(matter->getMass(subAtoms(j / 3)) *
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

  ParametersLoadAccess::hessian_options(params).resume = true;
  ParametersLoadAccess::hessian_options(params).checkpoint_path = ckpt;
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

TEST_CASE_METHOD(HessianScratch,
                 "Hessian central fd_scheme produces finite symmetric H",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::hessian_options(params).fd_scheme = "central";
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
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

TEST_CASE_METHOD(HessianScratch,
                 "Hessian on Pt frozen layers system handles mixed fixed/free",
                 "[hessian][morse_pt]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::MORSE_PT, params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(ptPos.string());

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

TEST_CASE("Trivial modes are counted from the structure's symmetries",
          "[hessian]") {
  REQUIRE(trivialModeCountIsPhysical(6, 0)); // free cluster
  REQUIRE(trivialModeCountIsPhysical(5, 0)); // linear molecule
  REQUIRE(trivialModeCountIsPhysical(3, 0)); // periodic cell
  REQUIRE(trivialModeCountIsPhysical(0, 4)); // fixed atoms pin everything
  REQUIRE_FALSE(trivialModeCountIsPhysical(4, 0));
  REQUIRE_FALSE(trivialModeCountIsPhysical(3, 4));
}

TEST_CASE_METHOD(HessianScratch,
                 "removeZeroFreqs drops the modes under the zero threshold",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(std::string("reactant.con"));
  const long n = 3 * matter->numberOfAtoms();
  const double zero = params.hessian_options().zero_freq_value;
  // A free cluster's spectrum: six modes under the threshold, the rest above.
  VectorXd freqs = VectorXd::LinSpaced(n, 1.0, 2.0);
  for (int i = 0; i < 6; ++i) {
    freqs(5 * i) = 0.1 * zero * (i % 2 == 0 ? 1.0 : -1.0);
  }
  Hessian hess(params, matter.get());
  const VectorXd kept = hess.removeZeroFreqs(freqs);
  REQUIRE(kept.size() == n - 6);
  REQUIRE((kept.array().abs() > zero).all());
  // A spectrum that is not 3N long is not the whole structure's; it is left
  // alone.
  REQUIRE(hess.removeZeroFreqs(freqs.head(n - 1)).size() == n - 1);
}

TEST_CASE("Colored FD Hessian matches serial central difference",
          "[hessian][color]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::main_options(params).finiteDifference = 1e-5;
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  ParametersLoadAccess::hessian_options(params).fd_scheme = "central";

  // Nearest-neighbor chain: cutoff 1.5, spacing 1.2. Atoms 0 and 3 are
  // not adjacent and share a color of the squared cutoff graph.
  constexpr double kCutoff = 1.5;
  auto pot = std::make_shared<RgpotAdapter<rgpot::LJPot>>(
      PotType::LJ, params, rgpot::LJConfig{.cutoff = kCutoff});
  REQUIRE(pot->hessianRange() == kCutoff);

  auto matter = std::make_shared<Matter>(pot, params);
  matter->resize(4);
  for (long i = 0; i < 4; ++i) {
    matter->setAtomicNr(i, 1);
    matter->setMass(i, 1.0);
  }
  AtomMatrix pos(4, 3);
  pos.setZero();
  pos(0, 0) = 0.0;
  pos(1, 0) = 1.2;
  pos(2, 0) = 2.4;
  pos(3, 0) = 3.6;
  matter->setPositions(pos);
  matter->setCell(Matrix3d::Identity() * 40.0);
  matter->setPeriodic(false);

  VectorXi all(4);
  all << 0, 1, 2, 3;

  const std::vector<int> colors =
      eonc::colorMobileCutoffGraph(*matter, all, pot->hessianRange());
  REQUIRE(colors.size() == 4);
  bool shared = false;
  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      if (colors[static_cast<size_t>(i)] != colors[static_cast<size_t>(j)]) {
        continue;
      }
      REQUIRE(matter->distance(i, j) > kCutoff);
      shared = true;
    }
  }
  REQUIRE(shared);
  int nColors = 0;
  for (int c : colors) {
    nColors = std::max(nColors, c + 1);
  }
  REQUIRE(nColors == 3);

  pot->forceCallCounter.store(0);
  Hessian hess(params, matter.get());
  MatrixXd H = hess.getHessian(matter.get(), all);
  REQUIRE(H.rows() == 12);
  // Central difference: two evaluations per direction per color, plus the
  // undisplaced gradient. Strictly below one column per coordinate.
  const auto calls = pot->forceCallCounter.load();
  REQUIRE(calls == static_cast<size_t>(1 + 2 * 3 * nColors));
  REQUIRE(calls < static_cast<size_t>(1 + 2 * 12));

  const double dr = params.main_options().finiteDifference;
  Matter probe(*matter);
  MatrixXd serial = MatrixXd::Zero(12, 12);
  for (int i = 0; i < 12; ++i) {
    AtomMatrix disp = AtomMatrix::Zero(4, 3);
    disp(all(i / 3), i % 3) = dr;
    probe.setPositions(pos + disp);
    AtomMatrix fp = probe.getForces();
    probe.setPositions(pos - disp);
    AtomMatrix fm = probe.getForces();
    for (int j = 0; j < 12; ++j) {
      const double dF = fp(all(j / 3), j % 3) - fm(all(j / 3), j % 3);
      serial(i, j) = -dF / (2.0 * dr);
      const double effMass =
          std::sqrt(matter->getMass(all(j / 3)) * matter->getMass(all(i / 3)));
      serial(i, j) = eonc::safemath::safe_div(serial(i, j), effMass, 0.0);
    }
  }
  for (int i = 0; i < 12; ++i) {
    for (int j = 0; j < i; ++j) {
      serial(i, j) = (serial(i, j) + serial(j, i)) / 2.0;
      serial(j, i) = serial(i, j);
    }
  }

  for (int i = 0; i < 12; ++i) {
    for (int j = 0; j < 12; ++j) {
      REQUIRE_THAT(H(i, j), Catch::Matchers::WithinAbs(serial(i, j), 1e-8));
    }
  }
}

TEST_CASE("Fourth-order stencil matches a complex-step oracle", "[hessian]") {
  // Polynomial oracle. Complex-step stays here; Potential::force is real.
  const auto f = [](std::complex<double> z) {
    return z * z * z * z + std::complex<double>(0.3, 0.0) * z * z;
  };
  const double x = 0.8;
  const double hcs = 1e-8;
  const double oracle = std::imag(f(std::complex<double>(x, hcs))) / hcs;
  const double dr = 0.05;
  const auto sample = [&](double z) {
    VectorXd v(1);
    v(0) = std::real(f(std::complex<double>(z, 0.0)));
    return v;
  };
  const VectorXd f0 = sample(x);
  const VectorXd fourth = fdForceDerivative(
      FdScheme::Fourth, dr, f0, sample(x + dr), sample(x - dr),
      sample(x + 2.0 * dr), sample(x - 2.0 * dr));
  const VectorXd central = fdForceDerivative(
      FdScheme::Central, dr, f0, sample(x + dr), sample(x - dr), f0, f0);
  REQUIRE(parseFdScheme("fourth_order") == FdScheme::Fourth);
  REQUIRE(parseFdScheme("CENTRAL4") == FdScheme::Fourth);
  REQUIRE(parseFdScheme("central") == FdScheme::Central);
  REQUIRE(parseFdScheme("nope") == FdScheme::OneSided);
  REQUIRE_THAT(fourth(0), Catch::Matchers::WithinAbs(oracle, 1e-10));
  REQUIRE(std::abs(central(0) - oracle) > 1e-3);
  REQUIRE(std::abs(fourth(0) - oracle) < std::abs(central(0) - oracle));
}

TEST_CASE("Colored fourth-order FD matches the serial stencil",
          "[hessian][color]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::main_options(params).finiteDifference = 1e-4;
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  ParametersLoadAccess::hessian_options(params).fd_scheme = "fourth";
  ParametersLoadAccess::lanczos_options(params).max_iterations = 3;
  ParametersLoadAccess::davidson_options(params).max_iterations = 3;

  constexpr double kCutoff = 1.5;
  auto pot = std::make_shared<RgpotAdapter<rgpot::LJPot>>(
      PotType::LJ, params, rgpot::LJConfig{.cutoff = kCutoff});
  auto matter = std::make_shared<Matter>(pot, params);
  matter->resize(4);
  for (long i = 0; i < 4; ++i) {
    matter->setAtomicNr(i, 1);
    matter->setMass(i, 1.0);
  }
  AtomMatrix pos(4, 3);
  pos.setZero();
  pos(0, 0) = 0.0;
  pos(1, 0) = 1.2;
  pos(2, 0) = 2.4;
  pos(3, 0) = 3.6;
  matter->setPositions(pos);
  matter->setCell(Matrix3d::Identity() * 40.0);
  matter->setPeriodic(false);

  VectorXi all(4);
  all << 0, 1, 2, 3;
  const std::vector<int> colors =
      eonc::colorMobileCutoffGraph(*matter, all, pot->hessianRange());
  int nColors = 0;
  for (int c : colors) {
    nColors = std::max(nColors, c + 1);
  }
  REQUIRE(nColors == 3);

  pot->forceCallCounter.store(0);
  Hessian hess(params, matter.get());
  MatrixXd H = hess.getHessian(matter.get(), all);
  REQUIRE(H.rows() == 12);
  const auto calls = pot->forceCallCounter.load();
  REQUIRE(calls == static_cast<size_t>(1 + 4 * 3 * nColors));

  const double dr = params.main_options().finiteDifference;
  Matter probe(*matter);
  const AtomMatrix force0 = probe.getForces();
  MatrixXd serial = MatrixXd::Zero(12, 12);
  for (int i = 0; i < 12; ++i) {
    AtomMatrix disp = AtomMatrix::Zero(4, 3);
    disp(all(i / 3), i % 3) = dr;
    probe.setPositions(pos + disp);
    const AtomMatrix fp = probe.getForces();
    probe.setPositions(pos - disp);
    const AtomMatrix fm = probe.getForces();
    disp(all(i / 3), i % 3) = 2.0 * dr;
    probe.setPositions(pos + disp);
    const AtomMatrix fp2 = probe.getForces();
    probe.setPositions(pos - disp);
    const AtomMatrix fm2 = probe.getForces();
    const AtomMatrix slope =
        fdForceDerivative(FdScheme::Fourth, dr, force0, fp, fm, fp2, fm2);
    for (int j = 0; j < 12; ++j) {
      serial(i, j) = -slope(all(j / 3), j % 3);
      const double effMass =
          std::sqrt(matter->getMass(all(j / 3)) * matter->getMass(all(i / 3)));
      serial(i, j) = eonc::safemath::safe_div(serial(i, j), effMass, 0.0);
    }
  }
  for (int i = 0; i < 12; ++i) {
    for (int j = 0; j < i; ++j) {
      serial(i, j) = (serial(i, j) + serial(j, i)) / 2.0;
      serial(j, i) = serial(i, j);
    }
  }
  for (int i = 0; i < 12; ++i) {
    for (int j = 0; j < 12; ++j) {
      REQUIRE_THAT(H(i, j), Catch::Matchers::WithinAbs(serial(i, j), 1e-8));
    }
  }

  AtomMatrix direction = AtomMatrix::Zero(4, 3);
  direction(0, 0) = 1.0;
  Lanczos lanczos(matter, params, pot);
  lanczos.compute(matter, direction);
  REQUIRE(std::isfinite(lanczos.getEigenvalue()));
  Davidson davidson(matter, params, pot);
  davidson.compute(matter, direction);
  REQUIRE(std::isfinite(davidson.getEigenvalue()));
}

TEST_CASE_METHOD(HessianScratch,
                 "modes.con carries each mode as a displacements section",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  Matter matter(pot, params);
  matter.con2matter(std::string("reactant.con"));
  const long n = matter.numberOfAtoms();

  // Two mobile atoms and the unit mass-weighted modes: mode k moves one
  // coordinate, so its Cartesian mode is that coordinate with norm 1.
  VectorXi mobile(2);
  mobile << 1, 4;
  VectorXd eigenvalues(6);
  eigenvalues << -0.5, 1e-9, 0.25, 1.0, 2.0, 4.0;
  const MatrixXd modes = MatrixXd::Identity(6, 6);
  REQUIRE(writeNormalModes(matter, mobile, eigenvalues, modes, "modes.con"));

  const auto frames = readcon::read_all_frames("modes.con");
  REQUIRE(frames.size() == 6);
  for (size_t k = 0; k < frames.size(); ++k) {
    REQUIRE(frames[k].has_displacements());
    std::vector<double> d(static_cast<size_t>(3 * n));
    REQUIRE(frames[k].copy_displacements(d.data(), d.size()) ==
            readcon::RKR_STATUS_SUCCESS);
    const long atom = mobile(static_cast<long>(k) / 3);
    for (long i = 0; i < 3 * n; ++i) {
      const double want =
          (i == 3 * atom + static_cast<long>(k) % 3) ? 1.0 : 0.0;
      REQUIRE_THAT(d[static_cast<size_t>(i)],
                   Catch::Matchers::WithinAbs(want, 1e-12));
    }
    const auto md = nlohmann::json::parse(frames[k].metadata_json());
    const double lambda = eigenvalues(static_cast<long>(k));
    const double hw =
        std::copysign(tunneling::kHbar * std::sqrt(std::abs(lambda)), lambda);
    REQUIRE_THAT(md.at("mode_eigenvalue").get<double>(),
                 Catch::Matchers::WithinRel(lambda, 1e-12));
    REQUIRE_THAT(md.at("hbar_omega").get<double>(),
                 Catch::Matchers::WithinRel(hw, 1e-12));
    REQUIRE_THAT(md.at("wavenumber").get<double>(),
                 Catch::Matchers::WithinRel(hw * kEvToWavenumber, 1e-12));
  }
  // The imaginary mode reads negative.
  REQUIRE(nlohmann::json::parse(frames[0].metadata_json())
              .at("wavenumber")
              .get<double>() < 0.0);
}

TEST_CASE("cartesianMode divides by sqrt(mass) and normalizes", "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  Matter matter(pot, params);
  matter.resize(2);
  matter.setMass(0, 1.0);
  matter.setMass(1, 4.0);
  VectorXi atoms(2);
  atoms << 0, 1;
  VectorXd q = VectorXd::Zero(6);
  q(0) = 1.0; // atom 0, x
  q(3) = 1.0; // atom 1, x
  const auto x = cartesianMode(matter, atoms, q);
  // x = (1/1, 1/2) along x, then unit norm.
  const double norm = std::sqrt(1.0 + 0.25);
  REQUIRE_THAT(x[0], Catch::Matchers::WithinAbs(1.0 / norm, 1e-12));
  REQUIRE_THAT(x[3], Catch::Matchers::WithinAbs(0.5 / norm, 1e-12));
  REQUIRE_THROWS_AS(cartesianMode(matter, atoms, VectorXd::Zero(3)),
                    std::invalid_argument);
}

namespace {
// LJ without a reported cutoff, so neither side takes the colored path;
// Batched reports batch support and takes the batched path through the
// default forceBatch loop.
template <bool Batched> class WrappedLJ final : public eonc::Potential {
public:
  WrappedLJ(std::shared_ptr<eonc::Potential> inner, const Parameters &p)
      : eonc::Potential(PotType::LJ, p),
        inner_(std::move(inner)) {}
  void force(long n, const double *r, const int *z, double *f, double *u,
             double *var, const double *box) override {
    inner_->force(n, r, z, f, u, var, box);
  }
  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return Batched;
  }

private:
  std::shared_ptr<eonc::Potential> inner_;
};
} // namespace

TEST_CASE_METHOD(HessianScratch,
                 "Batched FD Hessian matches the serial one on LJ",
                 "[hessian]") {
  for (const char *scheme : {"one_sided", "central", "fourth"}) {
    CAPTURE(scheme);
    Parameters params;
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    ParametersLoadAccess::hessian_options(params).fd_scheme = scheme;
    auto lj = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(PotType::LJ, params));
    auto serial = std::make_shared<WrappedLJ<false>>(lj, params);
    auto batched = std::make_shared<WrappedLJ<true>>(lj, params);
    Matter serialM(serial, params), batchedM(batched, params);
    serialM.con2matter(std::string("reactant.con"));
    batchedM.con2matter(std::string("reactant.con"));
    VectorXi sub(3);
    sub << 0, 4, 7;
    Hessian hs(params, &serialM), hb(params, &batchedM);
    hs.writeHessianFile(false);
    hb.writeHessianFile(false);
    const MatrixXd a = hs.getHessian(&serialM, sub);
    const MatrixXd b = hb.getHessian(&batchedM, sub);
    REQUIRE(a.rows() == 9);
    REQUIRE(b.rows() == 9);
    for (long i = 0; i < 9; ++i) {
      for (long j = 0; j < 9; ++j) {
        REQUIRE_THAT(b(i, j), Catch::Matchers::WithinAbs(a(i, j), 1e-12));
      }
    }
  }
}

TEST_CASE_METHOD(HessianScratch,
                 "evaluateTogether batches the dirty systems and matches "
                 "single evaluations",
                 "[matter][batch]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto lj = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto batched = std::make_shared<WrappedLJ<true>>(lj, params);
  Matter a(batched, params), b(batched, params), ref(lj, params);
  a.con2matter(std::string("reactant.con"));
  ref.con2matter(std::string("reactant.con"));
  b = a;
  AtomMatrix shifted = b.getPositions();
  shifted(3, 0) += 0.05;
  b.setPositions(shifted);
  REQUIRE(a.needsForceUpdate());
  REQUIRE(b.needsForceUpdate());
  const size_t before = batched->forceCallCounter;
  Matter *const both[] = {&a, &b};
  eonc::evaluateTogether(*batched, both);
  REQUIRE_FALSE(a.needsForceUpdate());
  REQUIRE_FALSE(b.needsForceUpdate());
  // One batch, two systems: the counter advances by the batch size.
  REQUIRE(batched->forceCallCounter == before + 2);
  REQUIRE_THAT(a.getPotentialEnergy(),
               Catch::Matchers::WithinRel(ref.getPotentialEnergy(), 1e-12));
  Matter refB(lj, params);
  refB.con2matter(std::string("reactant.con"));
  refB.setPositions(shifted);
  REQUIRE_THAT(b.getPotentialEnergy(),
               Catch::Matchers::WithinRel(refB.getPotentialEnergy(), 1e-12));
  for (long i = 0; i < b.numberOfAtoms(); ++i) {
    for (int c = 0; c < 3; ++c) {
      REQUIRE_THAT(b.getForces()(i, c),
                   Catch::Matchers::WithinAbs(refB.getForces()(i, c), 1e-10));
    }
  }
  // A clean set costs nothing.
  eonc::evaluateTogether(*batched, both);
  REQUIRE(batched->forceCallCounter == before + 2);
}

namespace {

/// Has config().cutoff and no HessianRangeOf specialization.
struct CutoffOnlyKernel {
  struct Config {
    double cutoff = 4.0;
  };
  explicit CutoffOnlyKernel(const Config &c)
      : cfg(c) {}
  [[nodiscard]] const Config &config() const noexcept { return cfg; }
  [[nodiscard]] rgpot::PotCaps caps() const noexcept { return {}; }
  void forceImpl(const rgpot::ForceInput &, rgpot::ForceOut *) const {}
  void forceBatchImpl(const rgpot::ForceBatch &) const {}
  Config cfg{};
};

/// EAM-like chain: density is a pair sum inside rc, embedding is 1/2 rho^2.
/// A block between atoms in (rc, 2 rc) is carried by the shared neighbour
/// and is zero past 2 rc.
class EmbeddingChain final : public eonc::Potential {
public:
  EmbeddingChain(const Parameters &params, double rc, double declared)
      : eonc::Potential(PotType::LJ, params),
        rc_(rc),
        declared_(declared) {}

  void force(long n, const double *r, const int *, double *f, double *u,
             double *variance, const double *) override {
    if (variance != nullptr) {
      *variance = 0.0;
    }
    std::fill(f, f + 3L * n, 0.0);
    std::vector<double> rho(static_cast<std::size_t>(n), 0.0);
    const auto weight = [&](double rij) {
      if (!(rij < rc_)) {
        return 0.0;
      }
      const double x = 1.0 - rij / rc_;
      return x * x;
    };
    const auto slope = [&](double rij) {
      if (!(rij < rc_)) {
        return 0.0;
      }
      const double x = 1.0 - rij / rc_;
      return -2.0 * x / rc_;
    };
    for (long i = 0; i < n; ++i) {
      for (long j = i + 1; j < n; ++j) {
        const double dx = r[3 * i] - r[3 * j];
        const double dy = r[3 * i + 1] - r[3 * j + 1];
        const double dz = r[3 * i + 2] - r[3 * j + 2];
        const double rij = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double w = weight(rij);
        rho[static_cast<std::size_t>(i)] += w;
        rho[static_cast<std::size_t>(j)] += w;
      }
    }
    double energy = 0.0;
    for (double rhoI : rho) {
      energy += 0.5 * rhoI * rhoI;
    }
    *u = energy;
    for (long i = 0; i < n; ++i) {
      for (long j = i + 1; j < n; ++j) {
        const double dx = r[3 * i] - r[3 * j];
        const double dy = r[3 * i + 1] - r[3 * j + 1];
        const double dz = r[3 * i + 2] - r[3 * j + 2];
        const double rij = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (!(rij > 0.0) || !(rij < rc_)) {
          continue;
        }
        const double coeff = (rho[static_cast<std::size_t>(i)] +
                              rho[static_cast<std::size_t>(j)]) *
                             slope(rij) / rij;
        f[3 * i] -= coeff * dx;
        f[3 * i + 1] -= coeff * dy;
        f[3 * i + 2] -= coeff * dz;
        f[3 * j] += coeff * dx;
        f[3 * j + 1] += coeff * dy;
        f[3 * j + 2] += coeff * dz;
      }
    }
  }

  [[nodiscard]] double hessianRange() const noexcept override {
    return declared_;
  }

private:
  double rc_;
  double declared_;
};

struct ChainHessian {
  MatrixXd H;
  std::size_t calls = 0;
};

struct ChainRequest {
  double declared = 0.0;
  bool fixOne = false;
};

ChainHessian embeddingChainHessian(const Parameters &params,
                                   const ChainRequest &request) {
  constexpr int kN = 8;
  // Planar zig-zag, bond 0.75, 120 degree angle. rc is 1, so atom i and
  // i+2 sit in (rc, 2 rc) and share a neighbour. The y offset makes the
  // yy and xy blocks nonzero.
  constexpr double kBond = 0.75;
  const double step = kBond * std::sqrt(3.0) / 2.0;
  const double height = kBond / 2.0;
  auto pot = std::make_shared<EmbeddingChain>(params, 1.0, request.declared);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->resize(kN);
  AtomMatrix pos(kN, 3);
  pos.setZero();
  for (int i = 0; i < kN; ++i) {
    matter->setAtomicNr(i, 1);
    matter->setMass(i, 1.0);
    pos(i, 0) = step * static_cast<double>(i);
    pos(i, 1) = (i % 2 == 0) ? 0.0 : height;
  }
  matter->setPositions(pos);
  matter->setCell(Matrix3d::Identity() * 40.0);
  matter->setPeriodic(false);
  if (request.fixOne) {
    matter->setFixed(0, 1);
  }
  VectorXi all(kN);
  for (int i = 0; i < kN; ++i) {
    all(i) = i;
  }
  pot->forceCallCounter.store(0);
  Hessian hess(params, matter.get());
  hess.writeHessianFile(false);
  ChainHessian out;
  out.H = hess.getHessian(matter.get(), all);
  out.calls = pot->forceCallCounter.load();
  return out;
}

/// Forwards forces to a real kernel and reports a chosen Hessian range.
/// The stored cell is forwarded so EDIP still sees a box wide enough for
/// its cell-width check while the image itself is non-periodic.
class DeclaredRange final : public eonc::Potential {
public:
  DeclaredRange(std::shared_ptr<eonc::Potential> inner,
                const Parameters &params, double range)
      : eonc::Potential(inner->getType(), params),
        inner_(std::move(inner)),
        range_(range) {}

  void force(long n, const double *r, const int *atomicNrs, double *f,
             double *u, double *variance, const double *box) override {
    inner_->force(n, r, atomicNrs, f, u, variance, box);
  }

  [[nodiscard]] double hessianRange() const noexcept override { return range_; }

  [[nodiscard]] bool forwardsStoredCell() const noexcept override {
    return true;
  }

  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return false;
  }

private:
  std::shared_ptr<eonc::Potential> inner_;
  double range_;
};

struct KernelProbe {
  const char *name;
  PotType type;
  int z;
  int zGuest;
  double bond;
  /// Perpendicular rise of the zig-zag. 0 selects a 120 degree angle.
  double rise;
  double interactionCutoff;
  bool reachesPastCutoff;
};

double zigRise(const KernelProbe &probe) {
  return probe.rise > 0.0 ? probe.rise : probe.bond / 2.0;
}

double zigStep(const KernelProbe &probe) {
  const double rise = zigRise(probe);
  return std::sqrt(probe.bond * probe.bond - rise * rise);
}

struct KernelReport {
  int atoms = 0;
  int colors = 0;
  std::size_t serialCalls = 0;
  std::size_t coloredCalls = 0;
  double maxAbs = 0.0;
  double largestNonzero = 0.0;
  double coloredErr = 0.0;
  double shortErr = 0.0;
};

double blockAbs(const MatrixXd &H, int atomI, int atomJ) {
  double m = 0.0;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      m = std::max(m, std::abs(H(3 * atomI + a, 3 * atomJ + b)));
    }
  }
  return m;
}

ChainHessian kernelHessian(const Parameters &params,
                           const std::shared_ptr<eonc::Potential> &pot, int n,
                           double step, double height, int z, int zGuest) {
  auto matter = std::make_shared<Matter>(pot, params);
  matter->resize(n);
  AtomMatrix pos(n, 3);
  pos.setZero();
  for (int i = 0; i < n; ++i) {
    const bool guest = zGuest != 0 && (i == 1 || i == n / 2);
    matter->setAtomicNr(i, guest ? zGuest : z);
    matter->setMass(i, 1.0);
    pos(i, 0) = step * static_cast<double>(i);
    pos(i, 1) = (i % 2 == 0) ? 0.0 : height;
  }
  matter->setPositions(pos);
  // Wide enough for EDIP's cell-width check. Diagonal, so the box the
  // kernel sees is the same matrix Fortran stores as columns.
  matter->setCell(Matrix3d::Identity() * 200.0);
  matter->setPeriodic(false);
  VectorXi all(n);
  for (int i = 0; i < n; ++i) {
    all(i) = i;
  }
  pot->forceCallCounter.store(0);
  Hessian hess(params, matter.get());
  hess.writeHessianFile(false);
  ChainHessian out;
  out.H = hess.getHessian(matter.get(), all);
  out.calls = pot->forceCallCounter.load();
  return out;
}

/// Spring cluster inside a range that covers every pair. Batching is on
/// so a fallback off the coloured path is visible: the coloured loop
/// never calls forceBatch.
class BatchCluster final : public eonc::Potential {
public:
  explicit BatchCluster(const Parameters &params)
      : eonc::Potential(PotType::LJ, params) {}

  void force(long n, const double *r, const int *, double *f, double *u,
             double *variance, const double *) override {
    if (variance != nullptr) {
      *variance = 0.0;
    }
    std::fill(f, f + 3L * n, 0.0);
    double energy = 0.0;
    for (long i = 0; i < n; ++i) {
      for (long j = i + 1; j < n; ++j) {
        const double dx = r[3 * i] - r[3 * j];
        const double dy = r[3 * i + 1] - r[3 * j + 1];
        const double dz = r[3 * i + 2] - r[3 * j + 2];
        const double rij = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (!(rij > 0.0) || !(rij < 5.0)) {
          continue;
        }
        const double stretch = rij - 1.2;
        energy += 0.5 * stretch * stretch;
        const double coeff = stretch / rij;
        f[3 * i] -= coeff * dx;
        f[3 * i + 1] -= coeff * dy;
        f[3 * i + 2] -= coeff * dz;
        f[3 * j] += coeff * dx;
        f[3 * j + 1] += coeff * dy;
        f[3 * j + 2] += coeff * dz;
      }
    }
    *u = energy;
  }

  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override {
    ++batchCalls;
    eonc::Potential::forceBatch(nSystems, nAtoms, positions, atomicNrs, forces,
                                energies, variances, boxes);
  }

  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return true;
  }

  [[nodiscard]] double hessianRange() const noexcept override { return 30.0; }

  std::size_t batchCalls = 0;
};

template <class RPot>
void requireManyBodyRange(PotType ptype, const Parameters &params,
                          bool enabled) {
  RgpotAdapter<RPot> pot(ptype, params);
  using Trait = eonc::HessianRangeOf<RPot>;
  if (enabled) {
    REQUIRE(Trait::locality == eonc::HessianLocality::AtomCentred);
    REQUIRE(Trait::interactionCutoff > 0.0);
    const double want = 2.0 * Trait::interactionCutoff;
    REQUIRE(want > Trait::interactionCutoff);
    REQUIRE(pot.hessianRange() == want);
  } else {
    REQUIRE(pot.hessianRange() == 0.0);
  }
}

} // namespace

TEST_CASE("Pair kernels opt in at the interaction cutoff", "[hessian][color]") {
  Parameters params;
  auto lj = eonc::helpers::makePotential(PotType::LJ, params);
  REQUIRE(lj->hessianRange() == 15.0);
  auto cluster = eonc::helpers::makePotential(PotType::LJCLUSTER, params);
  REQUIRE(cluster->hessianRange() == 15.0);
  auto morse = eonc::helpers::makePotential(PotType::MORSE_PT, params);
  REQUIRE(morse->hessianRange() == 9.5);
  ParametersLoadAccess::zbl_options(params).cut_inner = 1.5;
  ParametersLoadAccess::zbl_options(params).cut_global = 2.25;
  auto zbl = eonc::helpers::makePotential(PotType::ZBL, params);
  REQUIRE(zbl->hessianRange() == 2.25);
}

TEST_CASE("Many-body rgpot kernels report a disabled or doubled range",
          "[hessian][color]") {
  Parameters params;
  requireManyBodyRange<rgpot::fortranpots::EAMAlPot>(PotType::EAM_AL, params,
                                                     true);
  requireManyBodyRange<rgpot::fortranpots::SWPot>(PotType::SW_SI, params, true);
  requireManyBodyRange<rgpot::fortranpots::TersoffPot>(PotType::TERSOFF_SI,
                                                       params, true);
  requireManyBodyRange<rgpot::fortranpots::EDIPPot>(PotType::EDIP, params,
                                                    true);
  requireManyBodyRange<rgpot::fortranpots::LenoskyPot>(PotType::LENOSKY_SI,
                                                       params, true);
  requireManyBodyRange<rgpot::fortranpots::FeHePot>(PotType::FEHE, params,
                                                    true);
  requireManyBodyRange<rgpot::fortranpots::CuH2Pot>(PotType::CUH2, params,
                                                    true);
  requireManyBodyRange<rgpot::fortranpots::WaterHPot>(PotType::TIP4P_H, params,
                                                      false);
  // Factory arm agrees with the trait for the kernels jobs actually build.
  auto eam = eonc::helpers::makePotential(PotType::EAM_AL, params);
  REQUIRE(eam->hessianRange() ==
          2.0 * eonc::HessianRangeOf<
                    rgpot::fortranpots::EAMAlPot>::interactionCutoff);
  auto sw = eonc::helpers::makePotential(PotType::SW_SI, params);
  REQUIRE(sw->hessianRange() ==
          eonc::HessianRangeOf<rgpot::fortranpots::SWPot>::range(
              rgpot::fortranpots::SWPot{}));
}

TEST_CASE("A cutoff field without an opt-in is not colored",
          "[hessian][color]") {
  static_assert(eonc::HessianRangeOf<CutoffOnlyKernel>::locality ==
                eonc::HessianLocality::Unknown);
  static_assert(eonc::HessianRangeOf<CutoffOnlyKernel>::interactionCutoff ==
                0.0);
#if EON_RGPOT_HAS_UMA
  static_assert(eonc::HessianRangeOf<rgpot::UmaPot>::locality ==
                eonc::HessianLocality::MessagePassing);
  static_assert(eonc::HessianRangeOf<rgpot::UmaPot>::interactionCutoff == 0.0);
#endif
#if EON_RGPOT_HAS_WATER_POTS
  static_assert(eonc::HessianRangeOf<rgpot::TIP4PPot>::locality ==
                eonc::HessianLocality::MolecularSite);
  static_assert(eonc::HessianRangeOf<rgpot::TIP4PPot>::interactionCutoff ==
                0.0);
  static_assert(eonc::HessianRangeOf<rgpot::SPCEPot>::locality ==
                eonc::HessianLocality::MolecularSite);
  static_assert(eonc::HessianRangeOf<rgpot::TIP4PPtPot>::locality ==
                eonc::HessianLocality::MolecularSite);
#endif
  Parameters params;
  RgpotAdapter<CutoffOnlyKernel> pot(PotType::LJ, params,
                                     CutoffOnlyKernel::Config{.cutoff = 4.0});
  REQUIRE(pot.kernel().config().cutoff == 4.0);
  REQUIRE(pot.hessianRange() == 0.0);
}

void checkEmbeddingChain(const Parameters &params, std::size_t serialCalls,
                         const char *label) {
  const ChainHessian serial = embeddingChainHessian(params, ChainRequest{0.0});
  const ChainHessian colored = embeddingChainHessian(params, ChainRequest{2.0});
  const ChainHessian shortRange =
      embeddingChainHessian(params, ChainRequest{1.0});
  REQUIRE(serial.H.rows() == 24);
  REQUIRE(colored.H.rows() == 24);
  REQUIRE(shortRange.H.rows() == 24);
  REQUIRE(serial.calls == serialCalls);
  REQUIRE(colored.calls < serial.calls);
  REQUIRE(shortRange.calls < serial.calls);

  double maxAbs = 0.0;
  double errColored = 0.0;
  double errShort = 0.0;
  double offAxis = 0.0;
  for (int i = 0; i < 24; ++i) {
    for (int j = 0; j < 24; ++j) {
      maxAbs = std::max(maxAbs, std::abs(serial.H(i, j)));
      errColored =
          std::max(errColored, std::abs(colored.H(i, j) - serial.H(i, j)));
      errShort =
          std::max(errShort, std::abs(shortRange.H(i, j) - serial.H(i, j)));
      if (i % 3 != 0 || j % 3 != 0) {
        offAxis = std::max(offAxis, std::abs(serial.H(i, j)));
      }
    }
  }
  std::fprintf(stderr,
               "EMBED %s serial_calls=%zu colored_calls=%zu short_calls=%zu "
               "maxAbs=%.6e colored_err=%.6e short_err=%.6e off_axis=%.6e\n",
               label, serial.calls, colored.calls, shortRange.calls, maxAbs,
               errColored, errShort, offAxis);
  REQUIRE(maxAbs > 0.0);
  REQUIRE(offAxis > 1e-6);
  REQUIRE(errColored < 1e-8);
  REQUIRE(errShort > 1e-2 * maxAbs);
}

TEST_CASE("Embedding-chain coloring matches serial only at twice the cutoff",
          "[hessian][color]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::main_options(params).finiteDifference = 1e-5;
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  ParametersLoadAccess::hessian_options(params).fd_scheme = "central";
  // Eight atoms, 24 coordinates, central stencil.
  checkEmbeddingChain(params, static_cast<std::size_t>(1 + 2 * 24), "central");

  ParametersLoadAccess::hessian_options(params).fd_scheme = "fourth";
  checkEmbeddingChain(params, static_cast<std::size_t>(1 + 4 * 24), "fourth");

  // removeNetForce would couple every atom, except a fixed atom lets the
  // coloured path run.
  ParametersLoadAccess::hessian_options(params).fd_scheme = "central";
  ParametersLoadAccess::main_options(params).removeNetForce = true;
  const ChainHessian serial =
      embeddingChainHessian(params, ChainRequest{0.0, true});
  const ChainHessian colored =
      embeddingChainHessian(params, ChainRequest{2.0, true});
  const ChainHessian shortRange =
      embeddingChainHessian(params, ChainRequest{1.0, true});
  REQUIRE(serial.calls == static_cast<std::size_t>(1 + 2 * 24));
  REQUIRE(colored.calls < serial.calls);
  double maxAbs = 0.0;
  double errColored = 0.0;
  double errShort = 0.0;
  for (int i = 0; i < 24; ++i) {
    for (int j = 0; j < 24; ++j) {
      maxAbs = std::max(maxAbs, std::abs(serial.H(i, j)));
      errColored =
          std::max(errColored, std::abs(colored.H(i, j) - serial.H(i, j)));
      errShort =
          std::max(errShort, std::abs(shortRange.H(i, j) - serial.H(i, j)));
    }
  }
  std::fprintf(stderr,
               "EMBED fixed serial_calls=%zu colored_calls=%zu short_calls=%zu "
               "maxAbs=%.6e colored_err=%.6e short_err=%.6e\n",
               serial.calls, colored.calls, shortRange.calls, maxAbs,
               errColored, errShort);
  REQUIRE(errColored < 1e-8);
  REQUIRE(errShort > 1e-2 * maxAbs);
}

TEST_CASE("Opted-in kernels match the serial matrix at the declared range",
          "[hessian][color]") {
  const KernelProbe probes[] = {
      {"EAM_AL", PotType::EAM_AL, 13, 0, 2.86, 0.0, 7.1, true},
      {"SW_SI", PotType::SW_SI, 14, 0, 2.35, 0.0, 1.8 * 2.0951, true},
      {"TERSOFF_SI", PotType::TERSOFF_SI, 14, 0, 2.35, 0.0, 3.0, true},
      {"EDIP", PotType::EDIP, 14, 0, 2.35, 0.0, 3.1213820, true},
      // A 120 degree chain puts atom i+2 inside the 4.5 pair cutoff.
      // A flatter zig-zag keeps the 2.35 bond and puts i+2 in (4.5, 9).
      {"LENOSKY_SI", PotType::LENOSKY_SI, 14, 0, 2.35, 0.45, 4.5, true},
      // Fe-Fe density cuts at 4.2. A taller zig-zag keeps the 2.48 bond
      // and puts a shared neighbour of an (rc, 2 rc) pair inside 4.2.
      {"FEHE", PotType::FEHE, 26, 2, 2.48, 1.40, 5.4, true},
      {"CUH2", PotType::CUH2, 29, 1, 2.55, 0.0, 6.1, true},
      {"ZBL", PotType::ZBL, 14, 0, 1.50, 0.0, 2.5, false},
  };

  for (const KernelProbe &probe : probes) {
    INFO(probe.name);
    Parameters params;
    ParametersLoadAccess::potential_options(params).potential = probe.type;
    ParametersLoadAccess::main_options(params).finiteDifference = 1e-5;
    ParametersLoadAccess::main_options(params).removeNetForce = false;
    ParametersLoadAccess::hessian_options(params).fd_scheme = "central";
    if (probe.type == PotType::ZBL) {
      ParametersLoadAccess::zbl_options(params).cut_inner = 2.0;
      ParametersLoadAccess::zbl_options(params).cut_global = 2.5;
    }
    auto real = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(probe.type, params));
    const double declared = real->hessianRange();
    if (probe.reachesPastCutoff) {
      REQUIRE(declared == 2.0 * probe.interactionCutoff);
    } else {
      REQUIRE(declared == probe.interactionCutoff);
    }

    const double step = zigStep(probe);
    const double height = zigRise(probe);
    const int kMax = static_cast<int>(std::floor(declared / step));
    const int n = std::max(8, 2 * kMax + 4);
    auto serialPot = std::make_shared<DeclaredRange>(real, params, 0.0);
    auto coloredPot = std::make_shared<DeclaredRange>(real, params, declared);
    const ChainHessian serial = kernelHessian(params, serialPot, n, step,
                                              height, probe.z, probe.zGuest);
    const ChainHessian colored = kernelHessian(params, coloredPot, n, step,
                                               height, probe.z, probe.zGuest);
    REQUIRE(serial.H.rows() == 3 * n);
    REQUIRE(colored.H.rows() == 3 * n);
    REQUIRE(serial.calls == static_cast<std::size_t>(1 + 2 * 3 * n));
    REQUIRE(colored.calls < serial.calls);

    auto geom = std::make_shared<Matter>(serialPot, params);
    geom->resize(n);
    AtomMatrix pos(n, 3);
    pos.setZero();
    for (int i = 0; i < n; ++i) {
      pos(i, 0) = step * static_cast<double>(i);
      pos(i, 1) = (i % 2 == 0) ? 0.0 : height;
    }
    VectorXi all(n);
    for (int i = 0; i < n; ++i) {
      all(i) = i;
    }
    geom->setPositions(pos);
    geom->setCell(Matrix3d::Identity() * 200.0);
    geom->setPeriodic(false);
    const std::vector<int> colors =
        eonc::colorMobileCutoffGraph(*geom, all, declared);
    int nColors = 0;
    for (int c : colors) {
      nColors = std::max(nColors, c + 1);
    }
    REQUIRE(nColors >= 2);
    REQUIRE(nColors < n);

    double maxAbs = 0.0;
    double coloredErr = 0.0;
    for (int i = 0; i < 3 * n; ++i) {
      for (int j = 0; j < 3 * n; ++j) {
        maxAbs = std::max(maxAbs, std::abs(serial.H(i, j)));
        coloredErr =
            std::max(coloredErr, std::abs(colored.H(i, j) - serial.H(i, j)));
      }
    }
    REQUIRE(maxAbs > 0.0);
    REQUIRE(coloredErr <= 1e-8 * maxAbs);

    double largestNonzero = 0.0;
    double pastCutoff = 0.0;
    double pastRange = 0.0;
    bool sawPastCutoff = false;
    bool sawPastRange = false;
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        const double dx = pos(i, 0) - pos(j, 0);
        const double dy = pos(i, 1) - pos(j, 1);
        const double dist = std::sqrt(dx * dx + dy * dy);
        const double block = blockAbs(serial.H, i, j);
        if (block > 0.0) {
          largestNonzero = std::max(largestNonzero, dist);
        }
        if (dist > probe.interactionCutoff) {
          sawPastCutoff = true;
          pastCutoff = std::max(pastCutoff, block);
        }
        if (dist > declared) {
          sawPastRange = true;
          pastRange = std::max(pastRange, block);
        }
      }
    }
    REQUIRE(sawPastRange);
    REQUIRE(pastRange == 0.0);
    if (probe.reachesPastCutoff) {
      REQUIRE(sawPastCutoff);
      REQUIRE(pastCutoff > 1e-6 * maxAbs);
    }
    REQUIRE(largestNonzero > 0.0);
    REQUIRE(largestNonzero <= declared);

    double shortErr = 0.0;
    if (probe.name == std::string("SW_SI") ||
        probe.name == std::string("EAM_AL")) {
      auto shortPot = std::make_shared<DeclaredRange>(real, params,
                                                      probe.interactionCutoff);
      const ChainHessian shortH = kernelHessian(params, shortPot, n, step,
                                                height, probe.z, probe.zGuest);
      REQUIRE(shortH.calls < serial.calls);
      for (int i = 0; i < 3 * n; ++i) {
        for (int j = 0; j < 3 * n; ++j) {
          shortErr =
              std::max(shortErr, std::abs(shortH.H(i, j) - serial.H(i, j)));
        }
      }
      REQUIRE(shortErr > 1e-2 * maxAbs);
    }

    std::fprintf(stderr,
                 "KERNEL %s atoms=%d colors=%d serial_calls=%zu "
                 "colored_calls=%zu maxAbs=%.6e largest_nonzero=%.6f "
                 "colored_err=%.6e short_err=%.6e\n",
                 probe.name, n, nColors, serial.calls, colored.calls, maxAbs,
                 largestNonzero, coloredErr, shortErr);
  }
}

TEST_CASE("Coloring falls back when every mobile atom needs its own color",
          "[hessian][color]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::main_options(params).finiteDifference = 1e-5;
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  ParametersLoadAccess::hessian_options(params).fd_scheme = "central";

  constexpr int kN = 4;
  auto place = [kN](const std::shared_ptr<Matter> &matter) {
    matter->resize(kN);
    for (long i = 0; i < kN; ++i) {
      matter->setAtomicNr(i, 1);
      matter->setMass(i, 1.0);
    }
    AtomMatrix pos(kN, 3);
    pos.setZero();
    pos(1, 0) = 1.0;
    pos(2, 0) = 0.4;
    pos(2, 1) = 0.7;
    pos(3, 0) = 0.5;
    pos(3, 2) = 0.6;
    matter->setPositions(pos);
    matter->setCell(Matrix3d::Identity() * 40.0);
    matter->setPeriodic(false);
  };

  auto lj = std::make_shared<RgpotAdapter<rgpot::LJPot>>(
      PotType::LJ, params, rgpot::LJConfig{.cutoff = 20.0});
  REQUIRE(lj->hessianRange() == 20.0);
  auto ljMatter = std::make_shared<Matter>(lj, params);
  place(ljMatter);
  VectorXi all(kN);
  all << 0, 1, 2, 3;
  const std::vector<int> colors =
      eonc::colorMobileCutoffGraph(*ljMatter, all, lj->hessianRange());
  int nColors = 0;
  for (int c : colors) {
    nColors = std::max(nColors, c + 1);
  }
  REQUIRE(nColors == kN);

  lj->forceCallCounter.store(0);
  Hessian ljHess(params, ljMatter.get());
  ljHess.writeHessianFile(false);
  const MatrixXd ljH = ljHess.getHessian(ljMatter.get(), all);
  REQUIRE(ljH.rows() == 3 * kN);
  // Diameter is below the range, so colouring saves nothing and the
  // per-column loop runs: one base gradient and two central samples
  // on each of 12 coordinates.
  REQUIRE(lj->forceCallCounter.load() == static_cast<std::size_t>(1 + 2 * 12));

  auto batch = std::make_shared<BatchCluster>(params);
  auto batchMatter = std::make_shared<Matter>(batch, params);
  place(batchMatter);
  const std::vector<int> batchColors =
      eonc::colorMobileCutoffGraph(*batchMatter, all, batch->hessianRange());
  int batchNColors = 0;
  for (int c : batchColors) {
    batchNColors = std::max(batchNColors, c + 1);
  }
  REQUIRE(batchNColors == kN);
  batch->forceCallCounter.store(0);
  batch->batchCalls = 0;
  Hessian batchHess(params, batchMatter.get());
  batchHess.writeHessianFile(false);
  const MatrixXd batchH = batchHess.getHessian(batchMatter.get(), all);
  REQUIRE(batchH.rows() == 3 * kN);
  REQUIRE(batch->batchCalls > 0);
  REQUIRE(batch->forceCallCounter.load() ==
          static_cast<std::size_t>(1 + 2 * 12));
}

} /* namespace tests */

namespace tests {
TEST_CASE_METHOD(
    HessianScratch,
    "Atomic checkpoint replacement preserves valid data on failure",
    "[hessian][checkpoint]") {
  const std::string path = "atomic.ckpt";
  const std::string original = "eon_hess_ckpt 1 1\n1.2345678901234567\n";
  {
    std::ofstream out(path);
    out << original;
  }
  auto contents = [&]() {
    std::ifstream in(path);
    return std::string(std::istreambuf_iterator<char>(in), {});
  };
  SECTION("serializer failure") {
    REQUIRE_FALSE(eonc::helpers::atomicWriteFile(path, [](FILE *out) {
      std::fputs("incomplete", out);
      return false;
    }));
    REQUIRE(contents() == original);
  }
  SECTION("interrupted serializer") {
    REQUIRE_THROWS_AS(eonc::helpers::atomicWriteFile(
                          path,
                          [](FILE *out) -> bool {
                            std::fputs("incomplete", out);
                            throw std::runtime_error("interrupted write");
                          }),
                      std::runtime_error);
    REQUIRE(contents() == original);
  }
  SECTION("invalid destination") {
    REQUIRE_FALSE(eonc::helpers::atomicWriteFile("missing/checkpoint",
                                                 [](FILE *) { return true; }));
    REQUIRE(contents() == original);
  }
  SECTION("replacement failure") {
    std::filesystem::create_directory("destination");
    REQUIRE_FALSE(eonc::helpers::atomicWriteFile("destination", [](FILE *out) {
      return std::fputs("complete", out) >= 0;
    }));
    REQUIRE(std::filesystem::is_directory("destination"));
    REQUIRE(contents() == original);
  }
  SECTION("successful replacement") {
    const std::string replacement = "eon_hess_ckpt 1 1\n9.8765432109876543\n";
    REQUIRE(eonc::helpers::atomicWriteFile(path, [&](FILE *out) {
      return std::fputs(replacement.c_str(), out) >= 0;
    }));
    REQUIRE(contents() == replacement);
  }
  for (const auto &entry : std::filesystem::directory_iterator(".")) {
    REQUIRE(entry.path().filename().string().find(".tmp.") ==
            std::string::npos);
  }
}

TEST_CASE_METHOD(HessianScratch, "Final Hessian retains every binary64 entry",
                 "[hessian][checkpoint]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  Matter matter(pot, params);
  matter.con2matter(std::string("reactant.con"));
  VectorXi atoms(2);
  atoms << 0, 1;
  Hessian hessian(params, &matter);
  const MatrixXd matrix = hessian.getHessian(&matter, atoms);
  REQUIRE(matrix.rows() == 6);
  std::ifstream input("hessian.dat");
  REQUIRE(input.good());
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      double value = 0;
      REQUIRE(static_cast<bool>(input >> value));
      REQUIRE(value == matrix(i, j));
    }
  }
}
} // namespace tests

namespace tests {
namespace {
class InterruptedLJ final : public eonc::Potential {
public:
  InterruptedLJ(std::shared_ptr<eonc::Potential> inner, const Parameters &p)
      : eonc::Potential(PotType::LJ, p),
        inner_(std::move(inner)) {}
  int calls = 0;
  bool interrupt = true;
  void force(long n, const double *r, const int *z, double *f, double *u,
             double *var, const double *box) override {
    if (++calls > 3 && interrupt) {
      throw std::runtime_error("interrupted finite-difference evaluation");
    }
    inner_->force(n, r, z, f, u, var, box);
  }

private:
  std::shared_ptr<eonc::Potential> inner_;
};
} // namespace

TEST_CASE_METHOD(HessianScratch,
                 "Hessian resume rejects changed finite-difference context",
                 "[hessian][checkpoint-context]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  ParametersLoadAccess::hessian_options(params).fd_scheme = "one_sided";
  ParametersLoadAccess::hessian_options(params).checkpoint_path =
      "context.ckpt";
  auto lj = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(PotType::LJ, params));
  auto potential = std::make_shared<InterruptedLJ>(lj, params);
  Matter matter(potential, params);
  matter.con2matter(std::string("reactant.con"));
  matter.setPositions(matter.getPositions());
  VectorXi atoms(2);
  atoms << 0, 1;
  Hessian interrupted(params, &matter);
  REQUIRE_THROWS_AS(interrupted.getHessian(&matter, atoms), std::runtime_error);
  REQUIRE(std::filesystem::exists("context.ckpt"));
  {
    std::ifstream checkpoint("context.ckpt");
    std::string magic;
    int size = 0, nextColumn = 0;
    REQUIRE(static_cast<bool>(checkpoint >> magic >> size >> nextColumn));
    REQUIRE(magic == "eon_hess_ckpt");
    REQUIRE(size == 6);
    REQUIRE(nextColumn == 2);
  }
  potential->interrupt = false;
  potential->calls = 0;
  int expectedResumeCalls = 7;
  SECTION("unchanged context") { expectedResumeCalls = 5; }
  SECTION("coordinates") {
    AtomMatrix positions = matter.getPositions();
    positions(0, 0) += 0.2;
    matter.setPositions(positions);
  }
  SECTION("displacement") {
    ParametersLoadAccess::main_options(params).finiteDifference *= 2;
  }
  SECTION("stencil") {
    ParametersLoadAccess::hessian_options(params).fd_scheme = "central";
    expectedResumeCalls = 13;
  }
  SECTION("mobile atom order") { atoms << 1, 0; }
  ParametersLoadAccess::hessian_options(params).resume = true;
  Hessian resumed(params, &matter);
  const MatrixXd actual = resumed.getHessian(&matter, atoms);
  REQUIRE(potential->calls == expectedResumeCalls);
  ParametersLoadAccess::hessian_options(params).resume = false;
  Hessian fresh(params, &matter);
  const MatrixXd expected = fresh.getHessian(&matter, atoms);
  REQUIRE(actual.rows() == 6);
  REQUIRE(expected.rows() == 6);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      REQUIRE_THAT(actual(i, j),
                   Catch::Matchers::WithinAbs(expected(i, j), 1e-12));
    }
  }
}
TEST_CASE("Linkcell cutoff rows preserve periodic atom-image semantics",
          "[Hessian][Linkcell][MinimumImageLibrary]") {
  const bool skew = GENERATE(false, true);
  const bool full = GENERATE(false, true);
  Matrix3d cell = Matrix3d::Identity();
  if (skew)
    cell << 1, 0, 0, .35, 1, 0, .1, .2, 1.1;
  AtomMatrix positions(3, 3);
  positions << .02, .03, .04, .93, .04, .05, .4, .7, .2;
  constexpr double cutoff = 1.35;
  using Key = std::array<int, 5>;
  auto canonical = [&](Key key) {
    if (full)
      return key;
    const Key reverse{key[1], key[0], -key[2], -key[3], -key[4]};
    return std::min(key, reverse);
  };
  std::map<Key, double> expected, actual, reference;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int a = -4; a <= 4; ++a)
        for (int b = -4; b <= 4; ++b)
          for (int c = -4; c <= 4; ++c) {
            if (i == j && a == 0 && b == 0 && c == 0)
              continue;
            const Eigen::RowVector3d d = positions.row(j) - positions.row(i) +
                                         a * cell.row(0) + b * cell.row(1) +
                                         c * cell.row(2);
            if (d.squaredNorm() < cutoff * cutoff)
              expected[canonical({i, j, a, b, c})] = d.squaredNorm();
          }
  const auto rows = eonc::periodicCutoffNeighbors(positions.data(), 3,
                                                  cell.data(), cutoff, full);
  for (const auto &pair : rows) {
    const Key key = canonical(
        {pair.first, pair.second, pair.shift[0], pair.shift[1], pair.shift[2]});
    REQUIRE(actual.emplace(key, pair.distanceSquared).second);
  }
  eonc::VesinNeighbors vesin;
  eonc::VesinNeighbors::Options options;
  options.cutoff = cutoff;
  options.full = full;
  options.return_shifts = true;
  options.return_distances = true;
  vesin.compute(positions.data(), 3, cell.data(), options);
  for (size_t p = 0; p < vesin.size(); ++p) {
    const auto shift = vesin.shift(p);
    const Key key =
        canonical({static_cast<int>(vesin.i(p)), static_cast<int>(vesin.j(p)),
                   shift[0], shift[1], shift[2]});
    REQUIRE(reference.emplace(key, std::pow(vesin.distance(p), 2)).second);
  }
  REQUIRE(actual.size() == expected.size());
  REQUIRE(reference.size() == expected.size());
  for (const auto &[key, d2] : expected) {
    REQUIRE(actual.count(key) == 1);
    REQUIRE(reference.count(key) == 1);
    REQUIRE_THAT(actual.at(key), Catch::Matchers::WithinAbs(d2, 1e-12));
    REQUIRE_THAT(reference.at(key), Catch::Matchers::WithinAbs(d2, 1e-12));
  }
}

} // namespace tests
