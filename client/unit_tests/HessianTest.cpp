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
#include "eon/ConFileIO.h"
#include "eon/Davidson.h"
#include "eon/FiniteDifference.h"
#include "eon/Lanczos.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/SafeMath.h"
#include "eon/Tunneling.h"
#include "eon/potentials/RgpotAdapter/RgpotAdapter.h"
#include "rgpot/LennardJones/LJPot.hpp"
#include <nlohmann/json.hpp>

#include <algorithm>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
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

TEST_CASE_METHOD(HessianScratch, "Hessian getFreqs returns finite eigenvalues",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
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

TEST_CASE_METHOD(HessianScratch,
                 "Hessian getFreqs rejects out-of-range atom indices",
                 "[hessian]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
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
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
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

  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
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

TEST_CASE_METHOD(HessianScratch,
                 "Hessian on Pt frozen layers system handles mixed fixed/free",
                 "[hessian][morse_pt]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::MORSE_PT;
  auto pot = eonc::helpers::makePotential(PotType::MORSE_PT, params);
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
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
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
  REQUIRE(pot->finiteCutoff() == kCutoff);

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
      eonc::colorMobileCutoffGraph(*matter, all, pot->finiteCutoff());
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
      eonc::colorMobileCutoffGraph(*matter, all, pot->finiteCutoff());
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
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
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
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
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

} /* namespace tests */
