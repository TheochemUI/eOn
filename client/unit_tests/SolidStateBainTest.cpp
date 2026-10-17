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
// Bain path of a Lennard-Jones crystal on the solid-state band.
//
// LJ with u0 1 eV and psi 1 Angstrom, unshifted, times the C^2 quintic
// switch over [2.0, 2.2] Angstrom (rgpot LJConfig cutoff 2.2,
// switch_width 0.2). bcc and fcc are the two-atom body-centred
// tetragonal cell diag(a, a, c). Both atoms sit on inversion centres, so
// the transformation is the cell strain alone, and the only lattice wave
// vector of this cell is k = 0: the tetragonal saddle of the scan is the
// saddle of this cell. A supercell admits the k = (1/2, 0, 1/2) shuffle,
// which is unstable at that saddle, and its band leaves the tetragonal
// path. The cell is narrower than twice the cutoff, so the potential
// here sums every lattice translation inside the cutoff instead of
// taking the minimum image.
//
// Reference: validation/bain_lj_reference.py in rgsaddle, a relaxed scan
// of the energy per atom over c/a with spacing 2.5e-4:
//   bcc     a = 1.228441371,                  E = -7.4148256256 eV/atom
//   saddle  c/a = 1.09331985,                 E = -7.3557899883 eV/atom
//   fcc     a = 1.100691058, c = 1.556612222, E = -7.6698410342 eV/atom
// The climbing image must reach the scan's c/a within the spacing, and
// its barrier within 1.2e-6 eV/atom of 0.0590356373, the energy change
// over one scan step at the saddle. rgsaddle's tests/c/solid_bain_lj.c
// runs the same band.
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/NudgedElasticBand.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>

namespace tests {
namespace {

static eonc::helpers::test::QuillTestLogger quillSetup;

constexpr double kCutoff = 2.2;
constexpr double kWidth = 0.2;
constexpr double kBccA = 1.228441371;
constexpr double kFccA = 1.100691058;
constexpr double kFccC = 1.556612222;
constexpr double kBccEnergy = -7.4148256256;
constexpr double kSaddleQ = 1.09331985;
constexpr double kBarrier = 0.0590356373;
constexpr double kScanSpacing = 2.5e-4;
constexpr double kBarrierTol = 1.2e-6;

// Switched LJ over every lattice translation inside the cutoff, with the
// Cauchy stress sigma = (1/V) dE/deps that the solid-state band reads.
struct LatticeLJ final : eonc::Potential {
  Matrix3d last{Matrix3d::Zero()};

  LatticeLJ()
      : Potential(eonc::PotType::LJ) {}

  [[nodiscard]] bool computesStress() const noexcept override { return true; }

  [[nodiscard]] Matrix3d cauchyStress() const override { return last; }

  static void pair(double r, double &phi, double &dphi) {
    const double sr6 = std::pow(1.0 / r, 6);
    const double v = 4.0 * (sr6 * sr6 - sr6);
    const double dv = 4.0 * (-12.0 * sr6 * sr6 + 6.0 * sr6) / r;
    const double t = std::clamp((r - (kCutoff - kWidth)) / kWidth, 0.0, 1.0);
    const double s = 1.0 + t * t * t * (-10.0 + t * (15.0 - 6.0 * t));
    const double ds = -30.0 * t * t * (1.0 - t) * (1.0 - t) / kWidth;
    phi = v * s;
    dphi = dv * s + v * ds;
  }

  void force(long nAtoms, const double *positions, const int * /*atomicNrs*/,
             double *forces, double *energy, double *variance,
             const double *box) override {
    // Row-major box; its rows are the lattice vectors.
    const Matrix3d h = Eigen::Map<const Matrix3d>(box);
    const Matrix3d inv = h.inverse();
    int reach[3];
    for (int k = 0; k < 3; ++k) {
      reach[k] =
          static_cast<int>(std::ceil(kCutoff * inv.col(k).norm() + 0.5));
      if (reach[k] > 64) {
        throw std::runtime_error("LatticeLJ: the cell has collapsed");
      }
    }
    double e = 0.0;
    Matrix3d w = Matrix3d::Zero();
    for (long i = 0; i < 3 * nAtoms; ++i) {
      forces[i] = 0.0;
    }
    for (long i = 0; i < nAtoms; ++i) {
      for (long j = i; j < nAtoms; ++j) {
        Eigen::RowVector3d d0(positions[3 * i] - positions[3 * j],
                              positions[3 * i + 1] - positions[3 * j + 1],
                              positions[3 * i + 2] - positions[3 * j + 2]);
        const Eigen::RowVector3d f = d0 * inv;
        d0 -= f.array().round().matrix() * h;
        const double weight = (i == j) ? 0.5 : 1.0;
        for (int a = -reach[0]; a <= reach[0]; ++a) {
          for (int b = -reach[1]; b <= reach[1]; ++b) {
            for (int c = -reach[2]; c <= reach[2]; ++c) {
              if (i == j && a == 0 && b == 0 && c == 0) {
                continue;
              }
              const Eigen::RowVector3d d =
                  d0 + a * h.row(0) + b * h.row(1) + c * h.row(2);
              const double r2 = d.squaredNorm();
              if (r2 >= kCutoff * kCutoff) {
                continue;
              }
              const double r = std::sqrt(r2);
              double phi = 0.0, dphi = 0.0;
              pair(r, phi, dphi);
              e += weight * phi;
              const double g = weight * dphi / r;
              for (int k = 0; k < 3; ++k) {
                forces[3 * i + k] -= g * d(k);
                forces[3 * j + k] += g * d(k);
              }
              w += g * d.transpose() * d;
            }
          }
        }
      }
    }
    *energy = e;
    if (variance != nullptr) {
      *variance = 0.0;
    }
    last = w / std::abs(h.determinant());
  }
};

void place(eonc::Matter &matter, double a, double c) {
  matter.resize(2);
  matter.setAtomicNr(0, 18);
  matter.setAtomicNr(1, 18);
  Matrix3d cell = Matrix3d::Zero();
  cell(0, 0) = a;
  cell(1, 1) = a;
  cell(2, 2) = c;
  matter.setCell(cell);
  matter.setPeriodic(true);
  AtomMatrix positions = AtomMatrix::Zero(2, 3);
  positions(1, 0) = 0.5 * a;
  positions(1, 1) = 0.5 * a;
  positions(1, 2) = 0.5 * c;
  matter.setPositions(positions);
}

TEST_CASE("LatticeLJ stress is the strain derivative of its energy",
          "[neb][solid_state][bain]") {
  LatticeLJ pot;
  Matrix3d h;
  h << 1.2, 0.0, 0.0, 0.03, 1.19, 0.0, -0.02, 0.01, 1.31;
  AtomMatrix x = AtomMatrix::Zero(2, 3);
  x.row(1) = 0.5 * h.colwise().sum() + Eigen::RowVector3d(0.01, 0.02, 0.03);
  VectorXi z(2);
  z << 18, 18;
  pot.get_ef(x, z, h);
  const Matrix3d sigma = pot.last;
  const double volume = std::abs(h.determinant());
  const int ij[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 0}, {2, 0}, {2, 1}};
  for (const auto &m : ij) {
    Matrix3d eps = Matrix3d::Zero();
    eps(m[0], m[1]) = 1e-6;
    const Matrix3d grow = Matrix3d::Identity() + eps;
    const Matrix3d shrink = Matrix3d::Identity() - eps;
    const double up = std::get<0>(pot.get_ef(x * grow, z, h * grow));
    const double down = std::get<0>(pot.get_ef(x * shrink, z, h * shrink));
    REQUIRE_THAT((up - down) / 2e-6 / volume,
                 Catch::Matchers::WithinAbs(sigma(m[0], m[1]), 1e-6));
  }
}

TEST_CASE("the solid-state band climbs the Bain path of an LJ crystal",
          "[neb][solid_state][bain]") {
  Parameters params;
  ParametersLoadAccess::main_options(params).parallel = false;
  auto &neb = ParametersLoadAccess::neb_options(params);
  neb.image_count = 7;
  neb.max_iterations = 5000;
  neb.opt_method = eonc::OptType::FIRE;
  neb.force_tolerance = 1e-4;
  neb.solid_state.enabled = true;
  neb.climbing_image.enabled = true;
  // Arm the climbing image once the band force is under 0.5 eV/Angstrom
  // (ci_after), rgsaddle's default.
  neb.climbing_image.trigger_force = 0.5;
  auto &opt = ParametersLoadAccess::optimizer_options(params);
  opt.method = eonc::OptType::FIRE;
  opt.max_move = 0.1;
  opt.max_iterations = 5000;
  opt.converged_force = 1e-4;
  opt.convergence_metric = "max_component";

  auto pot = std::make_shared<LatticeLJ>();
  auto reactant = std::make_shared<eonc::Matter>(pot, params);
  auto product = std::make_shared<eonc::Matter>(pot, params);
  place(*reactant, kBccA, kBccA);
  place(*product, kFccA, kFccC);
  REQUIRE_THAT(reactant->getPotentialEnergy() / 2.0,
               Catch::Matchers::WithinAbs(kBccEnergy, 1e-9));
  const size_t callsBefore = pot->forceCallCounter.load();

  eonc::NudgedElasticBand band(reactant, product, params, pot);
  const auto status = band.compute();

  const double e0 = band.path.front()->getPotentialEnergy();
  for (long i = 0; i <= band.numImages + 1; ++i) {
    const auto &image = *band.path[static_cast<size_t>(i)];
    const Matrix3d h = image.getCell();
    std::printf("image %ld c/a %.8f v %.10f E - E_bcc %.10f%s\n", i,
                h(2, 2) / h(0, 0), std::abs(h.determinant()) / 2.0,
                (image.getPotentialEnergy() - e0) / 2.0,
                i == band.climbingImage ? " climbing" : "");
  }
  const long ci = band.climbingImage;
  const Matrix3d hci = band.path[static_cast<size_t>(ci)]->getCell();
  const double q = hci(2, 2) / hci(0, 0);
  const double barrier =
      (band.path[static_cast<size_t>(ci)]->getPotentialEnergy() - e0) / 2.0;
  std::printf("eOn Bain band: status %d iterations %ld force_calls %zu "
              "barrier %.10f eV/atom c/a %.8f\n",
              static_cast<int>(status), band.lastIteration(),
              pot->forceCallCounter.load() - callsBefore, barrier, q);

  REQUIRE(status == eonc::NudgedElasticBand::NEBStatus::GOOD);
  REQUIRE(ci >= 1);
  REQUIRE(ci <= band.numImages);
  REQUIRE(std::abs(q - kSaddleQ) <= kScanSpacing);
  REQUIRE(std::abs(barrier - kBarrier) <= kBarrierTol);

  // No rotation, no shear drift, a square base, and no atom off its
  // site: up to the lattice translation that wrapping into the cell
  // adds, every image keeps the reactant's fractional coordinates.
  const AtomMatrix s0 =
      reactant->getPositions() * reactant->getCell().inverse();
  for (long i = 0; i <= band.numImages + 1; ++i) {
    const auto &image = *band.path[static_cast<size_t>(i)];
    const Matrix3d h = image.getCell();
    REQUIRE(h(0, 1) == 0.0);
    REQUIRE(h(0, 2) == 0.0);
    REQUIRE(h(1, 2) == 0.0);
    REQUIRE(std::abs(h(1, 0)) < 1e-10);
    REQUIRE(std::abs(h(2, 0)) < 1e-10);
    REQUIRE(std::abs(h(2, 1)) < 1e-10);
    REQUIRE(std::abs(h(0, 0) - h(1, 1)) < 1e-10);
    const AtomMatrix ds = image.getPositions() * h.inverse() - s0;
    REQUIRE((ds.array() - ds.array().round()).abs().maxCoeff() < 1e-6);
  }
}

} // namespace
} // namespace tests
