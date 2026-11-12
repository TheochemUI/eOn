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
#include <iomanip>
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

// E = -x^2 + s y^2 + z^2 + 1/2 (h_xx - 5)^2. s = +1 is index 1 at the
// origin; s = -1 is index 2. The cell spring is positive. A second kind
// is the two-atom surface whose only negative curvature is a rigid
// translation, or a relative mode once that translation is removed.
class IndexSurface final : public eonc::Potential {
public:
  enum class Kind { Saddle, Double, Translated };

  IndexSurface(Kind kind, bool negativeRelative)
      : Potential(eonc::PotType::LJ), kind_(kind),
        negativeRelative_(negativeRelative) {}

  [[nodiscard]] bool computesStress() const noexcept override { return true; }
  [[nodiscard]] Matrix3d cauchyStress() const override { return sigma_; }
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  void force(long N, const double *R, const int *, double *F, double *U,
             double *, const double *box) override {
    Matrix3d h = Matrix3d::Zero();
    for (int i = 0; i < 9; ++i) {
      h(i / 3, i % 3) = box[i];
    }
    energyOf(R, N, h, U, F);
    sigma_.setZero();
    const double volume = std::abs(h.determinant());
    if (!(volume > 0.0)) {
      return;
    }
    const double step = 1e-7;
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        Matrix3d eps = Matrix3d::Zero();
        eps(a, b) = step;
        const Matrix3d plus = h * (Matrix3d::Identity() + eps);
        const Matrix3d minus = h * (Matrix3d::Identity() - eps);
        std::vector<double> rPlus(static_cast<size_t>(3 * N));
        std::vector<double> rMinus(static_cast<size_t>(3 * N));
        const Matrix3d deformPlus = Matrix3d::Identity() + eps;
        const Matrix3d deformMinus = Matrix3d::Identity() - eps;
        for (long atom = 0; atom < N; ++atom) {
          Eigen::RowVector3d row(R[3 * atom], R[3 * atom + 1], R[3 * atom + 2]);
          Eigen::RowVector3d up = row * deformPlus;
          Eigen::RowVector3d down = row * deformMinus;
          for (int c = 0; c < 3; ++c) {
            rPlus[static_cast<size_t>(3 * atom + c)] = up[c];
            rMinus[static_cast<size_t>(3 * atom + c)] = down[c];
          }
        }
        double ePlus = 0.0;
        double eMinus = 0.0;
        energyOf(rPlus.data(), N, plus, &ePlus, nullptr);
        energyOf(rMinus.data(), N, minus, &eMinus, nullptr);
        sigma_(a, b) = (ePlus - eMinus) / (2.0 * step) / volume;
      }
    }
  }

private:
  Kind kind_;
  bool negativeRelative_;
  Matrix3d sigma_{Matrix3d::Zero()};

  void energyOf(const double *R, long N, const Matrix3d &h, double *energy,
                double *force) const {
    if (kind_ == Kind::Translated) {
      REQUIRE(N == 2);
      const double s2 = std::sqrt(2.0);
      const double cx = (R[0] + R[3]) / s2;
      const double ux = (R[0] - R[3]) / s2;
      const double uy = (R[1] - R[4]) / s2;
      const double sign = negativeRelative_ ? -1.0 : 1.0;
      *energy = sign * 0.5 * ux * ux - 0.5 * cx * cx + 0.5 * uy * uy;
      if (force == nullptr) {
        return;
      }
      for (int i = 0; i < 6; ++i) {
        force[i] = 0.0;
      }
      const double dCx = -cx;
      const double dUx = sign * ux;
      const double dUy = uy;
      const double g0 = (dCx + dUx) / s2;
      const double g3 = (dCx - dUx) / s2;
      const double g1 = dUy / s2;
      const double g4 = -dUy / s2;
      force[0] = -g0;
      force[3] = -g3;
      force[1] = -g1;
      force[4] = -g4;
      return;
    }
    REQUIRE(N == 1);
    const double ySign = kind_ == Kind::Double ? -1.0 : 1.0;
    const double x = R[0] - 2.5;
    const double y = R[1] - 2.5;
    const double z = R[2] - 2.5;
    const double dh = h(0, 0) - 5.0;
    *energy = -x * x + ySign * y * y + z * z + 0.5 * dh * dh;
    if (force == nullptr) {
      return;
    }
    force[0] = 2.0 * x;
    force[1] = -2.0 * ySign * y;
    force[2] = -2.0 * z;
  }
};

void placeOrigin(eonc::Matter &matter, long atoms) {
  matter.resize(atoms);
  for (long i = 0; i < atoms; ++i) {
    matter.setAtomicNr(i, 1);
  }
  matter.setCell(Matrix3d::Identity() * 5.0);
  matter.setPeriodic(true);
  AtomMatrix positions = AtomMatrix::Zero(atoms, 3);
  if (atoms == 1) {
    positions.row(0).setConstant(2.5);
  }
  matter.setPositions(positions);
}

// Analytic Mishin EAM1 (Phys. Rev. B 63, 224106, Table I). The setfl
// file is a later tabulation; this is the form the reference saddle uses.
class MishinEam1 final : public eonc::Potential {
public:
  MishinEam1() : Potential(eonc::PotType::LJ) {}

  [[nodiscard]] bool computesStress() const noexcept override { return true; }
  [[nodiscard]] Matrix3d cauchyStress() const override { return sigma_; }
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  void force(long N, const double *R, const int *, double *F, double *U,
             double *, const double *box) override {
    Matrix3d cell = Matrix3d::Zero();
    for (int i = 0; i < 9; ++i) {
      cell(i / 3, i % 3) = box[i];
    }
    const double rc = 5.50679;
    const double hcut = 0.50037;
    const double lengths[3] = {cell.row(0).norm(), cell.row(1).norm(),
                               cell.row(2).norm()};
    const double minLen = std::min(lengths[0], std::min(lengths[1], lengths[2]));
    const int nimg = std::max(1, static_cast<int>(std::ceil(rc / minLen)));
    const double rc2 = rc * rc;
    struct Bond {
      int i;
      int j;
      Eigen::Vector3d d;
      double r;
    };
    std::vector<Bond> bonds;
    auto positive = [](int sx, int sy, int sz) {
      return sx > 0 || (sx == 0 && (sy > 0 || (sy == 0 && sz > 0)));
    };
    for (int i = 0; i < static_cast<int>(N); ++i) {
      const Eigen::Vector3d ri(R[3 * i], R[3 * i + 1], R[3 * i + 2]);
      for (int j = i; j < static_cast<int>(N); ++j) {
        const Eigen::Vector3d rj(R[3 * j], R[3 * j + 1], R[3 * j + 2]);
        for (int sx = -nimg; sx <= nimg; ++sx) {
          for (int sy = -nimg; sy <= nimg; ++sy) {
            for (int sz = -nimg; sz <= nimg; ++sz) {
              if (i == j && !positive(sx, sy, sz)) {
                continue;
              }
              const Eigen::Vector3d shift = sx * cell.row(0).transpose() +
                                            sy * cell.row(1).transpose() +
                                            sz * cell.row(2).transpose();
              const Eigen::Vector3d d = rj + shift - ri;
              const double r2 = d.squaredNorm();
              if (r2 > 1e-16 && r2 < rc2) {
                bonds.push_back(Bond{i, j, d, std::sqrt(r2)});
              }
            }
          }
        }
      }
    }
    auto weight = [&](double r) {
      const double x = (r - rc) / hcut;
      if (x >= 0.0) {
        return std::pair<double, double>{0.0, 0.0};
      }
      const double x2 = x * x;
      const double x4 = x2 * x2;
      const double den = 1.0 + x4;
      return std::pair<double, double>{x4 / den,
                                       4.0 * x * x2 / (den * den) / hcut};
    };
    auto morse = [](double r, double r0, double alpha) {
      const double z = -alpha * (r - r0);
      const double e1 = std::exp(z);
      const double e2 = std::exp(2.0 * z);
      return std::pair<double, double>{e2 - 2.0 * e1, 2.0 * alpha * (e1 - e2)};
    };
    auto pair = [&](double r) {
      if (r >= rc) {
        return std::pair<double, double>{0.0, 0.0};
      }
      const auto [m1, dm1] = morse(r, 0.83591, 2.97758);
      const auto [m2, dm2] = morse(r, 4.46867, 1.54927);
      const auto [psi, dpsi] = weight(r);
      const double core = 2.01458e2 * m1 + 6.59288e-3 * m2 + 0.86225e-2;
      double value = core * psi;
      double slope = (2.01458e2 * dm1 + 6.59288e-3 * dm2) * psi + core * dpsi;
      const double rs[3] = {2.24, 1.8, 1.2};
      const double strength[3] = {4.0, 40.0, 1.15e3};
      for (int k = 0; k < 3; ++k) {
        if (r < rs[k]) {
          const double gap = rs[k] - r;
          const double gap3 = gap * gap * gap;
          value -= strength[k] * gap * gap3;
          slope += 4.0 * strength[k] * gap3;
        }
      }
      return std::pair<double, double>{value, slope};
    };
    auto density = [&](double r) {
      if (r >= rc) {
        return std::pair<double, double>{0.0, 0.0};
      }
      const auto [psi, dpsi] = weight(r);
      const double d3 = r - (-2.19885);
      const double g1 = std::exp(-0.17394 * d3 * d3);
      const double dg1 = g1 * (-2.0 * 0.17394 * d3);
      const double arg = -5.35661e2 * (r - (-2.61984e2));
      double g2 = 0.0;
      double dg2 = 0.0;
      if (arg >= -700.0) {
        g2 = std::exp(arg);
        dg2 = -5.35661e2 * g2;
      }
      const double inner = 3.80362 * g1 + g2;
      const double dinner = 3.80362 * dg1 + dg2;
      return std::pair<double, double>{inner * psi, dinner * psi + inner * dpsi};
    };
    auto embed = [](double rho) {
      const double x = rho - 1.0;
      const double x2 = x * x;
      const double x3 = x2 * x;
      if (rho < 1.0) {
        const double value = -2.28235 + 0.5 * 1.35535 * x2 + (-1.27775) * x3 +
                             (-0.86074) * x2 * x2 + 1.78804 * x2 * x3 +
                             2.97571 * x3 * x3;
        const double slope = 1.35535 * x + 3.0 * (-1.27775) * x2 +
                             4.0 * (-0.86074) * x3 + 5.0 * 1.78804 * x2 * x2 +
                             6.0 * 2.97571 * x2 * x3;
        return std::pair<double, double>{value, slope};
      }
      const double num =
          -2.28235 + 0.5 * 1.35535 * x2 + (-1.27775) * x3 + 0.4 * x2 * x2;
      const double den = 1.0 + 0.3 * x3;
      const double dnum = 1.35535 * x + 3.0 * (-1.27775) * x2 + 4.0 * 0.4 * x3;
      const double dden = 3.0 * 0.3 * x2;
      return std::pair<double, double>{
          num / den, (dnum * den - num * dden) / (den * den)};
    };
    std::vector<double> rho(static_cast<size_t>(N), 0.0);
    for (const Bond &bond : bonds) {
      const double w = density(bond.r).first;
      rho[static_cast<size_t>(bond.i)] += w;
      if (bond.i != bond.j) {
        rho[static_cast<size_t>(bond.j)] += w;
      } else {
        rho[static_cast<size_t>(bond.i)] += w;
      }
    }
    double energy = 0.0;
    std::vector<double> slope(static_cast<size_t>(N), 0.0);
    for (int i = 0; i < static_cast<int>(N); ++i) {
      const auto [value, derivative] = embed(rho[static_cast<size_t>(i)]);
      energy += value;
      slope[static_cast<size_t>(i)] = derivative;
    }
    for (int i = 0; i < static_cast<int>(N); ++i) {
      F[3 * i] = 0.0;
      F[3 * i + 1] = 0.0;
      F[3 * i + 2] = 0.0;
    }
    Matrix3d virial = Matrix3d::Zero();
    for (const Bond &bond : bonds) {
      const auto [v, dv] = pair(bond.r);
      const double drho = density(bond.r).second;
      const double invR = 1.0 / bond.r;
      if (bond.i != bond.j) {
        energy += v;
        const double deriv = dv + (slope[static_cast<size_t>(bond.i)] +
                                    slope[static_cast<size_t>(bond.j)]) *
                                       drho;
        const Eigen::Vector3d f = deriv * bond.d * invR;
        for (int c = 0; c < 3; ++c) {
          F[3 * bond.i + c] += f[c];
          F[3 * bond.j + c] -= f[c];
        }
        virial += deriv * (bond.d * bond.d.transpose()) * invR;
      } else {
        energy += v;
        const double deriv =
            dv + 2.0 * slope[static_cast<size_t>(bond.i)] * drho;
        virial += deriv * (bond.d * bond.d.transpose()) * invR;
      }
    }
    *U = energy;
    const double volume = std::abs(cell.determinant());
    sigma_.setZero();
    if (volume > 0.0) {
      sigma_ = virial / volume;
    }
  }

private:
  Matrix3d sigma_{Matrix3d::Zero()};
};

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

TEST_CASE("Mishin EAM1 crystal matches the analytic cohesive energy",
          "[solid_state][eam_cu]") {
  auto params = baseParams();
  auto pot = std::make_shared<MishinEam1>();
  eonc::Matter crystal(pot, params);
  crystal.resize(4);
  const double a = 3.615;
  const double basis[4][3] = {{0, 0, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5},
                              {0.5, 0.5, 0}};
  AtomMatrix positions(4, 3);
  for (int b = 0; b < 4; ++b) {
    crystal.setAtomicNr(b, 29);
    positions(b, 0) = basis[b][0] * a;
    positions(b, 1) = basis[b][1] * a;
    positions(b, 2) = basis[b][2] * a;
  }
  crystal.setCell(Matrix3d::Identity() * a);
  crystal.setPeriodic(true);
  crystal.setPositions(positions);
  const double perAtom = crystal.getPotentialEnergy() / 4.0;
  INFO("E/atom " << std::setprecision(17) << perAtom);
  REQUIRE(perAtom == Catch::Approx(-3.5402183104868774).margin(1e-12));
  REQUIRE(crystal.getForces().norm() < 1e-6);
  REQUIRE(crystal.cauchyStress().cwiseAbs().maxCoeff() < 1e-4);
  const auto path = mishinFile();
  if (std::filesystem::exists(path)) {
    auto tabulated = baseParams();
    ParametersLoadAccess::potential_options(tabulated).eamFile = path.string();
    auto setfl = std::make_shared<eonc::SetflEam>(tabulated, path.string());
    constexpr int nFit = 4;
    eonc::Matter fileCrystal(setfl, tabulated);
    fileCrystal.resize(4 * nFit * nFit * nFit);
    AtomMatrix big(4 * nFit * nFit * nFit, 3);
    int placed = 0;
    for (int ix = 0; ix < nFit; ++ix) {
      for (int iy = 0; iy < nFit; ++iy) {
        for (int iz = 0; iz < nFit; ++iz) {
          for (int b = 0; b < 4; ++b, ++placed) {
            fileCrystal.setAtomicNr(placed, 29);
            big(placed, 0) = (ix + basis[b][0]) * a;
            big(placed, 1) = (iy + basis[b][1]) * a;
            big(placed, 2) = (iz + basis[b][2]) * a;
          }
        }
      }
    }
    fileCrystal.setCell(Matrix3d::Identity() * (nFit * a));
    fileCrystal.setPeriodic(true);
    fileCrystal.setPositions(big);
    const double setflPerAtom = fileCrystal.getPotentialEnergy() / big.rows();
    std::ofstream note("/tmp/cu-setfl.txt");
    note << std::setprecision(17) << setflPerAtom << "\n";
    INFO("setfl E/atom " << std::setprecision(17) << setflPerAtom);
  }
}

TEST_CASE("index check rejects an index-2 point and accepts index 1",
          "[solid_state][dimer][index_joint]") {
  auto params = baseParams();
  ParametersLoadAccess::solid_state_options(params).enabled = true;
  auto mask = ParametersLoadAccess::solid_state_options(params).mask;
  mask.fill(0);
  mask[0] = 1;
  ParametersLoadAccess::solid_state_options(params).mask = mask;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 1e-8;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 6;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.05;

  auto load = [&](IndexSurface::Kind kind) {
    auto pot = std::make_shared<IndexSurface>(kind, false);
    auto matter = std::make_shared<eonc::Matter>(pot, params);
    placeOrigin(*matter, 1);
    return std::pair{pot, matter};
  };
  {
    auto [pot, matter] = load(IndexSurface::Kind::Saddle);
    eonc::solid::Chart chart(*matter, eonc::solid::optionsFrom(params));
    const auto accepted = eonc::solid::indexAt(*matter, chart);
    INFO("accepted " << accepted.values[0] << " " << accepted.values[1]);
    REQUIRE(accepted.isFirstOrder());
    REQUIRE(accepted.negative == 1);
    REQUIRE(accepted.values[0] < -1.5);
    REQUIRE(accepted.values[1] > 0.5);
    REQUIRE(accepted.evaluations > 0);
  }
  {
    auto [pot, matter] = load(IndexSurface::Kind::Double);
    eonc::solid::Chart chart(*matter, eonc::solid::optionsFrom(params));
    const auto rejected = eonc::solid::indexAt(*matter, chart);
    INFO("rejected " << rejected.values[0] << " " << rejected.values[1]);
    REQUIRE(rejected.negative == 2);
    REQUIRE_FALSE(rejected.isFirstOrder());
    REQUIRE(rejected.values[0] < -1.5);
    REQUIRE(rejected.values[1] < -1.5);
  }
  {
    auto [pot, matter] = load(IndexSurface::Kind::Saddle);
    AtomMatrix mode(1, 3);
    mode.setZero();
    mode(0, 0) = 1.0;
    eonc::MinModeSaddleSearch search(matter, mode, 0.0, params, pot);
    const int status = search.run();
    INFO("saddle status " << status << " n " << search.negativeModes()
                          << " x " << matter->getPositions());
    REQUIRE(status == eonc::MinModeSaddleSearch::STATUS_GOOD);
    REQUIRE(search.negativeModes() == 1);
    REQUIRE_FALSE(search.saddleRejected());
    const Eigen::RowVector3d centre =
        matter->getPositions().row(0) - Eigen::RowVector3d::Constant(2.5);
    REQUIRE(centre.norm() < 1e-5);
    REQUIRE(search.variableCell());
  }
  {
    auto frozen = params;
    auto frozenMask = ParametersLoadAccess::solid_state_options(frozen).mask;
    frozenMask.fill(0);
    ParametersLoadAccess::solid_state_options(frozen).mask = frozenMask;
    auto pot = std::make_shared<IndexSurface>(IndexSurface::Kind::Saddle, false);
    auto matter = std::make_shared<eonc::Matter>(pot, frozen);
    placeOrigin(*matter, 1);
    auto before = eonc::solid::indexAtoms(*matter);
    INFO("direct n " << before.negative << " v0 "
                     << (before.values.empty() ? 0.0 : before.values.front())
                     << " v1 "
                     << (before.values.size() < 2 ? 0.0 : before.values[1]));
    REQUIRE(before.negative == 1);
    AtomMatrix mode(1, 3);
    mode.setZero();
    mode(0, 0) = 1.0;
    eonc::MinModeSaddleSearch search(matter, mode, 0.0, frozen, pot);
    const int status = search.run();
    INFO("frozen status " << status << " n " << search.negativeModes());
    REQUIRE(status == eonc::MinModeSaddleSearch::STATUS_GOOD);
    REQUIRE_FALSE(search.variableCell());
    REQUIRE(search.negativeModes() == 1);
    REQUIRE(search.getIndexForceCalls() > 0);
    const Eigen::RowVector3d frozenCentre =
        matter->getPositions().row(0) - Eigen::RowVector3d::Constant(2.5);
    REQUIRE(frozenCentre.norm() < 1e-5);
  }
  {
    auto [pot, matter] = load(IndexSurface::Kind::Double);
    AtomMatrix mode(1, 3);
    mode.setZero();
    mode(0, 0) = 1.0;
    eonc::MinModeSaddleSearch search(matter, mode, 0.0, params, pot);
    const int status = search.run();
    const double moved =
        (matter->getPositions().row(0) - Eigen::RowVector3d::Constant(2.5))
            .norm();
    INFO("double status " << status << " n " << search.negativeModes()
                          << " moved " << moved);
    const bool stillThere = moved < 1e-6;
    if (stillThere) {
      REQUIRE(status == eonc::MinModeSaddleSearch::STATUS_BAD_NO_BARRIER);
      REQUIRE(search.saddleRejected());
      REQUIRE(search.negativeModes() == 2);
    } else {
      const bool kept =
          status != eonc::MinModeSaddleSearch::STATUS_GOOD ||
          search.negativeModes() == 1;
      REQUIRE(kept);
    }
  }
}

TEST_CASE("index check does not count a rigid translation",
          "[solid_state][dimer][index_translation]") {
  auto params = baseParams();
  ParametersLoadAccess::main_options(params).removeNetForce = false;
  auto bare = std::make_shared<IndexSurface>(IndexSurface::Kind::Translated, false);
  eonc::Matter matter(bare, params);
  placeOrigin(matter, 2);
  matter.setPeriodic(false);
  eonc::solid::IndexCheck projected;
  auto hidden = eonc::solid::indexAtoms(matter, projected);
  INFO("projected " << hidden.negative << " lowest "
                    << (hidden.values.empty() ? 0.0 : hidden.values.front()));
  REQUIRE(hidden.negative == 0);
  eonc::solid::IndexCheck raw = projected;
  raw.project = false;
  auto visible = eonc::solid::indexAtoms(matter, raw);
  INFO("raw " << visible.negative << " lowest "
              << (visible.values.empty() ? 0.0 : visible.values.front()));
  REQUIRE(visible.negative == 1);
  REQUIRE(visible.values.front() < -0.5);

  auto saddlePot =
      std::make_shared<IndexSurface>(IndexSurface::Kind::Translated, true);
  eonc::Matter saddle(saddlePot, params);
  placeOrigin(saddle, 2);
  saddle.setPeriodic(false);
  auto one = eonc::solid::indexAtoms(saddle, projected);
  INFO("relative " << one.negative << " lowest "
                   << (one.values.empty() ? 0.0 : one.values.front()));
  REQUIRE(one.isFirstOrder());
  REQUIRE(one.negative == 1);
  REQUIRE(one.values.front() < -0.5);
}

TEST_CASE("copper dumbbell saddle follows the shear-aligned switch",
          "[solid_state][eam_cu][dimer]") {
  auto params = baseParams();
  ParametersLoadAccess::optimizer_options(params).method = eonc::OptType::LBFGS;
  ParametersLoadAccess::optimizer_options(params).converged_force = 2e-3;
  ParametersLoadAccess::optimizer_options(params).convergence_metric = "max_component";
  ParametersLoadAccess::optimizer_options(params).max_iterations = 600;
  ParametersLoadAccess::optimizer_options(params).max_move = 0.05;
  ParametersLoadAccess::saddle_search_options(params).converged_force = 2e-3;
  ParametersLoadAccess::saddle_search_options(params).max_iterations = 2500;
  ParametersLoadAccess::saddle_search_options(params).max_energy = 1.2;
  ParametersLoadAccess::dimer_options(params).torque_min = 1e-3;
  ParametersLoadAccess::dimer_options(params).rotations_max = 8;
  auto pot = std::make_shared<MishinEam1>();

  constexpr double kA = 3.615;
  constexpr int kN = 3;
  constexpr double kPi = 3.14159265358979323846;
  const double length = kN * kA;
  const double basis[4][3] = {{0, 0, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5},
                              {0.5, 0.5, 0}};
  std::vector<Eigen::RowVector3d> perfect;
  for (int ix = 0; ix < kN; ++ix) {
    for (int iy = 0; iy < kN; ++iy) {
      for (int iz = 0; iz < kN; ++iz) {
        for (int b = 0; b < 4; ++b) {
          perfect.emplace_back((ix + basis[b][0]) * kA, (iy + basis[b][1]) * kA,
                               (iz + basis[b][2]) * kA);
        }
      }
    }
  }
  REQUIRE(perfect.size() == 108);
  auto build = [&](const Eigen::Vector3d &direction) {
    std::vector<Eigen::RowVector3d> atoms(perfect.begin() + 1, perfect.end());
    atoms.push_back(-0.5 * 2.05 * direction.transpose());
    atoms.push_back(0.5 * 2.05 * direction.transpose());
    return atoms;
  };
  const auto along100 = build(Eigen::Vector3d(1.0, 0.0, 0.0));
  const auto along110 =
      build(Eigen::Vector3d(1.0, 1.0, 0.0).normalized());
  REQUIRE(along100.size() == 109);

  auto shearOf = [&](const std::vector<Eigen::RowVector3d> &atoms, double gamma) {
    Matrix3d deform = Matrix3d::Identity();
    deform(1, 0) = gamma;
    Matrix3d cell = (Matrix3d::Identity() * length) * deform;
    AtomMatrix positions(static_cast<long>(atoms.size()), 3);
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
      positions.row(i) = atoms[static_cast<size_t>(i)] * deform;
    }
    return std::pair<AtomMatrix, Matrix3d>{positions, cell};
  };
  auto load = [&](const AtomMatrix &positions, const Matrix3d &cell,
                  const eonc::Parameters &local) {
    auto matter = std::make_shared<eonc::Matter>(pot, local);
    matter->resize(positions.rows());
    for (int i = 0; i < positions.rows(); ++i) {
      matter->setAtomicNr(i, 29);
    }
    matter->setCell(cell);
    matter->setPeriodic(true);
    matter->setPositions(positions);
    return matter;
  };
  auto mic = [](const Eigen::Vector3d &delta, const Matrix3d &cell) {
    Eigen::RowVector3d frac = delta.transpose() * cell.inverse();
    for (int c = 0; c < 3; ++c) {
      frac[c] -= std::round(frac[c]);
    }
    return Eigen::Vector3d((frac * cell).transpose());
  };
  auto pairAxis = [&](const eonc::Matter &matter) {
    const AtomMatrix pos = matter.getPositions();
    const Matrix3d cell = matter.getCell();
    double best = 1e9;
    Eigen::Vector3d axis = Eigen::Vector3d::Zero();
    for (int i = 0; i < pos.rows(); ++i) {
      for (int j = i + 1; j < pos.rows(); ++j) {
        const Eigen::Vector3d d = mic(
            (pos.row(j) - pos.row(i)).transpose(), cell);
        if (d.norm() < best) {
          best = d.norm();
          axis = d;
        }
      }
    }
    return axis;
  };
  auto family = [](const Eigen::Vector3d &u, std::initializer_list<Eigen::Vector3d> dirs) {
    double best = 0.0;
    for (const auto &dir : dirs) {
      best = std::max(best, std::abs(u.normalized().dot(dir.normalized())));
    }
    return std::acos(std::clamp(best, 0.0, 1.0)) * 180.0 / kPi;
  };
  auto angles = [&](const Eigen::Vector3d &axis) {
    const double to100 =
        family(axis, {Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0),
                      Eigen::Vector3d(0, 0, 1)});
    const double to110 = family(axis, {{1, 1, 0}, {1, -1, 0}, {1, 0, 1},
                                      {1, 0, -1}, {0, 1, 1}, {0, 1, -1}});
    return std::pair<double, double>{to100, to110};
  };
  auto enthalpyOf = [](double energy, const Matrix3d &h0, const Matrix3d &h,
                       const Matrix3d &external) {
    Matrix3d eps = h0.inverse() * h - Matrix3d::Identity();
    eps(0, 1) = 0.0;
    eps(0, 2) = 0.0;
    eps(1, 2) = 0.0;
    double work = 0.0;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j <= i; ++j) {
        work += external(i, j) * eps(i, j);
      }
    }
    return energy + std::abs(h0.determinant()) * work;
  };

  const double fixedRef[] = {0.473916239370908, 0.24356279471709286,
                             0.14826689300639373, 0.06784288866265342,
                             0.00776114637324099};
  const double enthalpyRef[] = {0.34702642693292773, 0.1630119173183857,
                                0.08761266284551539, 0.0274332017463621};
  const double gammas[] = {0.0, 0.02, 0.03, 0.04, 0.05};
  std::ofstream table("/tmp/cu-saddle.txt");
  table << std::setprecision(17);
  table << "shear,cell,calls,index_calls,barrier_eV,enthalpy_eV,curvature,"
           "max_force,n_neg,reactant_n_neg,reactant_lowest,status,rejected,"
           "switch,result\n";
  int fixedSlot = 0;
  int variableSlot = 0;
  for (double gamma : gammas) {
    const auto [pos100, cell] = shearOf(along100, gamma);
    const auto [pos110, unusedCell] = shearOf(along110, gamma);
    auto minimum = load(pos100, cell, params);
    auto other = load(pos110, cell, params);
    const bool minOk = minimum->relax(true);
    INFO("minimum max force " << minimum->getForces().cwiseAbs().maxCoeff()
                              << " E " << minimum->getPotentialEnergy());
    REQUIRE(minOk);
    const bool otherOk = other->relax(true);
    INFO("110 max force " << other->getForces().cwiseAbs().maxCoeff());
    REQUIRE(otherOk);
    const double eMin = minimum->getPotentialEnergy();
    const Matrix3d sigmaMin = minimum->cauchyStress();
    const double volumeMin = std::abs(minimum->getCell().determinant());
    const double pMin = volumeMin * sigmaMin(1, 0);
    const Eigen::Vector3d minAxis = pairAxis(*minimum);
    const auto [min100, min110] = angles(minAxis);
    INFO("shear " << gamma << " min100 " << min100 << " min110 " << min110
                  << " E " << eMin << " stress_xy " << sigmaMin(1, 0));
    if (min100 > 25.0) {
      table << gamma << " fixed 0 0 0 " << sigmaMin(1, 0) << " " << pMin
            << " " << pMin
            << " 0 0 -1 " << min100
            << " nan nan minimum_left_<100>\n";
      table.flush();
      FAIL("the <100> minimum did not survive relaxation at shear " << gamma
           << " (angle from <100> is " << min100 << ")");
    }
    AtomMatrix connect(pos100.rows(), 3);
    for (int i = 0; i < pos100.rows(); ++i) {
      connect.row(i) = mic((other->getPositions().row(i) -
                            minimum->getPositions().row(i))
                               .transpose(),
                           cell)
                           .transpose();
    }
    const double norm = std::sqrt(connect.squaredNorm());
    REQUIRE(norm > 0.1);
    const AtomMatrix savedMin = minimum->getPositions();
    auto directionCurvature = [&](double distance) {
      const double fd = 1e-3;
      minimum->setPositions(savedMin + (distance / norm) * connect);
      const double energy = minimum->getPotentialEnergy();
      minimum->setPositions(savedMin + ((distance + fd) / norm) * connect);
      const double plus = minimum->getPotentialEnergy();
      minimum->setPositions(savedMin + ((distance - fd) / norm) * connect);
      const double minus = minimum->getPotentialEnergy();
      return (plus + minus - 2.0 * energy) / (fd * fd);
    };
    double nudge = 0.02;
    if (directionCurvature(nudge) > 0.0) {
      for (double distance = 0.025; distance <= 0.06; distance += 0.005) {
        if (directionCurvature(distance) < 0.0) {
          nudge = distance;
          break;
        }
      }
    }
    minimum->setPositions(savedMin);
    AtomMatrix nudged = savedMin + (nudge / norm) * connect;
    const Matrix3d h0 = minimum->getCell();

    auto runSearch = [&](bool variable) {
      auto local = params;
      ParametersLoadAccess::optimizer_options(local).max_move = 0.005;
      ParametersLoadAccess::solid_state_options(local).enabled = true;
      if (variable) {
        for (int i = 0; i < 9; ++i) {
          ParametersLoadAccess::solid_state_options(local)
              .stress[static_cast<size_t>(i)] = -sigmaMin(i / 3, i % 3);
        }
      } else {
        auto frozen = ParametersLoadAccess::solid_state_options(local).mask;
        frozen.fill(0);
        ParametersLoadAccess::solid_state_options(local).mask = frozen;
      }
      auto reactantMatter = load(savedMin, h0, local);
      const auto reactantBefore = pot->forceCallCounter.load();
      eonc::solid::IndexSpectrum reactantSpectrum;
      if (variable) {
        eonc::solid::Chart reactantChart(
            *reactantMatter, eonc::solid::optionsFrom(local));
        reactantSpectrum =
            eonc::solid::indexAt(*reactantMatter, reactantChart);
      } else {
        reactantSpectrum = eonc::solid::indexAtoms(*reactantMatter);
      }
      const int reactantCalls = static_cast<int>(
          pot->forceCallCounter.load() - reactantBefore);
      const int reactantNeg = reactantSpectrum.negative;
      const double reactantLowest = reactantSpectrum.values.empty()
                                        ? 0.0
                                        : reactantSpectrum.values.front();
      auto matter = load(nudged, h0, local);
      eonc::MinModeSaddleSearch search(matter, connect, eMin, local, pot);
      const auto before = pot->forceCallCounter.load();
      const int status = search.run();
      const int calls =
          static_cast<int>(pot->forceCallCounter.load() - before);
      const double eSad = matter->getPotentialEnergy();
      const double barrier = eSad - eMin;
      Matrix3d external = Matrix3d::Zero();
      if (variable) {
        external = -sigmaMin;
        external(0, 1) = 0.0;
        external(0, 2) = 0.0;
        external(1, 2) = 0.0;
      }
      const double enthalpy =
          enthalpyOf(eSad, h0, matter->getCell(), external) - eMin;
      const Eigen::Vector3d axis = pairAxis(*matter);
      const auto [sad100, sad110] = angles(axis);
      const int nNeg = search.negativeModes();
      const bool rejected = search.saddleRejected();
      const bool accepted =
          status == eonc::MinModeSaddleSearch::STATUS_GOOD && !rejected &&
          nNeg == 1 && search.getEigenvalue() < 0.0;
      const bool switched = accepted && min100 < 25.0 && sad110 < 15.0 &&
                            axis[0] * axis[1] > 0.0 && sad100 > min100 + 10.0;
      std::string result = "index 1";
      if (reactantNeg > 0) {
        result = "no barrier / reactant unstable";
      } else if (!accepted) {
        result = "no barrier / no saddle";
      }
      table << gamma << "," << (variable ? "variable" : "fixed") << ","
            << calls << "," << search.getIndexForceCalls() << "," << barrier
            << "," << enthalpy << "," << search.getEigenvalue() << ","
            << matter->getForces().cwiseAbs().maxCoeff() << "," << nNeg << ","
            << reactantNeg << "," << reactantLowest << "," << status << ","
            << static_cast<int>(rejected) << "," << static_cast<int>(switched)
            << "," << result << "\n";
      table.flush();
      INFO("shear " << gamma << " variable " << variable << " status " << status
                    << " calls " << calls << " index "
                    << search.getIndexForceCalls() << " reactant_calls "
                    << reactantCalls << " barrier " << barrier << " enthalpy "
                    << enthalpy << " n_neg " << nNeg << " reactant_n_neg "
                    << reactantNeg << " reactant_lowest " << reactantLowest
                    << " switch " << switched << " result " << result
                    << " sad100 " << sad100 << " sad110 " << sad110);
      if (reactantNeg > 0) {
        CHECK(result == "no barrier / reactant unstable");
      } else if (!variable) {
        CHECK(accepted);
        CHECK(barrier == Catch::Approx(fixedRef[fixedSlot]).margin(0.005));
        CHECK(switched);
        ++fixedSlot;
      } else if (variableSlot < 4) {
        CHECK(accepted);
        CHECK(enthalpy == Catch::Approx(enthalpyRef[variableSlot]).margin(0.005));
        CHECK(switched);
        ++variableSlot;
      } else {
        CHECK(accepted);
        CHECK_FALSE(switched);
        CHECK(std::abs(enthalpy) < 1e-3);
        ++variableSlot;
      }
    };
    runSearch(false);
    runSearch(true);
  }
}

TEST_CASE("solid state ini keys stay off unless set", "[solid_state]") {
  eonc::Parameters params;
  REQUIRE_FALSE(params.solid_state_options().enabled);
  REQUIRE(params.solid_state_options().weight == Catch::Approx(1.0));
  REQUIRE(params.solid_state_options().mask[1] == 1);
  REQUIRE(params.solid_state_options().max_strain_step == Catch::Approx(0.05));
  REQUIRE(params.solid_state_options().max_strain == Catch::Approx(0.20));
  REQUIRE(params.solid_state_options().index_check);
  REQUIRE(params.solid_state_options().flat == Catch::Approx(0.005));
  const char *text = "[Solid State]\n"
                     "enabled = true\n"
                     "weight = 1.5\n"
                     "pressure = 0.02\n"
                     "stress = 0 0.1 0 0.1 0 0 0 0 0\n"
                     "mask = 1 0 1 0 0 1\n"
                     "max_strain_step = 0.04\n"
                     "max_strain = 0.15\n"
                     "index_check = false\n"
                     "flat = 0\n"
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
  REQUIRE_FALSE(params.solid_state_options().index_check);
  REQUIRE(params.solid_state_options().flat == Catch::Approx(0.0));
  REQUIRE(params.potential_options().eamFile == "cu.setfl");
  REQUIRE_FALSE(eonc::solid::movesCell(eonc::Parameters{}));
}

} // namespace tests
