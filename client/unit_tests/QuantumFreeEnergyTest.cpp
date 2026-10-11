#include "eon/QuantumFreeEnergy.h"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Hessian.h"
#include "eon/Potential.h"
#include "eon/Tunneling.h"

#include <Eigen/Eigenvalues>
#include <cmath>

namespace {

struct CurvedValley final : eonc::Potential {
  double ka{2.0};
  double kp{8.0};
  double kz{3.0};
  double alpha{0.5};

  CurvedValley()
      : Potential(eonc::PotType::LJ) {}

  void force(long nAtoms, const double *positions, const int *, double *forces,
             double *energy, double *variance, const double *) override {
    const double x = positions[0];
    const double y = positions[1];
    const double z = positions[2];
    const double u = y - alpha * x * x;
    *energy = 0.5 * ka * x * x + 0.5 * kp * u * u + 0.5 * kz * z * z;
    if (variance != nullptr) {
      *variance = 0.0;
    }
    for (long i = 0; i < nAtoms * 3; ++i) {
      forces[i] = 0.0;
    }
    const double dVdx = ka * x + kp * u * (-2.0 * alpha * x);
    const double dVdy = kp * u;
    const double dVdz = kz * z;
    forces[0] = -dVdx;
    forces[1] = -dVdy;
    forces[2] = -dVdz;
  }
};

// Three unit masses joined pairwise by springs k of rest length d, under
// the minimum image of a cubic box when there is one. The equilateral
// minimum vibrates at omega^2 = 3k (breathing) and 3k/2 (twice); its three
// translations and three rotations carry no vibration.
struct SpringTriangle final : eonc::Potential {
  double k{2.0};
  double d{1.5};

  SpringTriangle()
      : Potential(eonc::PotType::LJ) {}

  void force(long nAtoms, const double *positions, const int *, double *forces,
             double *energy, double *variance, const double *box) override {
    *energy = 0.0;
    if (variance != nullptr) {
      *variance = 0.0;
    }
    for (long i = 0; i < 3 * nAtoms; ++i) {
      forces[i] = 0.0;
    }
    const double side = box != nullptr ? box[0] : 0.0;
    for (long i = 0; i < nAtoms; ++i) {
      for (long j = i + 1; j < nAtoms; ++j) {
        Eigen::Vector3d r;
        for (int c = 0; c < 3; ++c) {
          r(c) = positions[3 * j + c] - positions[3 * i + c];
          if (side > 0.0) {
            r(c) -= side * std::round(r(c) / side);
          }
        }
        const double len = r.norm();
        *energy += 0.5 * k * (len - d) * (len - d);
        const Eigen::Vector3d g = k * (len - d) * r / len;
        for (int c = 0; c < 3; ++c) {
          forces[3 * i + c] += g(c);
          forces[3 * j + c] -= g(c);
        }
      }
    }
  }
};

AtomMatrix equilateral(double d) {
  AtomMatrix r(3, 3);
  r << 0.0, 0.0, 0.0, d, 0.0, 0.0, 0.5 * d, 0.5 * std::sqrt(3.0) * d, 0.0;
  return r;
}

std::shared_ptr<eonc::Matter>
triangleImage(const std::shared_ptr<eonc::Potential> &pot,
              const eonc::Parameters &params, const AtomMatrix &r,
              bool periodic) {
  auto m = std::make_shared<eonc::Matter>(pot, params);
  m->resize(3);
  for (long i = 0; i < 3; ++i) {
    m->setAtomicNr(i, 1);
    m->setMass(i, 1.0);
  }
  m->setCell(Eigen::Matrix3d::Identity() * 10.0);
  m->setPeriodic(periodic);
  m->setPositions(r);
  return m;
}

Eigen::Matrix3d analyticHessian(double x, double alpha, double ka, double kp,
                                double kz) {
  Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
  const double slope = 2.0 * alpha * x;
  hessian(0, 0) = ka + kp * slope * slope;
  hessian(1, 1) = kp;
  hessian(2, 2) = kz;
  hessian(0, 1) = hessian(1, 0) = -kp * slope;
  return hessian;
}

double directSum(const Eigen::Matrix3d &hessian, const Eigen::Vector3d &tangent,
                 double temperature) {
  const Eigen::Vector3d unit = tangent.normalized();
  const Eigen::Matrix3d projector =
      Eigen::Matrix3d::Identity() - unit * unit.transpose();
  const Eigen::Matrix3d projected = projector * hessian * projector;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(projected);
  const Eigen::Vector3d values = solver.eigenvalues();
  long drop = 0;
  double nearest = std::abs(values(0));
  for (long i = 1; i < 3; ++i) {
    if (std::abs(values(i)) < nearest) {
      nearest = std::abs(values(i));
      drop = i;
    }
  }
  double sum = 0.0;
  for (long i = 0; i < 3; ++i) {
    if (i == drop || !(values(i) > 0.0)) {
      continue;
    }
    const double hw = eonc::tunneling::kHbar * std::sqrt(values(i));
    sum += eonc::harmonicModeFreeEnergy(hw, temperature);
  }
  return sum;
}

} // namespace

TEST_CASE("curved valley free energy matches a direct sum over normal modes",
          "[quantum][neb]") {
  constexpr double kX = 0.4;
  constexpr double kAlpha = 0.5;
  constexpr double kA = 2.0;
  constexpr double kP = 8.0;
  constexpr double kZ = 3.0;
  const Eigen::Matrix3d hessian = analyticHessian(kX, kAlpha, kA, kP, kZ);
  const Eigen::Vector3d tangent(1.0, 2.0 * kAlpha * kX, 0.0);
  constexpr double kTemperature = 300.0;
  const double summed = directSum(hessian, tangent, kTemperature);
  const double profile =
      eonc::perpendicularHarmonicFreeEnergy(hessian, tangent, kTemperature);
  REQUIRE(profile == Catch::Approx(summed).margin(1e-10));

  const double zpe =
      eonc::perpendicularHarmonicFreeEnergy(hessian, tangent, 0.0);
  REQUIRE(zpe == Catch::Approx(directSum(hessian, tangent, 0.0)).margin(1e-12));
  constexpr double kValleyEnergy = 0.5 * kA * kX * kX;
  REQUIRE(kValleyEnergy + zpe ==
          Catch::Approx(kValleyEnergy + directSum(hessian, tangent, 0.0))
              .margin(1e-12));

  eonc::Parameters params;
  eonc::ParametersLoadAccess::hessian_options(params).fd_scheme = "central";
  eonc::ParametersLoadAccess::main_options(params).finiteDifference = 1e-5;
  auto pot = std::make_shared<CurvedValley>();
  eonc::Matter matter(pot, params);
  matter.resize(1);
  matter.setAtomicNr(0, 1);
  matter.setMass(0, 1.0);
  matter.setCell(Eigen::Matrix3d::Identity() * 20.0);
  matter.setPeriodic(false);
  AtomMatrix positions(1, 3);
  positions << kX, kAlpha * kX * kX, 0.0;
  matter.setPositions(positions);
  eonc::Hessian numerical(params, &matter);
  numerical.writeHessianFile(false);
  Eigen::VectorXi mobile(1);
  mobile << 0;
  const Eigen::MatrixXd weighted = numerical.getHessian(&matter, mobile);
  REQUIRE(weighted.rows() == 3);
  REQUIRE((weighted - hessian).norm() < 1e-6);
  REQUIRE(
      eonc::perpendicularHarmonicFreeEnergy(weighted, tangent, kTemperature) ==
      Catch::Approx(summed).margin(1e-5));
}

// A free cluster's band: the end images are minima and keep all three
// vibrations, the middle image loses the band direction, and none of the
// six rigid motions enters, periodic cell or not. Each kept eigenvalue on
// the middle image is the finite-difference Hessian's on the complement of
// the rigid motions and the tangent.
TEST_CASE("a free cluster's band free energy leaves out its rigid motions",
          "[quantum][neb]") {
  eonc::Parameters params;
  eonc::ParametersLoadAccess::hessian_options(params).fd_scheme = "central";
  eonc::ParametersLoadAccess::main_options(params).finiteDifference = 1e-5;
  auto spring = std::make_shared<SpringTriangle>();
  const double k = spring->k;
  const double d = spring->d;
  constexpr double kTemperature = 300.0;
  double minimum = 0.0;
  for (const double omega2 : {3.0 * k, 1.5 * k, 1.5 * k}) {
    minimum += eonc::harmonicModeFreeEnergy(
        eonc::tunneling::kHbar * std::sqrt(omega2), kTemperature);
  }

  AtomMatrix stretched = equilateral(d);
  stretched(1, 0) += 0.2;
  AtomMatrix moved = equilateral(d);
  moved.rowwise() += Eigen::RowVector3d(0.3, -0.1, 0.2);
  std::vector<std::shared_ptr<eonc::Matter>> path{
      triangleImage(spring, params, equilateral(d), false),
      triangleImage(spring, params, stretched, false),
      triangleImage(spring, params, moved, false)};
  // Two minima of one cluster differ by a rigid motion, so the middle
  // image takes the band's own tangent: atom 1 along x.
  std::vector<std::shared_ptr<AtomMatrix>> tangents(path.size());
  tangents[1] = std::make_shared<AtomMatrix>(AtomMatrix::Zero(3, 3));
  (*tangents[1])(1, 0) = 1.0;
  const std::vector<double> free =
      eonc::quantumFreeEnergies(path, tangents, kTemperature, params);
  REQUIRE(free.size() == 3);
  CAPTURE(free[0], free[2], minimum);
  REQUIRE(free[0] == Catch::Approx(minimum).margin(1e-6));
  REQUIRE(free[2] == Catch::Approx(minimum).margin(1e-6));

  eonc::Hessian numerical(params, path[1].get());
  numerical.writeHessianFile(false);
  Eigen::VectorXi mobile(3);
  mobile << 0, 1, 2;
  const Eigen::MatrixXd h = numerical.getHessian(path[1].get(), mobile);
  Eigen::MatrixXd removed = Eigen::MatrixXd::Zero(9, 7);
  Eigen::RowVector3d com = stretched.colwise().mean();
  for (long a = 0; a < 3; ++a) {
    const Eigen::Vector3d x = (stretched.row(a) - com).transpose();
    for (int c = 0; c < 3; ++c) {
      removed(3 * a + c, c) = 1.0;
      Eigen::Vector3d axis = Eigen::Vector3d::Zero();
      axis(c) = 1.0;
      removed.block(3 * a, 3 + c, 3, 1) = axis.cross(x);
    }
  }
  removed(3, 6) = 1.0;
  const Eigen::VectorXd kept = eonc::complementEigenvalues(h, removed);
  REQUIRE(kept.size() == 2);
  double middle = path[1]->getPotentialEnergy();
  for (long i = 0; i < kept.size(); ++i) {
    REQUIRE(kept(i) > 0.0);
    middle += eonc::harmonicModeFreeEnergy(
        eonc::tunneling::kHbar * std::sqrt(kept(i)), kTemperature);
  }
  REQUIRE(free[1] == Catch::Approx(middle).margin(1e-10));

  // The same minimum straddling a face of a periodic box still turns
  // freely.
  AtomMatrix across = equilateral(d);
  across.col(0).array() -= 0.5 * d;
  for (long a = 0; a < 3; ++a) {
    if (across(a, 0) < 0.0) {
      across(a, 0) += 10.0;
    }
  }
  std::vector<std::shared_ptr<eonc::Matter>> boxed{
      triangleImage(spring, params, across, true),
      triangleImage(spring, params, across, true)};
  const std::vector<double> inBox = eonc::quantumFreeEnergies(
      boxed, std::vector<std::shared_ptr<AtomMatrix>>(2), kTemperature, params);
  REQUIRE(inBox[0] == Catch::Approx(minimum).margin(1e-6));
}

TEST_CASE("the complement spectrum drops one value per independent direction",
          "[quantum][neb]") {
  Eigen::MatrixXd h = Eigen::MatrixXd::Zero(3, 3);
  h.diagonal() << 1.0, 2.0, 3.0;
  REQUIRE(eonc::complementEigenvalues(h, Eigen::MatrixXd(3, 0)).size() == 3);
  Eigen::MatrixXd removed = Eigen::MatrixXd::Zero(3, 3);
  removed(0, 0) = 1.0;
  removed(0, 1) = 2.0;
  const Eigen::VectorXd two = eonc::complementEigenvalues(h, removed);
  REQUIRE(two.size() == 2);
  REQUIRE(two(0) == Catch::Approx(2.0));
  REQUIRE(two(1) == Catch::Approx(3.0));
  REQUIRE_THROWS_AS(eonc::complementEigenvalues(h, Eigen::MatrixXd::Ones(2, 1)),
                    std::invalid_argument);
}
