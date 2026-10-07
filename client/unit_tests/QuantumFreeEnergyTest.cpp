#include "catch2/catch_amalgamated.hpp"
#include "eon/Hessian.h"
#include "eon/Potential.h"
#include "eon/QuantumFreeEnergy.h"
#include "eon/Tunneling.h"

#include <Eigen/Eigenvalues>
#include <cmath>

namespace {

struct CurvedValley final : eonc::Potential {
  double ka{2.0};
  double kp{8.0};
  double kz{3.0};
  double alpha{0.5};

  CurvedValley() : Potential(eonc::PotType::LJ) {}

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
  const Eigen::Matrix3d hessian =
      analyticHessian(kX, kAlpha, kA, kP, kZ);
  const Eigen::Vector3d tangent(1.0, 2.0 * kAlpha * kX, 0.0);
  constexpr double kTemperature = 300.0;
  const double summed = directSum(hessian, tangent, kTemperature);
  const double profile =
      eonc::perpendicularHarmonicFreeEnergy(hessian, tangent, kTemperature);
  REQUIRE(profile == Catch::Approx(summed).margin(1e-10));

  const double zpe = eonc::perpendicularHarmonicFreeEnergy(hessian, tangent, 0.0);
  REQUIRE(zpe == Catch::Approx(directSum(hessian, tangent, 0.0)).margin(1e-12));
  constexpr double kValleyEnergy = 0.5 * kA * kX * kX;
  REQUIRE(kValleyEnergy + zpe ==
          Catch::Approx(kValleyEnergy + directSum(hessian, tangent, 0.0))
              .margin(1e-12));

  Parameters params;
  ParametersLoadAccess::hessian_options(params).fd_scheme = "central";
  ParametersLoadAccess::main_options(params).finiteDifference = 1e-5;
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
  REQUIRE(eonc::perpendicularHarmonicFreeEnergy(weighted, tangent, kTemperature) ==
          Catch::Approx(summed).margin(1e-5));
}
