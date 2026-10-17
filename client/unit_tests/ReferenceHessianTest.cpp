/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"
#include "eon/ReferenceHessian.h"

#include <filesystem>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <limits>

namespace tests {
static eonc::helpers::test::QuillTestLogger logger;

TEST_CASE("reference Hessians retain signed curvature and bind the complete geometry",
          "[instanton][hessian][reference]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::sharePotential(eonc::helpers::makePotential(PotType::LJ, params));
  Matter reference(pot, params);
  REQUIRE(eonc::io::io_ok(reference.con2matter(std::string("reactant.con"))));
  for (long i = 0; i < reference.numberOfAtoms(); ++i) {
    reference.setFixed(i, i == 0 ? 0 : 1);
  }
  const auto directory = std::filesystem::temp_directory_path() /
                         ("eon-reference-hessian-" + std::to_string(
                             std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directories(directory);
  struct Cleanup {
    std::filesystem::path path;
    ~Cleanup() { std::error_code error; std::filesystem::remove_all(path, error); }
  } cleanup{directory};
  const std::string filename = (directory / "reference.dat").string();
  Matter recorded(reference);
  MatrixXd matrix(3, 3);
  matrix << 2.0, 0.5, 0.0, 0.5, -1.0, 0.1, 0.0, 0.1, 3.0;
  bool reject = false;
  SECTION("matching signed matrix") {}
  SECTION("mobile position") {
    auto positions = recorded.getPositions();
    positions(0, 0) += 0.02;
    recorded.setPositions(positions);
    reject = true;
  }
  SECTION("frozen position") {
    auto positions = recorded.getPositions();
    positions(1, 0) += 0.02;
    recorded.setPositions(positions);
    reject = true;
  }
  SECTION("mass metric") {
    for (long i = 0; i < recorded.numberOfAtoms(); ++i) {
      recorded.setMass(i, 2.0 * reference.getMass(i));
    }
    reject = true;
  }
  SECTION("constraints") {
    recorded.setFixed(1, 0);
    reject = true;
  }
  SECTION("cell") {
    auto cell = recorded.getCell();
    cell(0, 0) += 0.1;
    recorded.setCell(cell);
    reject = true;
  }
  SECTION("asymmetric matrix") {
    matrix(0, 1) += 0.01;
    reject = true;
  }
  REQUIRE(eonc::io::io_ok(recorded.matter2con(filename + ".con", false)));
  {
    std::ofstream output(filename);
    output << std::setprecision(std::numeric_limits<double>::max_digits10) << matrix << '\n';
  }
  if (reject) {
    REQUIRE_THROWS(eonc::tunneling::readReferenceHessian(filename, reference, 3));
  } else {
    const auto restored = eonc::tunneling::readReferenceHessian(filename, reference, 3);
    REQUIRE(restored == matrix);
    REQUIRE(restored(1, 1) == -1.0);
    REQUIRE_THROWS(eonc::tunneling::readReferenceHessian(filename, reference, 6));
    {
      std::ofstream output(filename, std::ios::app);
      output << " 4.0\n";
    }
    REQUIRE_THROWS(eonc::tunneling::readReferenceHessian(filename, reference, 3));
    {
      std::ofstream output(filename);
      output << "1 2\n";
    }
    REQUIRE_THROWS(eonc::tunneling::readReferenceHessian(filename, reference, 3));
  }
}
} // namespace tests
