/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#include "eon/CyclicBlockBath.h"
#include "catch2/catch_amalgamated.hpp"

using eonc::tunneling::CyclicBlockBath;

TEST_CASE("cyclic block bath solves the full noncommuting physical matrix",
          "[instanton][physical-bath][cyclic]") {
  for (int count : {1, 2, 3, 5}) {
    for (double spring : {0.0, 0.17, 3.2}) {
      std::vector<Eigen::MatrixXd> blocks;
      Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(2 * count, 2 * count);
      for (int j = 0; j < count; ++j) {
        Eigen::Matrix2d transform;
        transform << 1 + 0.1 * j, 0.2 + 0.13 * j, -0.1 * j, 0.8;
        blocks.push_back(transform.transpose() * transform +
                         0.3 * Eigen::Matrix2d::Identity());
        dense.block<2, 2>(2 * j, 2 * j) += blocks.back();
      }
      for (int j = 0; j < count; ++j) {
        const int k = (j + 1) % count;
        dense.block<2, 2>(2 * j, 2 * j) += spring * Eigen::Matrix2d::Identity();
        dense.block<2, 2>(2 * k, 2 * k) += spring * Eigen::Matrix2d::Identity();
        dense.block<2, 2>(2 * j, 2 * k) -= spring * Eigen::Matrix2d::Identity();
        dense.block<2, 2>(2 * k, 2 * j) -= spring * Eigen::Matrix2d::Identity();
      }
      Eigen::MatrixXd rhs(2 * count, 3);
      for (Eigen::Index i = 0; i < rhs.rows(); ++i)
        for (Eigen::Index k = 0; k < rhs.cols(); ++k)
          rhs(i, k) = std::sin(0.3 + i + 0.2 * k);
      const CyclicBlockBath factor(blocks, spring);
      const Eigen::MatrixXd solution = factor.solve(rhs);
      REQUIRE((dense * solution - rhs).norm() < 2e-13);
      REQUIRE((solution - dense.llt().solve(rhs)).norm() < 2e-13);
      Eigen::LLT<Eigen::MatrixXd> reference(dense);
      const Eigen::MatrixXd lower = reference.matrixL();
      REQUIRE(std::abs(factor.logDeterminant() -
                       2 * lower.diagonal().array().log().sum()) < 2e-13);
      REQUIRE((factor.solve(Eigen::MatrixXd::Identity(2 * count, 2 * count)) -
               dense.inverse())
                  .norm() < 2e-13);
    }
  }
}

TEST_CASE("cyclic block bath rejects invalid physical blocks without flooring",
          "[instanton][physical-bath][cyclic]") {
  const Eigen::MatrixXd good = Eigen::MatrixXd::Identity(2, 2);
  REQUIRE_THROWS_AS(CyclicBlockBath({}, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(CyclicBlockBath({good}, -1), std::invalid_argument);
  Eigen::MatrixXd bad = good;
  bad(0, 0) = 0;
  REQUIRE_THROWS_AS(CyclicBlockBath({bad}, 1), std::invalid_argument);
  bad(0, 0) = -0.1;
  REQUIRE_THROWS_AS(CyclicBlockBath({good, bad}, 10), std::invalid_argument);
  bad = good;
  bad(0, 1) = 0.2;
  REQUIRE_THROWS_AS(CyclicBlockBath({bad}, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(CyclicBlockBath({good, Eigen::MatrixXd::Identity(3, 3)}, 1),
                    std::invalid_argument);
  const CyclicBlockBath factor({good, good}, 1);
  REQUIRE_THROWS_AS(factor.solve(good), std::invalid_argument);
}
