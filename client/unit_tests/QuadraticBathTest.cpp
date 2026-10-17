/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#include "eon/QuadraticBath.h"
#include "catch2/catch_amalgamated.hpp"

using eonc::tunneling::QuadraticBath;

TEST_CASE("physical bath elimination retains frequency response and determinant",
          "[instanton][quadratic-bath]") {
  Eigen::Matrix3d h;
  h << -1.0,0.4,0.2, 0.4,2.0,0.3, 0.2,0.3,4.0;
  const Eigen::Vector3d n=Eigen::Vector3d(1.0,0.2,-0.1).normalized();
  QuadraticBath bath(h,n);
  REQUIRE(bath.reactionBasis().cols()==1);
  REQUIRE(bath.bathCurvatures().size()==2);
  REQUIRE(bath.bathCurvatures().minCoeff()>0);
  for (double nu : {0.0,0.2,1.7,4.0}) {
    INFO(nu);
    const Eigen::Matrix3d full=h+nu*nu*Eigen::Matrix3d::Identity();
    const Eigen::MatrixXd reduced=bath.kernel(nu);
    REQUIRE(std::abs(full.determinant()-std::exp(bath.logBathDeterminant(nu))*reduced.determinant())<1e-10);
    const Eigen::MatrixXd projected=bath.reactionBasis().transpose()*full.inverse()*bath.reactionBasis();
    REQUIRE(projected.isApprox(reduced.inverse(),1e-12));
    const Eigen::MatrixXd memoryForm=bath.kernel(0)+nu*nu*Eigen::MatrixXd::Identity(1,1)+nu*bath.memory(nu);
    REQUIRE(reduced.isApprox(memoryForm,1e-13));
    REQUIRE(bath.memory(nu)(0,0)>=0);
  }
  const double betaHbar=2.3;
  const double expected=-(2*(betaHbar*bath.bathCurvatures().array().sqrt()/2).sinh()).log().sum();
  REQUIRE(std::abs(bath.logBathPartition(betaHbar)-expected)<1e-14);
}

TEST_CASE("physical bath promotes unstable and zero complementary directions",
          "[instanton][quadratic-bath]") {
  Eigen::Matrix4d h=Eigen::Vector4d(2.0,-3.0,0.0,5.0).asDiagonal();
  QuadraticBath bath(h,Eigen::Vector4d::UnitX());
  REQUIRE(bath.reactionBasis().cols()==3);
  REQUIRE(bath.bathCurvatures().size()==1);
  REQUIRE(std::abs(bath.bathCurvatures()[0]-5.0)<1e-14);
  REQUIRE(std::abs(bath.kernel(0).determinant())<1e-14);
  const Eigen::Matrix4d orthogonal=(Eigen::Matrix4d() << bath.reactionBasis(),bath.bathBasis()).finished();
  REQUIRE((orthogonal.transpose()*orthogonal).isApprox(Eigen::Matrix4d::Identity(),1e-13));
  QuadraticBath all(-Eigen::Matrix3d::Identity(),Eigen::Vector3d::UnitX());
  REQUIRE(all.reactionBasis().cols()==3);
  REQUIRE(all.bathCurvatures().size()==0);
  REQUIRE(all.logBathPartition(1.0)==0);
  REQUIRE(all.logBathDeterminant(1.0)==0);
  REQUIRE(all.kernel(0).isApprox(-Eigen::Matrix3d::Identity(),1e-14));
}

TEST_CASE("physical bath preserves rotated multi-coordinate reaction subspaces",
          "[instanton][quadratic-bath]") {
  Eigen::Matrix3d h;
  h << -1.0,0.4,0.2, 0.4,2.0,0.3, 0.2,0.3,4.0;
  Eigen::Matrix<double,3,2> r;
  r << 1,0,0,1,0,0;
  QuadraticBath bath(h,r);
  REQUIRE(bath.reactionBasis().cols()==2);
  const Eigen::MatrixXd expected=bath.reactionCurvature()
      -bath.coupling()*bath.bathCurvatures().cwiseInverse().asDiagonal()*bath.coupling().transpose();
  REQUIRE(bath.kernel(0).isApprox(expected,1e-14));
  REQUIRE(std::abs((h+0.25*Eigen::Matrix3d::Identity()).determinant()
      -std::exp(bath.logBathDeterminant(0.5))*bath.kernel(0.5).determinant())<1e-13);
  Eigen::Matrix3d rotation=Eigen::AngleAxisd(0.7,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
  QuadraticBath rotated(rotation*h*rotation.transpose(),rotation*r);
  REQUIRE(std::abs(bath.logBathDeterminant(0.5)-rotated.logBathDeterminant(0.5))<1e-13);
  REQUIRE(std::abs(bath.kernel(0.5).determinant()-rotated.kernel(0.5).determinant())<1e-13);
}

TEST_CASE("physical bath rejects malformed coordinates without changing curvature",
          "[instanton][quadratic-bath]") {
  Eigen::Matrix2d h=Eigen::Matrix2d::Identity();
  REQUIRE_THROWS_AS(QuadraticBath(h,Eigen::Vector2d(2,0)),std::invalid_argument);
  h(0,1)=0.1;
  REQUIRE_THROWS_AS(QuadraticBath(h,Eigen::Vector2d::UnitX()),std::invalid_argument);
  QuadraticBath bath(Eigen::Matrix2d::Identity(),Eigen::Vector2d::UnitX());
  REQUIRE_THROWS_AS(bath.memory(-1),std::invalid_argument);
  REQUIRE_THROWS_AS(bath.kernel(std::numeric_limits<double>::infinity()),std::invalid_argument);
  REQUIRE_THROWS_AS(bath.logBathPartition(0),std::invalid_argument);
}
