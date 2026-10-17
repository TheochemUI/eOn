/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#include "catch2/catch_amalgamated.hpp"
#include "eon/PhysicalBathModel.h"
#include <chrono>
#include <filesystem>
#include <fstream>

using namespace eonc::tunneling;

namespace {
PhysicalQuadraticBath variableBath() {
  PhysicalQuadraticBath bath;
  bath.coupling = Eigen::MatrixXd::Zero(2, 2);
  bath.curvatures = Eigen::Vector2d(2, 4);
  bath.source = [](const Eigen::VectorXd &q, Eigen::VectorXd &s,
                   Eigen::MatrixXd &j, std::vector<Eigen::MatrixXd> &h) {
    s = Eigen::Vector2d(q[0]*q[0]+0.2*q[1], std::sin(q[1])+0.3*q[0]);
    j.resize(2,2); j << 2*q[0],0.2,0.3,std::cos(q[1]);
    h.assign(2,Eigen::MatrixXd::Zero(2,2));
    h[0](0,0)=2; h[1](1,1)=-std::sin(q[1]);
  };
  bath.curvature = [](const Eigen::VectorXd &q,int order) {
    PhysicalQuadraticBath::CurvatureDerivatives state;
    state.value.resize(2,2);
    state.value << 2+std::exp(q[0]),0.25,0.25,3+std::exp(q[1]);
    if (order>=1) {
      state.first.assign(2,Eigen::MatrixXd::Zero(2,2));
      state.first[0](0,0)=std::exp(q[0]);
      state.first[1](1,1)=std::exp(q[1]);
    }
    if (order>=2) {
      state.second.assign(2,std::vector<Eigen::MatrixXd>(2,Eigen::MatrixXd::Zero(2,2)));
      state.second[0][0](0,0)=std::exp(q[0]);
      state.second[1][1](1,1)=std::exp(q[1]);
    }
    return state;
  };
  return bath;
}
}

TEST_CASE("variable physical bath preserves the full conditional Gaussian",
          "[instanton][physical-bath][variable]") {
  const auto bath=variableBath();
  const double spring=1.7;
  std::vector<Eigen::VectorXd> q={Eigen::Vector2d(.2,-.1),Eigen::Vector2d(-.3,.4),
      Eigen::Vector2d(.15,.25)};
  const auto actual=bath.evaluate(q,spring,true);
  const long n=q.size(),r=2,b=2;
  Eigen::MatrixXd d=Eigen::MatrixXd::Zero(n*b,n*b);
  Eigen::VectorXd s(n*b);
  Eigen::MatrixXd qq=Eigen::MatrixXd::Zero(n*r,n*r);
  Eigen::MatrixXd mixed=Eigen::MatrixXd::Zero(n*b,n*r);
  for (long j=0;j<n;++j) {
    Eigen::VectorXd sj; Eigen::MatrixXd jac;
    std::vector<Eigen::MatrixXd> hs;
    bath.source(q[j],sj,jac,hs);
    s.segment(j*b,b)=sj;
    const auto e=bath.curvature(q[j],2);
    d.block(j*b,j*b,b,b)+=e.value;
    const long k=(j+1)%n;
    d.block(j*b,j*b,b,b)+=spring*Eigen::Matrix2d::Identity();
    d.block(k*b,k*b,b,b)+=spring*Eigen::Matrix2d::Identity();
    d.block(j*b,k*b,b,b)-=spring*Eigen::Matrix2d::Identity();
    d.block(k*b,j*b,b,b)-=spring*Eigen::Matrix2d::Identity();
    const auto &y=actual.bathPositions[j];
    for (long a=0;a<b;++a) qq.block(j*r,j*r,r,r)+=y[a]*hs[a];
    for (long i=0;i<r;++i) {
      mixed.block(j*b,j*r+i,b,1)=jac.col(i)+e.first[i]*y;
      for (long k=0;k<r;++k) qq(j*r+i,j*r+k)+=.5*y.dot(e.second[i][k]*y);
    }
  }
  const Eigen::VectorXd y=-d.llt().solve(s);
  REQUIRE(std::abs(actual.energy+.5*s.dot(d.llt().solve(s)))<1e-13);
  for (long j=0;j<n;++j)
    REQUIRE((actual.bathPositions[j]-y.segment(j*b,b)).norm()<1e-13);
  REQUIRE((actual.hessian-(qq-mixed.transpose()*d.llt().solve(mixed))).norm()<1e-12);
  REQUIRE((actual.solveBath(mixed)-d.llt().solve(mixed)).norm()<1e-13);
  Eigen::LLT<Eigen::MatrixXd> factor(d);
  const Eigen::MatrixXd lower=factor.matrixL();
  REQUIRE(std::abs(actual.bathLogDeterminant-2*lower.diagonal().array().log().sum())<1e-13);
  const double step=2e-5;
  for (long k=0;k<n*r;++k) {
    auto plus=q,minus=q;
    plus[k/r][k%r]+=step; minus[k/r][k%r]-=step;
    const auto p=bath.evaluate(plus,spring),m=bath.evaluate(minus,spring);
    REQUIRE(std::abs((p.energy-m.energy)/(2*step)-actual.gradient[k/r][k%r])<1e-8);
    for (long j=0;j<n;++j)
      REQUIRE(((p.gradient[j]-m.gradient[j])/(2*step)-
               actual.hessian.block(j*r,k,r,1)).norm()<2e-8);
  }
  const auto single=bath.evaluate({q[0]},spring,true);
  const auto full=bath.fullHessian(q[0],Eigen::Matrix2d::Identity());
  REQUIRE((full.topLeftCorner(2,2)-full.topRightCorner(2,2)*
           full.bottomRightCorner(2,2).llt().solve(full.bottomLeftCorner(2,2))-
           Eigen::Matrix2d::Identity()-single.hessian).norm()<1e-12);
  auto invalid=bath;
  invalid.curvature=[](const Eigen::VectorXd&,int) {
    PhysicalQuadraticBath::CurvatureDerivatives state;
    state.value=-Eigen::Matrix2d::Identity(); return state;
  };
  REQUIRE_THROWS_AS(invalid.evaluate(q,spring),std::invalid_argument);
}

TEST_CASE("convex bath model preserves endpoint curvature jets in one basis",
          "[instanton][physical-bath][variable][input]") {
  const auto path=std::filesystem::temp_directory_path()/
      ("eon-bath-curvature-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".json");
  struct Cleanup { std::filesystem::path path; ~Cleanup(){std::filesystem::remove(path);} } cleanup{path};
  nlohmann::json data={{"format","eon.physical_quadratic_bath.v1"},{"callback","frozen"},
      {"center",{0,0,0,0}},{"reaction_basis",{{1,0},{0,1},{0,0},{0,0}}},
      {"bath_basis",{{0,0},{0,0},{1,0},{0,1}}},{"coupling",{{0,0},{0,0}}},
      {"curvatures",{2,4}},{"qualification","quadratic input regression"},
      {"stationary_tolerance",1e-8},{"curvature_model",{{"kind","convex_endpoint_blend"},
      {"coordinate",{1,0}},{"endpoint",{{3,.4},{.4,2.5}}}}}};
  auto read=[&]() {std::ofstream stream(path); stream<<data; stream.close();
    return PhysicalBathModel::read(path.string(),4);};
  const auto model=read();
  for (double p:{-.2,0.,1.,1.2}) {
    const auto state=model.bath->curvature(Eigen::Vector2d(p,.3),2);
    REQUIRE(state.first[0].norm()==0);
    REQUIRE(state.second[0][0].norm()==0);
    Eigen::Matrix2d expected;
    if (p<=0) expected=Eigen::Vector2d(2,4).asDiagonal();
    else expected<<3,.4,.4,2.5;
    REQUIRE((state.value-expected).norm()<1e-14);
  }
  const double p=.37,step=1e-5;
  const auto state=model.bath->curvature(Eigen::Vector2d(p,0),2);
  const auto plus=model.bath->curvature(Eigen::Vector2d(p+step,0),1);
  const auto minus=model.bath->curvature(Eigen::Vector2d(p-step,0),1);
  REQUIRE(((plus.value-minus.value)/(2*step)-state.first[0]).norm()<2e-8);
  REQUIRE(((plus.first[0]-minus.first[0])/(2*step)-state.second[0][0]).norm()<2e-8);
  data["curvature_model"]["endpoint"]={{-1,0},{0,1}};
  REQUIRE_THROWS_AS(read(),std::invalid_argument);
}
