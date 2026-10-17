/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace eonc::tunneling {

/// V(q,b)=V_fixed(q)+source(q)^T b+b^T diag(curvatures)b/2.
/// The bath is physical and constant, not an Ohmic ring thermostat.
struct PhysicalQuadraticBath {
  Eigen::MatrixXd coupling; ///< retained x bath, source=C^T q by default
  Eigen::VectorXd curvatures; ///< strictly positive physical bath eigenvalues
  /// With V_relaxed=V_fixed-source^T E^-1 source/2, restore its counterterm.
  bool relaxedPotential = false;
  /// Source, bath x retained Jacobian, one retained Hessian per bath mode.
  using Source = std::function<void(const Eigen::VectorXd &, Eigen::VectorXd &,
      Eigen::MatrixXd &, std::vector<Eigen::MatrixXd> &)>;
  Source source;

  struct Evaluation {
    double energy = 0.0;
    std::vector<Eigen::VectorXd> gradient;
    Eigen::MatrixXd hessian;
    std::vector<Eigen::VectorXd> bathPositions;
    std::vector<Eigen::VectorXd> inverseBath;
    std::vector<Eigen::MatrixXd> sourceJacobians;
    double bathLogDeterminant = 0.0;
  };

  void validate(long dof) const {
    if (dof < 1 || coupling.rows()!=dof || coupling.cols()!=curvatures.size() ||
        !coupling.allFinite() || !curvatures.allFinite() ||
        (curvatures.array()<=0).any())
      throw std::invalid_argument("physical bath requires finite coupling and positive bath curvature");
  }

  /// Physical curvature at the conditional bath equilibrium. The negative
  /// eigenvalue of this full operator sets the dynamical crossover frequency.
  Eigen::MatrixXd fullHessian(const Eigen::VectorXd &q,
                             const Eigen::MatrixXd &potentialHessian) const {
    const long r=q.size(),b=curvatures.size();
    validate(r);
    if (potentialHessian.rows()!=r || potentialHessian.cols()!=r ||
        !potentialHessian.allFinite())
      throw std::invalid_argument("physical bath potential Hessian mismatch");
    Eigen::VectorXd s;
    Eigen::MatrixXd jac;
    std::vector<Eigen::MatrixXd> second;
    if (source) source(q,s,jac,second);
    else { s=coupling.transpose()*q;jac=coupling.transpose();second.assign(b,Eigen::MatrixXd::Zero(r,r)); }
    if (s.size()!=b || jac.rows()!=b || jac.cols()!=r || !s.allFinite() ||
        !jac.allFinite() || second.size()!=static_cast<size_t>(b))
      throw std::invalid_argument("physical bath source derivative mismatch");
    for (const auto &h:second)
      if (h.rows()!=r || h.cols()!=r || !h.allFinite() || !h.isApprox(h.transpose(),1e-12))
        throw std::invalid_argument("physical bath source Hessian mismatch");
    Eigen::MatrixXd out=Eigen::MatrixXd::Zero(r+b,r+b);
    out.topLeftCorner(r,r)=potentialHessian;
    if (relaxedPotential)
      out.topLeftCorner(r,r).noalias()+=jac.transpose()*curvatures.cwiseInverse().asDiagonal()*jac;
    else for (long k=0;k<b;++k) out.topLeftCorner(r,r).noalias()-=(s[k]/curvatures[k])*second[k];
    out.topRightCorner(r,b)=jac.transpose();
    out.bottomLeftCorner(b,r)=jac;
    out.bottomRightCorner(b,b)=curvatures.asDiagonal();
    if (!out.allFinite() || !out.isApprox(out.transpose(),1e-12))
      throw std::invalid_argument("physical bath full Hessian is not finite and symmetric");
    return out;
  }

  /// Exact finite primitive-ring Gaussian elimination, with c=(beta_N*hbar)^-2.
  Evaluation evaluate(const std::vector<Eigen::VectorXd> &q, double c,
                      bool needHessian = false) const {
    if (q.empty() || !(c>0) || !std::isfinite(c))
      throw std::invalid_argument("physical bath requires a nonempty ring and positive spring");
    const long n=q.size(), r=q.front().size(), b=curvatures.size();
    validate(r);
    Evaluation out;
    out.gradient.assign(n,Eigen::VectorXd::Zero(r));
    out.bathPositions.assign(n,Eigen::VectorXd::Zero(b));
    if (needHessian) out.hessian=Eigen::MatrixXd::Zero(n*r,n*r);
    std::vector<Eigen::VectorXd> s(n);
    std::vector<Eigen::MatrixXd> jac(n);
    std::vector<std::vector<Eigen::MatrixXd>> second(n);
    for (long j=0;j<n;++j) {
      if (q[j].size()!=r || !q[j].allFinite())
        throw std::invalid_argument("physical bath ring coordinate mismatch");
      if (source) source(q[j],s[j],jac[j],second[j]);
      else {
        s[j]=coupling.transpose()*q[j];
        jac[j]=coupling.transpose();
        second[j].assign(b,Eigen::MatrixXd::Zero(r,r));
      }
      if (s[j].size()!=b || jac[j].rows()!=b || jac[j].cols()!=r ||
          !s[j].allFinite() || !jac[j].allFinite() || second[j].size()!=static_cast<size_t>(b))
        throw std::invalid_argument("physical bath source derivative mismatch");
      for (const auto &h:second[j])
        if (h.rows()!=r || h.cols()!=r || !h.allFinite() || !h.isApprox(h.transpose(),1e-12))
          throw std::invalid_argument("physical bath source Hessian mismatch");
    }
    // inverse[lag] is diagonal in the fixed physical bath eigenbasis.
    std::vector<Eigen::VectorXd> inverse(n,Eigen::VectorXd::Zero(b));
    for (long k=0;k<n;++k) {
      const double sk=std::sin(std::numbers::pi*static_cast<double>(k)/n);
      const Eigen::VectorXd lambda=(curvatures.array()+4*c*sk*sk).matrix();
      out.bathLogDeterminant+=lambda.array().log().sum();
      for (long lag=0;lag<n;++lag)
        inverse[lag].array()+=std::cos(2*std::numbers::pi*k*static_cast<double>(lag)/n)/n/lambda.array();
    }
    for (long j=0;j<n;++j) {
      Eigen::VectorXd y=Eigen::VectorXd::Zero(b);
      for (long l=0;l<n;++l) y.array()+=inverse[(j-l+n)%n].array()*s[l].array();
      out.bathPositions[j]=-y;
      if (relaxedPotential) y.array()-=s[j].array()/curvatures.array();
      out.energy-=0.5*s[j].dot(y);
      out.gradient[j].noalias()-=jac[j].transpose()*y;
      if (needHessian) {
        for (long l=0;l<n;++l) {
          Eigen::VectorXd weights=inverse[(j-l+n)%n];
          if (relaxedPotential && j==l) weights.array()-=curvatures.cwiseInverse().array();
          out.hessian.block(j*r,l*r,r,r).noalias()-=jac[j].transpose()*weights.asDiagonal()*jac[l];
        }
        for (long a=0;a<b;++a) out.hessian.block(j*r,j*r,r,r).noalias()-=y[a]*second[j][a];
      }
    }
    if (needHessian) {
      out.inverseBath=inverse;
      out.sourceJacobians=jac;
    }
    return out;
  }
};
} // namespace eonc::tunneling
