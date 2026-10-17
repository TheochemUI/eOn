/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#pragma once

#include "PhysicalQuadraticBath.h"
#include "DerivativeRbfModel.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <memory>
#include <string>

namespace eonc::tunneling {

struct PhysicalBathModel {
  Eigen::VectorXd center;
  Eigen::MatrixXd reactionBasis, bathBasis;
  std::shared_ptr<PhysicalQuadraticBath> bath;
  std::shared_ptr<const DerivativeRbfModel> potentialModel;
  std::string qualification;
  double stationaryTolerance=0;

  static Eigen::VectorXd vector(const nlohmann::json &value) {
    if (!value.is_array()) throw std::invalid_argument("physical bath model vector is not an array");
    Eigen::VectorXd out(value.size());
    for (size_t i=0;i<value.size();++i) out[i]=value.at(i).get<double>();
    if (!out.allFinite()) throw std::invalid_argument("physical bath model vector is not finite");
    return out;
  }
  static Eigen::MatrixXd matrix(const nlohmann::json &value,long columns=-1) {
    if (!value.is_array()) throw std::invalid_argument("physical bath model matrix is not an array");
    const long cols=value.empty() ? std::max(0L,columns) : static_cast<long>(value.at(0).size());
    Eigen::MatrixXd out(value.size(),cols);
    for (size_t i=0;i<value.size();++i) {
      if (!value.at(i).is_array() || value.at(i).size()!=static_cast<size_t>(cols))
        throw std::invalid_argument("physical bath model matrix has inconsistent rows");
      for (long j=0;j<cols;++j) out(i,j)=value.at(i).at(j).get<double>();
    }
    if (!out.allFinite()) throw std::invalid_argument("physical bath model matrix is not finite");
    return out;
  }
  static PhysicalBathModel read(const std::string &filename,long dimension) {
    std::ifstream stream(filename);
    if (!stream) throw std::runtime_error("cannot read physical bath model "+filename);
    nlohmann::json data; stream>>data;
    if (data.at("format")!="eon.physical_quadratic_bath.v1" || data.at("callback")!="frozen")
      throw std::invalid_argument("physical bath model requires v1 and frozen callback convention");
    PhysicalBathModel out;
    out.center=vector(data.at("center"));
    out.reactionBasis=matrix(data.at("reaction_basis"));
    out.bathBasis=matrix(data.at("bath_basis"));
    out.qualification=data.at("qualification").get<std::string>();
    out.stationaryTolerance=data.at("stationary_tolerance").get<double>();
    out.bath=std::make_shared<PhysicalQuadraticBath>();
    out.bath->coupling=matrix(data.at("coupling"));
    out.bath->curvatures=vector(data.at("curvatures"));
    const long r=out.reactionBasis.cols(),b=out.bathBasis.cols();
    if (out.center.size()!=dimension || out.reactionBasis.rows()!=dimension ||
        out.bathBasis.rows()!=dimension || r+b!=dimension || r<1 ||
        out.bath->curvatures.size()!=b || out.qualification.empty() ||
        !(out.stationaryTolerance>0) || !std::isfinite(out.stationaryTolerance))
      throw std::invalid_argument("physical bath model dimension or qualification mismatch");
    Eigen::MatrixXd basis(dimension,dimension); basis<<out.reactionBasis,out.bathBasis;
    if (!(basis.transpose()*basis).isApprox(Eigen::MatrixXd::Identity(dimension,dimension),1e-10))
      throw std::invalid_argument("physical bath model basis is not complete and orthonormal");
    // Reactant coordinates define q=0 in the physical mass-weighted input.
    if (out.center.norm()>1e-10)
      throw std::invalid_argument("physical bath model center must match the input reactant");
    out.bath->validate(r);
    if (data.contains("curvature_model")) {
      const auto &spec = data.at("curvature_model");
      if (spec.at("kind") != "convex_endpoint_blend")
        throw std::invalid_argument("unsupported physical bath curvature model");
      const Eigen::MatrixXd endpoint = matrix(spec.at("endpoint"));
      const Eigen::VectorXd coordinate = vector(spec.at("coordinate"));
      if (endpoint.rows() != b || endpoint.cols() != b ||
          coordinate.size() != r ||
          !endpoint.isApprox(endpoint.transpose(), 1e-12))
        throw std::invalid_argument("physical bath curvature model dimension or symmetry mismatch");
      Eigen::LLT<Eigen::MatrixXd> positive(endpoint);
      if (positive.info() != Eigen::Success)
        throw std::invalid_argument("physical bath curvature endpoint is not positive definite");
      const Eigen::MatrixXd initial = out.bath->curvatures.asDiagonal();
      const Eigen::MatrixXd change = endpoint - initial;
      out.bath->curvature = [initial, change, coordinate](const Eigen::VectorXd &q, int order) {
        PhysicalQuadraticBath::CurvatureDerivatives result;
        const double p = coordinate.dot(q);
        double weight = 0, first = 0, second = 0;
        if (p >= 1) weight = 1;
        else if (p > 0) {
          weight = p*p*p*(10+p*(-15+6*p));
          first = 30*p*p*(1-p)*(1-p);
          second = 60*p*(1-p)*(1-2*p);
        }
        result.value = initial + weight * change;
        if (order >= 1) {
          result.first.resize(coordinate.size());
          for (long i = 0; i < coordinate.size(); ++i)
            result.first[i] = first * coordinate[i] * change;
        }
        if (order >= 2) {
          result.second.resize(coordinate.size());
          for (long i = 0; i < coordinate.size(); ++i) {
            result.second[i].resize(coordinate.size());
            for (long j = 0; j < coordinate.size(); ++j)
              result.second[i][j] = second * coordinate[i] * coordinate[j] * change;
          }
        }
        return result;
      };
    }
    if (data.contains("source_model")) {
      auto model=std::make_shared<const DerivativeRbfModel>(DerivativeRbfModel::fromJson(data.at("source_model")));
      if (model->dimension()!=r || model->outputs()!=b)
        throw std::invalid_argument("physical bath source interpolant dimension mismatch");
      out.bath->source=[model](const Eigen::VectorXd &q,Eigen::VectorXd &s,Eigen::MatrixXd &j,std::vector<Eigen::MatrixXd> &h) {
        const auto value=model->evaluate(q);
        s=value.value;j=value.jacobian;h=value.hessians;
      };
    }
    if (data.contains("potential_model")) {
      out.potentialModel=std::make_shared<const DerivativeRbfModel>(DerivativeRbfModel::fromJson(data.at("potential_model")));
      if (out.potentialModel->dimension()!=r || out.potentialModel->outputs()!=1)
        throw std::invalid_argument("physical bath potential interpolant dimension mismatch");
    }
    return out;
  }

  Eigen::VectorXd project(const Eigen::VectorXd &q) const {
    return reactionBasis.transpose()*(q-center);
  }
  Eigen::VectorXd frozen(const Eigen::VectorXd &q) const {
    return center+reactionBasis*q;
  }
  std::vector<Eigen::VectorXd> reconstruct(const std::vector<Eigen::VectorXd> &q,double spring) const {
    const auto eliminated=bath->evaluate(q,spring);
    std::vector<Eigen::VectorXd> out(q.size());
    for (size_t j=0;j<q.size();++j)
      out[j]=frozen(q[j])+bathBasis*eliminated.bathPositions[j];
    return out;
  }
};
} // namespace eonc::tunneling
