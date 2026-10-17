/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#pragma once
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace eonc::tunneling {

struct RbfEvaluation {
  Eigen::VectorXd value;
  Eigen::MatrixXd jacobian;
  std::vector<Eigen::MatrixXd> hessians;
};

/// Gaussian derivative interpolant with an explicit quadratic baseline.
/// Orders differentiate the site's argument; returned derivatives act on q.
class DerivativeRbfModel {
public:
  static DerivativeRbfModel fromJson(const nlohmann::json &data) {
    if (data.at("kind").get<std::string>() != "derivative_rbf")
      throw std::invalid_argument("unsupported reduced interpolant kind");
    DerivativeRbfModel model;
    model.scales_ = vector(data.at("length_scales"));
    model.constant_ = vector(data.at("constant"));
    const Eigen::Index d=model.scales_.size(), outputs=model.constant_.size();
    if (d<1 || outputs<1 || !(model.scales_.array()>0).all())
      throw std::invalid_argument("interpolant dimensions and scales must be positive");
    model.sites_=matrix(data.at("sites"),d);
    model.weights_=matrix(data.at("weights"),outputs);
    model.linear_=data.contains("linear") ? matrix(data.at("linear"),d)
                                          : Eigen::MatrixXd::Zero(outputs,d);
    if (model.weights_.rows()!=model.sites_.rows() || model.linear_.rows()!=outputs)
      throw std::invalid_argument("interpolant coefficient dimensions differ");
    const auto &orders=data.at("derivative_orders");
    if (orders.size()!=static_cast<size_t>(model.sites_.rows()))
      throw std::invalid_argument("interpolant derivative count differs");
    for (const auto &row:orders) {
      if (row.size()!=static_cast<size_t>(d))
        throw std::invalid_argument("interpolant derivative dimension differs");
      std::vector<int> values;
      int total=0;
      for (const auto &entry:row) {
        const int order=entry.get<int>();
        if (order<0 || order>2)
          throw std::invalid_argument("interpolant site derivative order must be 0, 1 or 2");
        values.push_back(order); total+=order;
      }
      if (total>2) throw std::invalid_argument("interpolant site derivative total exceeds two");
      model.orders_.push_back(std::move(values));
    }
    model.quadratic_.assign(static_cast<size_t>(outputs),Eigen::MatrixXd::Zero(d,d));
    if (data.contains("quadratic")) {
      if (data.at("quadratic").size()!=static_cast<size_t>(outputs))
        throw std::invalid_argument("interpolant quadratic output dimension differs");
      for (Eigen::Index i=0;i<outputs;++i) {
        auto h=matrix(data.at("quadratic").at(i),d);
        if (h.rows()!=d || !h.isApprox(h.transpose(),1e-12))
          throw std::invalid_argument("interpolant quadratic baseline must be symmetric");
        model.quadratic_[static_cast<size_t>(i)]=std::move(h);
      }
    }
    return model;
  }

  Eigen::Index dimension() const { return scales_.size(); }
  Eigen::Index outputs() const { return constant_.size(); }

  RbfEvaluation evaluate(const Eigen::VectorXd &q) const {
    if (q.size()!=dimension() || !q.allFinite())
      throw std::invalid_argument("interpolant coordinate dimension or finiteness differs");
    RbfEvaluation result{constant_+linear_*q,linear_,quadratic_};
    for (Eigen::Index k=0;k<outputs();++k) {
      result.value[k]+=.5*q.dot(quadratic_[static_cast<size_t>(k)]*q);
      result.jacobian.row(k)+=(quadratic_[static_cast<size_t>(k)]*q).transpose();
    }
    std::vector<int> derivative(static_cast<size_t>(dimension()),0);
    for (Eigen::Index site=0;site<sites_.rows();++site) {
      result.value.noalias()+=kernel(q,site,derivative)*weights_.row(site).transpose();
      for (Eigen::Index i=0;i<dimension();++i) {
        derivative[static_cast<size_t>(i)]++;
        result.jacobian.col(i).noalias()+=kernel(q,site,derivative)*weights_.row(site).transpose();
        for (Eigen::Index j=0;j<dimension();++j) {
          derivative[static_cast<size_t>(j)]++;
          const double value=kernel(q,site,derivative);
          for (Eigen::Index output=0;output<outputs();++output)
            result.hessians[static_cast<size_t>(output)](i,j)+=value*weights_(site,output);
          derivative[static_cast<size_t>(j)]--;
        }
        derivative[static_cast<size_t>(i)]--;
      }
    }
    return result;
  }

private:
  static Eigen::VectorXd vector(const nlohmann::json &data) {
    if (!data.is_array()) throw std::invalid_argument("interpolant vector must be an array");
    Eigen::VectorXd out(data.size());
    for (Eigen::Index i=0;i<out.size();++i) out[i]=data.at(i).get<double>();
    if (!out.allFinite()) throw std::invalid_argument("non-finite interpolant vector");
    return out;
  }
  static Eigen::MatrixXd matrix(const nlohmann::json &data,Eigen::Index columns) {
    if (!data.is_array()) throw std::invalid_argument("interpolant matrix must be an array");
    Eigen::MatrixXd out(data.size(),columns);
    for (Eigen::Index i=0;i<out.rows();++i) {
      if (data.at(i).size()!=static_cast<size_t>(columns))
        throw std::invalid_argument("interpolant matrix column count differs");
      for (Eigen::Index j=0;j<columns;++j) out(i,j)=data.at(i).at(j).get<double>();
    }
    if (!out.allFinite()) throw std::invalid_argument("non-finite interpolant matrix");
    return out;
  }
  static double hermite(int n,double z) {
    if (n==0) return 1.;
    if (n==1) return z;
    double previous=1.,current=z;
    for (int k=1;k<n;++k) {
      const double next=z*current-k*previous;
      previous=current; current=next;
    }
    return current;
  }
  double kernel(const Eigen::VectorXd &q,Eigen::Index site,
                const std::vector<int> &inputOrder) const {
    const Eigen::VectorXd z=(q-sites_.row(site).transpose()).array()/scales_.array();
    double value=std::exp(-.5*z.squaredNorm());
    for (Eigen::Index j=0;j<dimension();++j) {
      const int input=inputOrder[static_cast<size_t>(j)];
      const int total=input+orders_[static_cast<size_t>(site)][static_cast<size_t>(j)];
      value*=hermite(total,z[j])/std::pow(scales_[j],total);
      if (input%2) value=-value;
    }
    return value;
  }
  Eigen::VectorXd scales_,constant_;
  Eigen::MatrixXd sites_,weights_,linear_;
  std::vector<std::vector<int>> orders_;
  std::vector<Eigen::MatrixXd> quadratic_;
};
} // namespace eonc::tunneling
