/*
** This file is part of eOn.
** SPDX-License-Identifier: BSD-3-Clause
** Copyright (c) 2010--present, eOn Development Team
*/
#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace eonc::tunneling {

/// Exact factorization of diag(E_j) plus the cyclic primitive spring matrix.
/// Every physical bath block is positive definite; no spectral floor is used.
class CyclicBlockBath {
public:
  CyclicBlockBath(const std::vector<Eigen::MatrixXd> &blocks, double spring)
      : count_(blocks.size()),
        width_(blocks.empty() ? 0 : blocks.front().rows()) {
    if (count_ == 0 || width_ == 0 || !std::isfinite(spring) || spring < 0)
      throw std::invalid_argument("invalid cyclic bath dimensions or spring");
    for (const auto &block : blocks) {
      if (block.rows() != width_ || block.cols() != width_ ||
          !block.allFinite() ||
          (block - block.transpose()).norm() >
              64 * std::numeric_limits<double>::epsilon() *
                  std::max(1.0, block.norm()))
        throw std::invalid_argument(
            "invalid or nonsymmetric physical bath block");
      Eigen::LLT<Eigen::MatrixXd> physical(block);
      if (physical.info() != Eigen::Success)
        throw std::invalid_argument(
            "physical bath block is not positive definite");
    }
    factors_.resize(count_);
    next_.resize(count_ - 1);
    boundary_.resize(count_ - 1);
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(width_, width_);
    std::vector<Eigen::MatrixXd> diagonal = blocks;
    if (count_ > 1)
      for (auto &block : diagonal)
        block += 2 * spring * identity;
    Eigen::MatrixXd coupling = -(count_ == 2 ? 2 : 1) * spring * identity;
    for (std::size_t j = 0; j + 1 < count_; ++j) {
      factor(j, diagonal[j]);
      boundary_[j] = factors_[j].solve(coupling);
      diagonal.back() -= coupling.transpose() * boundary_[j];
      if (j + 2 < count_) {
        const Eigen::MatrixXd adjacent = -spring * identity;
        next_[j] = factors_[j].solve(adjacent);
        diagonal[j + 1] -= adjacent.transpose() * next_[j];
        Eigen::MatrixXd direct = j + 3 == count_
                                     ? Eigen::MatrixXd(-spring * identity)
                                     : Eigen::MatrixXd::Zero(width_, width_);
        coupling = direct - adjacent.transpose() * boundary_[j];
      } else {
        next_[j] = Eigen::MatrixXd::Zero(width_, width_);
      }
    }
    factor(count_ - 1, diagonal.back());
  }

  Eigen::MatrixXd solve(const Eigen::MatrixXd &rhs) const {
    if (rhs.rows() != static_cast<Eigen::Index>(count_) * width_ ||
        !rhs.allFinite())
      throw std::invalid_argument("invalid cyclic bath right hand side");
    Eigen::MatrixXd result = rhs;
    const Eigen::Index last = static_cast<Eigen::Index>(count_ - 1) * width_;
    for (std::size_t j = 0; j + 1 < count_; ++j) {
      const Eigen::Index row = static_cast<Eigen::Index>(j) * width_;
      const Eigen::MatrixXd current = result.middleRows(row, width_);
      if (j + 2 < count_)
        result.middleRows(row + width_, width_) -=
            next_[j].transpose() * current;
      result.middleRows(last, width_) -= boundary_[j].transpose() * current;
      result.middleRows(row, width_) = factors_[j].solve(current);
    }
    result.middleRows(last, width_) =
        factors_.back().solve(result.middleRows(last, width_).eval());
    for (std::size_t j = count_ - 1; j-- > 0;) {
      const Eigen::Index row = static_cast<Eigen::Index>(j) * width_;
      result.middleRows(row, width_) -=
          boundary_[j] * result.middleRows(last, width_);
      if (j + 2 < count_)
        result.middleRows(row, width_) -=
            next_[j] * result.middleRows(row + width_, width_);
    }
    return result;
  }

  double logDeterminant() const { return logDeterminant_; }

private:
  void factor(std::size_t index, const Eigen::MatrixXd &matrix) {
    factors_[index].compute((matrix + matrix.transpose()) * 0.5);
    if (factors_[index].info() != Eigen::Success)
      throw std::runtime_error(
          "cyclic bath Schur block is not positive definite");
    const Eigen::MatrixXd lower = factors_[index].matrixL();
    logDeterminant_ += 2 * lower.diagonal().array().log().sum();
    if (!std::isfinite(logDeterminant_))
      throw std::runtime_error("nonfinite cyclic bath determinant");
  }

  std::size_t count_;
  Eigen::Index width_;
  std::vector<Eigen::LLT<Eigen::MatrixXd>> factors_;
  std::vector<Eigen::MatrixXd> next_, boundary_;
  double logDeterminant_ = 0;
};

} // namespace eonc::tunneling
