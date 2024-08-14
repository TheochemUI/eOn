/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once
#include "Eigen.h"
#include "Optimizer.h"
#include <spdlog/sinks/basic_file_sink.h>
namespace eonc {
/**
 * @file
 * @ingroup Optimizers
 *
 * \brief Direct optimization for energy minimization
 *
 * The conjugate gradient method is an algorithm for the numerical solution of
 * particular stystems of linear equations, namely symetric,positive-definite
 * ones.
 *
 */

/**
 * Decleration of the Conjugate Gradients optimizer
 */

class ConjugateGradients : public Optimizer<ConjugateGradients> {
public:
  struct Params final {
    size_t max_iter_before_reset{0};
    bool no_overshooting{false};
    bool knock_out_max_move{false};
    ScalarType line_convergence{0.1};
    bool line_search{false};
    size_t max_line_serach_iter{10};
    ScalarType finite_diff{0.01};
    bool saddle_bowl_breakout{false};
  };

public:
  //! Conjugate Gradients optimizer constructor
  /*!
   * \param std::shared_ptr<ObjectiveFunction> m_objf that tells the optimizer
   * how to run \param std::shared_ptr<Parameters> m_params defined by the
   * config.init file
   */
  ConjugateGradients(const ObjectiveFunction &a_objf,
                     const ConjugateGradients::Params &p_a)
      : Optimizer(a_objf),
        m_p{p_a},
        m_directionOld{(a_objf.getPositions()).setZero()},
        m_forceOld{(a_objf.getPositions()).setZero()}, // use setZero instead
        m_cg_i{0} {
    if (spdlog::get("cg")) {
      m_log = spdlog::get("cg");
    } else {
      m_log = spdlog::basic_logger_mt("cg", "_cg.log", true);
    }
    m_log->set_pattern("[%l] [CG] %v");
  }
  //! Conjugant Gradient deconstructor
  ~ConjugateGradients() = default;

  //! Calls the next step in the algorithm
  /**
   * Either calls the single_step or line_search method depending on the
   * parameters \return whether or not the algorithm has converged
   */
  bool stepImpl(ScalarType);
  bool runOptImpl(size_t, ScalarType);
  //! Gets the direction of the next step
  VectorType getStep();

private:
  const ConjugateGradients::Params m_p;
  //! Current step direction of the conjugate gradient
  VectorType m_direction;
  //! Algorithms previous step direction
  VectorType m_directionOld;
  //! Normalised version of the current direction vector
  VectorType m_directionNorm;
  //! Current force vector
  VectorType m_force;
  //! Previous force vector
  VectorType m_forceOld;
  std::shared_ptr<spdlog::logger> m_log;

  //! Counts the number of discrete steps until algorithm convergence
  size_t m_cg_i;

  //! Steps the conjugate gradient
  /**
   *  Checks for convergence based on the change in displacement or direction
   */
  int single_step(double a_maxMove);
  //! Steps the conjugate gradient
  /**
   * Checks for convergence based on the ratio of the projected force along a
   * line to the norm of the total force
   */
  int line_search(double a_maxMove);
};

} // namespace eonc
