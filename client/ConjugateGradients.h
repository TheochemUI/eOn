#ifndef CG_H
#define CG_H

#include "Eigen.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "Optimizer.h"
#include "Parameters.h"

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

class ConjugateGradients : public Optimizer {
public:
  //! Conjugate Gradients optimizer constructor
  /*!
   * \param std::shared_ptr<ObjectiveFunction> m_objf that tells the optimizer
   * how to run \param std::shared_ptr<Parameters> m_params defined by the
   * config.init file
   */
  ConjugateGradients(std::shared_ptr<ObjectiveFunction> a_objf,
                     std::shared_ptr<Parameters> a_params)
      : m_objf{a_objf},
        m_params{a_params},
        m_directionOld{a_objf->getPositions() * 0.0},
        m_forceOld{a_objf->getPositions() * 0.0}, // use setZero instead
        m_cg_i{0} {
    if (spdlog::get("cg")) {
      m_log = spdlog::get("cg");
    } else {
      m_log = spdlog::basic_logger_st("cg", "_cg.log", true);
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
  int step(double maxMove);
  //! Runs the conjugate gradient
  /**
   * \todo method should also return an error code and message if the algorithm
   * errors out \return algorithm convergence
   */
  int run(int maxIterations, double maxMove);
  //! Gets the direction of the next step
  Eigen::VectorXd getStep();

private:
  //! An objective function relating a certain job method to the conjugate
  //! gradient optimizer
  std::shared_ptr<ObjectiveFunction> m_objf;
  //! Parameters set by the config.init file
  std::shared_ptr<Parameters> m_params;

  //! Current step direction of the conjugate gradient
  Eigen::VectorXd m_direction;
  //! Algorithms previous step direction
  Eigen::VectorXd m_directionOld;
  //! Normalised version of the current direction vector
  Eigen::VectorXd m_directionNorm;
  //! Current force vector
  Eigen::VectorXd m_force;
  //! Previous force vector
  Eigen::VectorXd m_forceOld;
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

#endif
