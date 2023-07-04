#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "Eigen.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"

/** @defgroup Optimizers
 *
 * \brief ClientEON methods for optimizing atomic structures
 *
 * This page provides links to all the available optimizers that can be run by
 * the ClientEON, as well as documentation on the optimizer class.
 *
 */

/**
 * @file
 * @ingroup Optimizers
 *
 * \brief The optimizer class is used to serve as an abstract class for all
 * optimizers, as well as to call an optimizer at runtime based off of the
 * passed in paramters.
 *
 * The set of optimizers are methods for optimizing atomic structures, solving
 * unconstrained energy minimization. Only a certain set of job types that
 * ClientEON runs can take advantage of the numeric optimizer (SEE OVERVIEW) and
 * are documented in their own files accordingly.
 *
 */

/**
 * Declaration of the optimizer class
 */

class Optimizer {
private:
  const OptType otype;

protected:
  std::shared_ptr<Parameters> m_params;

public:
  Optimizer(std::shared_ptr<Parameters> a_params)
      : otype{a_params->potential},
        m_params{a_params} {}
  //! optimizer deconstructor
  virtual ~Optimizer(){};
  //! Template for stepping the optimizer, returns convergence
  virtual int step(double maxMove) = 0;
  //! Template for runnning the optimizer; uses a series of steps, checking for
  //! convergence each time
  virtual int run(int maxIterations, double maxMove) = 0;
  //! Grabs the correct optimizer as specified by the parameters
  /*!
   * \param *objf an ref ObjectiveFunction that tells the optimizer how to run
   * \param *parameters defined by the config.init file
   */
  static std::unique_ptr<Optimizer>
  getOptimizer(std::shared_ptr<ObjectiveFunction> a_objf,
               std::shared_ptr<Parameters> a_params, bool a_refine = false);
};

#endif
