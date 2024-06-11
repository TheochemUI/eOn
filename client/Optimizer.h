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
  const OptType m_otype;

protected:
  //! Parameters set by the config.init file
  std::shared_ptr<Parameters> m_params;
  //! An objective function relating a certain job method to the conjugate
  //! gradient optimizer
  std::shared_ptr<ObjectiveFunction> m_objf;

public:
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf,
            std::shared_ptr<Parameters> a_params)
      : m_otype{a_params->optMethod},
        m_params{a_params},
        m_objf{a_objf} {
    SPDLOG_WARN("You should explicitly set an optimizer while constructing the "
                "optimizer!!\n Defaulting to opt_method from the parameters");
  }
  Optimizer(std::shared_ptr<ObjectiveFunction> a_objf, OptType a_optype,
            std::shared_ptr<Parameters> a_params)
      : m_otype{a_optype},
        m_params{a_params},
        m_objf{a_objf} {}
  //! optimizer deconstructor
  virtual ~Optimizer(){};
  //! Template for stepping the optimizer, returns convergence
  virtual int step(double a_maxMove) = 0;
  //! Template for running the optimizer; uses a series of steps, checking for
  //! convergence each time
  virtual int run(size_t a_maxIterations, double a_maxMove) = 0;
};

namespace helpers::create {
std::unique_ptr<Optimizer> mkOptim(std::shared_ptr<ObjectiveFunction> a_objf,
                                   OptType a_otype,
                                   std::shared_ptr<Parameters> a_params);
}

#endif
