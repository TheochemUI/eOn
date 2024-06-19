//----------------------------------------------------------------------------

#ifndef SADDLESEARCHJOB_H
#define SADDLESEARCHJOB_H

#include "Job.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "Parameters.h"

/**
 * @file
 * @ingroup Jobs
 *
 * \brief Finds transition states by finding saddle points on the potential
 * energy surface
 *
 * The saddle seach job implements a ref MinModeSaddleSearch, as well as an \ref
 * Optimizer "optimizer". A saddle search is initiated by making a local
 * displacement of atoms from their position at the minimum of the current
 * state. It can either be run using the ref Dimer or the ref Lanczos min-mode
 * method.
 *
 * The saddle seach job also optionally generates and works with a GP
 * approximated potential energy system with the GPR Dimer minimum mode method
 */

/**
 * Declaration of the Saddle Search job
 */

class SaddleSearchJob : public Job {
public:
  //! Saddle Search job constructor
  /*!
   * \param *params defined by the config.init file
   */
  SaddleSearchJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)),
        fCallsSaddle{0} {
    log = spdlog::get("combi");
  }
  //! Saddle Search Job Deconstructor
  ~SaddleSearchJob(void) = default;
  //! Kicks off the Saddle Search
  std::vector<std::string> run(void);

private:
  //! Runs the correct saddle search; also checks if the run was successful
  int doSaddleSearch();
  //! Logs the run status and makes sure the run was successful
  void printEndState(int status);
  //! Writes the results from the run to file
  void saveData(int status);

  //! Container for the results of the run
  std::vector<std::string> returnFiles;

  //! Initializes a ref MinModeSaddleSearch
  std::unique_ptr<MinModeSaddleSearch> saddleSearch;
  //! Initial configuration.
  std::shared_ptr<Matter> initial;
  //! Configuration used during the saddle point search.
  std::shared_ptr<Matter> saddle;
  //! Configuration used during the saddle point search.
  std::shared_ptr<Matter> displacement;

  //! Force calls to find the saddle
  int fCallsSaddle;

  std::shared_ptr<spdlog::logger> log;
};

#endif
