
#ifndef PROCESSSEARCHJOB_H
#define PROCESSSEARCHJOB_H

#include "Job.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "Parameters.h"

/**
 * @file
 * @ingroup Jobs
 *
 * \brief Finds possible escape mecahnisms from a state.
 *
 * The process search job implements one of the following types of saddle
 * searches, as well as an \ref Optimizer "optimizer", as defined by the
 * config.init file.
 *
 * <ul>
 * <li> ref MinModeSaddleSearch
 * <li> ref BasinHopingSaddleSearch
 * <li> ref DynamicsSaddleSearch
 * <li> ref BiasedGradientSquaredDescent
 * </ul>
 *
 */

/**
 * Declaration of the Process Search job
 */

class ProcessSearchJob : public Job {
public:
  //! Process Search job Constructor
  /*!
   * \param *params defined by the config.init file
   */
  ProcessSearchJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)),
        fCallsSaddle{0},
        fCallsMin{0},
        fCallsPrefactors{0} {
    log = spdlog::get("combi");
  }
  //! Process Search job De-constructor
  ~ProcessSearchJob() = default;
  //! Kicks off the Process Search
  std::vector<std::string> run(void) override;

private:
  shared_ptr<spdlog::logger> log;
  //! Runs the correct saddle search; also checks if the run was successful
  int doProcessSearch(void);
  //! UNDEFINED
  /*!
   * Function is not defined in the cpp file and is marked for deletion.
   */
  VectorXi movedAtoms(void);
  //! Logs the run status and makes sure the run was successful
  void printEndState(int status);
  //! Writes the results from the run to file
  void saveData(int status);

  //! Container for the results of the run
  std::vector<std::string> returnFiles;

  //! Pulled from parameters
  /*!
   * sa ref BasinHoppingSaddleSearch DynamicsSaddleSearch
   * BiasedGradientSquaredDescent MinModeSaddleSearch
   */
  std::unique_ptr<SaddleSearchMethod> saddleSearch;
  //! Initial configuration
  std::shared_ptr<Matter> initial;
  //! Configuration used during the saddle point search
  std::shared_ptr<Matter> saddle;
  //! Configuration used during the saddle point search
  std::shared_ptr<Matter> displacement;
  //! First minimum from the saddle
  std::shared_ptr<Matter> min1;
  //! Second minimum from the saddle
  std::shared_ptr<Matter> min2;

  //! Array containg the value of the potential barriers from reactant to
  //! product and vice versa
  double barriersValues[2];
  //! Array containg the value of the prefactors of the matter from reactant to
  //! product and vice versa
  double prefactorsValues[2];

  //! Force calls to find the saddle
  size_t fCallsSaddle;
  //! Force calls to minimize
  size_t fCallsMin;
  //! Force calls to find the prefactors
  size_t fCallsPrefactors;
};

#endif
