#include "MonteCarloJob.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "MonteCarlo.h"

std::vector<std::string> MonteCarloJob::run(void) {
  string posInFilename("pos.con");
  string posOutFilename("out.con");

  if (params->checkpoint) {
    FILE *pos;
    pos = fopen("pos_cp.con", "r");
    if (pos != NULL) {
      posInFilename = "pos_cp.con";
      SPDLOG_LOGGER_DEBUG(log, "Resuming from checkpoint\n");
    } else {
      SPDLOG_LOGGER_DEBUG(log, "No checkpoint files found\n");
    }
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back(posOutFilename);

  auto matter = std::make_shared<Matter>(pot, params);
  matter->con2matter(posInFilename);

  // code will go
  MonteCarlo mc = MonteCarlo(matter, params);
  mc.run(params->monteCarloSteps, params->temperature,
         params->monteCarloStepSize);

  // FILE *fileResults;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  // fileResults = fopen(resultsFilename.c_str(), "wb");

  // fprintf(fileResults, "%d termination_reason\n", status);
  // fprintf(fileResults, "minimization job_type\n");
  // fprintf(fileResults, "%s potential_type\n",
  // helper_functions::getPotentialName(params->potential).c_str());
  // fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
  // if (status != STATUS_POTENTIAL_FAILED) {
  //     fprintf(fileResults, "%f potential_energy\n",
  //     pos->getPotentialEnergy());
  // }
  // fclose(fileResults);

  return returnFiles;
}
