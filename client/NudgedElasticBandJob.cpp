#include "NudgedElasticBandJob.h"
#include "ConjugateGradients.h"
#include "Potential.h"

using namespace std;

std::vector<std::string> NudgedElasticBandJob::run(void) {
  NudgedElasticBand::NEBStatus status;
  int f1;

  string reactantFilename = helper_functions::getRelevantFile("reactant.con");
  string productFilename = helper_functions::getRelevantFile("product.con");

  string transitionStateFilename = helper_functions::getRelevantFile("ts.con");
  bool tsInterpolate = false;
  auto transitionState = std::make_shared<Matter>(pot, params);
  FILE *fhTransitionState = fopen("ts.con", "r");
  if (fhTransitionState != nullptr) {
    tsInterpolate = true;
    fclose(fhTransitionState);
    transitionState->con2matter(transitionStateFilename);
  } else {
    transitionState = nullptr;
  }

  auto initial = std::make_shared<Matter>(pot, params);
  auto final_state = std::make_shared<Matter>(pot, params);

  initial->con2matter(reactantFilename);
  final_state->con2matter(productFilename);

  if (params->nebMinimEP && params->nebIpath == ""s) {
    SPDLOG_LOGGER_DEBUG(m_log, "Minimizing reactant");
    // TODO(rg): Maybe when we have even more parameters, false can be set by
    // the user too..
    initial->relax(false, params->writeMovies, params->checkpoint, "react_neb",
                   "react_neb");
    // TODO(rg): How do we report the total E/F now? Currently this is just the
    // total total, people might want "per-stage" totals (but they can also get
    // them from the log.)
    SPDLOG_LOGGER_DEBUG(m_log, "Minimized reactant in ");
    SPDLOG_LOGGER_DEBUG(m_log, "Minimizing product");
    final_state->relax(false, params->writeMovies, params->checkpoint,
                       "prod_neb", "prod_neb");
  }

  auto neb =
      std::make_unique<NudgedElasticBand>(initial, final_state, params, pot);

  if (tsInterpolate) {
    AtomMatrix reactantToTS = transitionState->pbc(
        transitionState->getPositions() - initial->getPositions());
    AtomMatrix TSToProduct = transitionState->pbc(
        final_state->getPositions() - transitionState->getPositions());
    for (int image = 1; image <= neb->numImages; image++) {
      int mid = neb->numImages / 2 + 1;
      if (image < mid) {
        double frac = ((double)image) / ((double)mid);
        neb->path[image]->setPositions(initial->getPositions() +
                                       frac * reactantToTS);
      } else if (image > mid) {
        double frac =
            (double)(image - mid) / (double)(neb->numImages - mid + 1);
        neb->path[image]->setPositions(transitionState->getPositions() +
                                       frac * TSToProduct);
      } else if (image == mid) {
        neb->path[image]->setPositions(transitionState->getPositions());
      }
    }
  }

  // f1 = Potential::fcalls;
  status = neb->compute();
  // fCallsNEB += Potential::fcalls - f1;

  if (status == NudgedElasticBand::NEBStatus::GOOD) {
    neb->printImageData();
    neb->findExtrema();
  }

  printEndState(status);
  saveData(status, neb.get());

  return returnFiles;
}

void NudgedElasticBandJob::saveData(NudgedElasticBand::NEBStatus status,
                                    NudgedElasticBand *neb) {
  FILE *fileResults, *fileNEB;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");

  fprintf(fileResults, "%d termination_reason\n", static_cast<int>(status));
  fprintf(
      fileResults, "%s potential_type\n",
      std::string{magic_enum::enum_name<PotType>(params->potential)}.c_str());
  // fprintf(fileResults, "%ld total_force_calls\n", Potential::fcalls);
  // fprintf(fileResults, "%ld force_calls_neb\n", fCallsNEB);
  fprintf(fileResults, "%f energy_reference\n",
          neb->path[0]->getPotentialEnergy());
  fprintf(fileResults, "%li number_of_images\n", neb->numImages);
  for (long i = 0; i <= neb->numImages + 1; i++) {
    fprintf(fileResults, "%f image%li_energy\n",
            neb->path[i]->getPotentialEnergy() -
                neb->path[0]->getPotentialEnergy(),
            i);
    fprintf(fileResults, "%f image%li_force\n",
            neb->path[i]->getForces().norm(), i);
    fprintf(fileResults, "%f image%li_projected_force\n",
            neb->projectedForce[i]->norm(), i);
  }
  fprintf(fileResults, "%li number_of_extrema\n", neb->numExtrema);
  for (long i = 0; i < neb->numExtrema; i++) {
    fprintf(fileResults, "%f extremum%li_position\n", neb->extremumPosition[i],
            i);
    fprintf(fileResults, "%f extremum%li_energy\n", neb->extremumEnergy[i], i);
  }

  fclose(fileResults);

  std::string nebFilename("neb.con");
  returnFiles.push_back(nebFilename);
  fileNEB = fopen(nebFilename.c_str(), "wb");
  for (long i = 0; i <= neb->numImages + 1; i++) {
    neb->path[i]->matter2con(fileNEB);
  }
  fclose(fileNEB);

  returnFiles.push_back("neb.dat");
  neb->printImageData(true);
}

void NudgedElasticBandJob::printEndState(NudgedElasticBand::NEBStatus status) {
  SPDLOG_LOGGER_DEBUG(m_log, "Final state: ");
  if (status == NudgedElasticBand::NEBStatus::GOOD)
    SPDLOG_LOGGER_DEBUG(m_log, "Nudged elastic band, successful.");
  else if (status == NudgedElasticBand::NEBStatus::BAD_MAX_ITERATIONS)
    SPDLOG_LOGGER_DEBUG(m_log, "Nudged elastic band, too many iterations.");
  else
    SPDLOG_LOGGER_WARN(m_log, "Unknown status: {}!", static_cast<int>(status));
  return;
}
