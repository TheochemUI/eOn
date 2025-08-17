#include "PointJob.h"
#include "Matter.h"

std::vector<std::string> PointJob::run(void) {
  std::vector<std::string> returnFiles;
  string posInFilename("pos.con");
  string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  auto pos = std::make_unique<Matter>(pot, params);
  pos->con2matter(posInFilename);

  SPDLOG_LOGGER_DEBUG(log, "Energy:         {:.12f}",
                      pos->getPotentialEnergy());
  std::stringstream freeForcesStream;
  freeForcesStream << pos->getForcesFree();
  SPDLOG_LOGGER_DEBUG(log, "(free) Forces:\n{}", freeForcesStream.str());
  SPDLOG_LOGGER_DEBUG(log, "Max atom force: {:.12f}", pos->maxForce());

  std::shared_ptr<spdlog::logger> fileLogger;
  fileLogger = spdlog::basic_logger_mt("point", "results.dat", true);

  fileLogger->set_pattern("%v");
  fileLogger->info("{:.12f} Energy", pos->getPotentialEnergy());
  fileLogger->info("{:.12f} Max_Force", pos->maxForce());

  spdlog::drop("point");
  fileLogger.reset();
  return returnFiles;
}
