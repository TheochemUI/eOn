#include "PointJob.h"
#include "Matter.h"

std::vector<std::string> PointJob::run(void) {
  std::vector<std::string> returnFiles;
  string posInFilename("pos.con");
  string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  auto pos = std::make_unique<Matter>(pot, params);
  pos->con2matter(posInFilename);

  SPDLOG_LOGGER_DEBUG(log, "Energy:         {}", pos->getPotentialEnergy());
  SPDLOG_LOGGER_DEBUG(log, "Max atom force: {}", pos->maxForce());

  std::shared_ptr<spdlog::logger> fileLogger;
  fileLogger = spdlog::basic_logger_mt("point", "results.dat", true);

  fileLogger->set_pattern("%v");
  fileLogger->info("{} Energy", pos->getPotentialEnergy());
  fileLogger->info("{} Max_Force", pos->maxForce());

  spdlog::drop("point");
  fileLogger.reset();
  return returnFiles;
}
