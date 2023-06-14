#include "GPSurrogateJob.h"
#include "Matter.h"
#include <memory>

std::vector<std::string> GPSurrogateJob::run(void) {
  std::vector<std::string> returnFiles;

  // Start working

  std::string reactantFilename = helper_functions::getRelevantFile("reactant.con");
  std::string productFilename = helper_functions::getRelevantFile("product.con");

  if (params->potential == PotType::PYSURROGATE){
    params->potential = params->true_pot;
  }

  auto initial = std::make_unique<Matter>(params);
  auto final_state = std::make_unique<Matter>(params);

  initial->con2matter(reactantFilename);
  final_state->con2matter(productFilename);

  std::string posInFilename("pos.con");
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  Matter *pos = new Matter(params);
  pos->con2matter(posInFilename);

  printf("Energy:         %f\n", pos->getPotentialEnergy());
  printf("Max atom force: %g\n", pos->maxForce());

  FILE *fileResults = fopen(resultsFilename.c_str(), "wb");
  fprintf(fileResults, "%f Energy\n", pos->getPotentialEnergy());
  fprintf(fileResults, "%f Max_Force\n", pos->maxForce());
  fclose(fileResults);

  return returnFiles;
}
