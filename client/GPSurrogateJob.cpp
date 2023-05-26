#include "GPSurrogateJob.h"
#include "Matter.h"

std::vector<std::string> GPSurrogateJob::run(void) {
  std::vector<std::string> returnFiles;

  // Start working

  string reactantFilename = helper_functions::getRelevantFile("reactant.con");
  string productFilename = helper_functions::getRelevantFile("product.con");

  Matter *initial = new Matter(params);
  Matter *final = new Matter(params);

  initial->con2matter(reactantFilename);
  final->con2matter(productFilename);

  string posInFilename("pos.con");
  string resultsFilename("results.dat");
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
