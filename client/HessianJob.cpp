#include "HessianJob.h"
#include "Hessian.h"
#include "Matter.h"
#include "Potential.h"

std::vector<std::string> HessianJob::run(void) {
  string matter_in("pos.con");

  std::vector<std::string> returnFiles;

  auto matter = std::make_unique<Matter>(pot, params);

  matter->con2matter(matter_in);

  Hessian hessian(params.get(), matter.get());
  long nAtoms = matter->numberOfAtoms();

  VectorXi moved(nAtoms);
  moved.setConstant(-1);

  int nMoved = 0;
  for (int i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      moved[nMoved] = i;
      nMoved++;
    }
  }
  moved = moved.head(nMoved);
  hessian.getFreqs(matter.get(), moved);

  FILE *fileResults;
  //    FILE *fileMode;

  std::string results_file("results.dat");

  returnFiles.push_back(results_file);

  fileResults = fopen(results_file.c_str(), "wb");

  // fprintf(fileResults, "%d force_calls\n", Potential::fcalls);
  fclose(fileResults);

  return returnFiles;
}
