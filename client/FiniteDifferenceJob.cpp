#include "FiniteDifferenceJob.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "Matter.h"

using namespace helper_functions;

std::vector<std::string> FiniteDifferenceJob::run(void) {
  // No bundling for this job, so bundleNumber is ignored.

  // Load the displacement con file and get the position.
  auto reactant = std::make_unique<Matter>(pot, params);
  reactant->con2matter("pos.con");
  AtomMatrix posA;
  posA = reactant->getPositions();

  double dRs[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 5e-3, 0.01, 0.05, 0.1, -1};

  AtomMatrix forceA;
  forceA = reactant->getForces();

  // Create a random displacement.
  long epicenter = EpiCenters::minCoordinatedEpiCenter(reactant.get(),
                                                       params->neighborCutoff);
  AtomMatrix displacement;
  displacement.resize(reactant->numberOfAtoms(), 3);
  displacement.setZero();
  printf("displacing atoms:");
  for (int i = 0; i < reactant->numberOfAtoms(); i++) {
    if (reactant->distance(epicenter, i) <= 3.3) {
      printf(" %i", i);
      for (int j = 0; j < 3; j++) {
        if (!reactant->getFixed(i)) {
          displacement(i, j) = randomDouble(1.0);
        }
      }
    }
  }
  printf("\n");
  displacement.normalize();

  // Loop over values of dimer dR and print the output to results.dat.
  FILE *results = fopen("results.dat", "w");
  fprintf(results, "%14s    %14s\n", "dR", "curvature");
  printf("%14s    %14s\n", "dR", "curvature");
  AtomMatrix posB;
  AtomMatrix forceB;
  double curvature = 0.0;
  for (int dRi = 0; dRs[dRi] != -1; dRi++) {
    posB = posA + displacement * dRs[dRi];
    reactant->setPositions(posB);
    forceB = reactant->getForces();
    curvature =
        ((forceB - forceA).array() * displacement.array()).sum() / dRs[dRi];
    fprintf(results, "%14.8f    %14.8f\n", dRs[dRi], curvature);
    printf("%14.8f    %14.8f\n", dRs[dRi], curvature);
    fflush(results);
  }
  fclose(results);

  std::vector<std::string> empty;
  return empty;
}
