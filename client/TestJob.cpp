#include "TestJob.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "Parameters.h"
#include "Potential.h"

#include <stdlib.h>

std::vector<std::string> TestJob::run(void) {
  checkPotentials();
  checkFullSearch();
  std::vector<std::string> empty;
  return empty;
}

void TestJob::checkFullSearch(void) {
  printf("\n---Beginning tests of saddle point search---\n");
  printf("Checks the potential energies of located configurations.\n");
  printf("Reported as OK if within a tolerance of: %f\n", tolerance);

  //    long status;
  bool ok = 1;
  double diffM1, diffM2, diffSP;

  string reactantFilename("reactant_test.con");
  string displacementFilename("displacement_test.con");
  string modeFilename("mode_test.dat");

  params->potential = "emt";

  auto initial = std::make_unique<Matter>(pot, params);
  //    displacement = std::make_unique<Matter>(pot, params);
  auto saddle = std::make_unique<Matter>(pot, params);
  auto min1 = std::make_unique<Matter>(pot, params);
  auto min2 = std::make_unique<Matter>(pot, params);
  auto matterTemp = std::make_unique<Matter>(pot, params);

  saddle->con2matter(displacementFilename);
  initial->con2matter(reactantFilename);
  *min1 = *min2 = *initial;

  printf("\n---Output for saddle point search start---\n");
  //    saddleSearch = new SaddleSearch();
  //    saddleSearch->initialize(initial, saddle, params);
  //    saddleSearch->loadMode(mode_passed);
  //    status = saddleSearch->locate();
  printf("---Output for saddle point search end---\n\n");

  printf("---Output relax from saddle point search start---\n");
  // relax from the saddle point located

  //    AtomMatrix posSaddle = saddleSearch->getSaddlePositions();
  AtomMatrix displacedPos;

  *min1 = *saddle;
  // XXX: the distance displaced from the saddle should be a parameter
  //    displacedPos = posSaddle - saddleSearch->getEigenMode() * 0.2;
  min1->setPositions(displacedPos);
  //   ConjugateGradients cgMin1(min1, params);
  //  cgMin1.fullRelax();
  //  fCallsMin += cgMin1.totalForceCalls;

  *min2 = *saddle;
  //    displacedPos = posSaddle + saddleSearch->getEigenMode() * 0.2;
  //    min2->setPositions(displacedPos);
  //    ConjugateGradients cgMin2(min2, params);
  //    cgMin2.fullRelax();
  //  fCallsMin += cgMin2.totalForceCalls;

  // If min2 corresponds to initial state swap min1 && min2
  if (!initial->compare(min1) && initial->compare(min2)) {
    *matterTemp = *min1;
    *min1 = *min2;
    *min2 = *matterTemp;
  }
  printf("---Output relax from saddle point search end---\n");

  // checking the energies of the obtained configurations
  diffM1 = abs(min1->getPotentialEnergy() - 45.737426);
  diffM2 = abs(min2->getPotentialEnergy() - 45.737433);
  diffSP = abs(saddle->getPotentialEnergy() - 46.284511);

  if ((diffM1 < tolerance) and (diffM2 < tolerance) and (diffSP < tolerance)) {
    ok *= 1;
    printf("OK: Saddle search structural energies\n");
  } else {
    if (tolerance < diffSP) {
      ok *= 0;
      printf("WARNING: Saddle point not within energy tolerance: %f\n", diffSP);
    }
    if (tolerance < diffM2) {
      ok *= 0;
      printf("WARNING: Minimum 2 not within energy tolerance: %f\n", diffM2);
    }
    if (tolerance < diffM1) {
      ok *= 0;
      printf("WARNING: Minimum 1 not within energy tolerance: %f\n", diffM2);
    }
  }

  // checking the structures of the obtained configurations
  diffM1 = abs((min1->getPositions()).row(384).norm() - 19.123375);
  diffM2 = abs((min2->getPositions()).row(384).norm() - 19.527995);
  diffSP = abs((saddle->getPositions()).row(384).norm() - 19.026709);

  if ((diffM1 < tolerance) and (diffM2 < tolerance) and (diffSP < tolerance)) {
    ok *= 1;
    printf("OK: Saddle search, adatom positions\n");
  } else {
    if (tolerance < diffSP) {
      ok *= 0;
      printf(
          "WARNING: Saddle point, adatom not within position tolerance: %f\n",
          diffSP);
    }
    if (tolerance < diffM2) {
      ok *= 0;
      printf("WARNING: Minimum 2, adatom not within position tolerance: %f\n",
             diffM2);
    }
    if (tolerance < diffM1) {
      ok *= 0;
      printf("WARNING: Minimum 1, adatom not within position tolerance: %f\n",
             diffM1);
    }
  }

  if (ok) {
    printf("Saddle search tests all good\n");
  } else {
    printf("Saddle search tests there were WARNINGS\n");
  }
  printf("SP done\n");

  delete initial;
  delete saddle;
  delete min1;
  delete min2;
  delete matterTemp;

  return;
}
/*
void TestJob::checkMode(void){

    string reactant_passed("reactant_test.con");
    params->potential = 1;  // always LJ in test
    Matter *saddle = std::make_unique<Matter>(pot, params);
    saddle->con2matter(reactant_passed);

    LowestEigenmode* lowestEigenmode = new Dimer(saddle, params);

}
*/
void TestJob::checkPotentials(void) {
  double energyDiff;
  double forceDiff;

  printf("\n---Beginning tests of potentials---\n");
  printf("Checks the potential energy and the max force.\n");
  printf("Reported as OK if within a tolerance of: %f\n\n", tolerance);

  energyDiff = getEnergyDiff(Potential::POT_LJ, -1475.984331);
  if (abs(energyDiff) > tolerance) {
    printf("WARNING: LJ energy difference: %f\n", energyDiff);
  } else {
    forceDiff = getForceDiff(Potential::POT_LJ, 2.007213);
    if (abs(forceDiff) > tolerance) {
      printf("WARNING: LJ force difference: %f\n", forceDiff);
    } else {
      printf("OK: LJ\n");
    }
  }

  energyDiff = getEnergyDiff(Potential::POT_EMT, 46.086312);
  if (abs(energyDiff) > tolerance) {
    printf("WARNING: EMT energy difference: %f\n", energyDiff);
  } else {
    forceDiff = getForceDiff(Potential::POT_EMT, 0.357493);
    if (abs(forceDiff) > tolerance) {
      printf("WARNING: EMT force difference: %f\n", forceDiff);
    } else {
      printf("OK: EMT\n");
    }
  }

  energyDiff = getEnergyDiff(Potential::POT_EDIP, -1033.250950);
  if (abs(energyDiff) > tolerance) {
    printf("WARNING: EDIP energy difference: %f\n", energyDiff);
  } else {
    forceDiff = getForceDiff(Potential::POT_EDIP, 7.080115);
    if (abs(forceDiff) > tolerance) {
      printf("WARNING: EDIP force difference: %f\n", forceDiff);
    } else {
      printf("OK: EDIP\n");
    }
  }

  energyDiff = getEnergyDiff(Potential::POT_TERSOFF_SI, -1035.809985);
  if (abs(energyDiff) > tolerance) {
    printf("WARNING: Tersoff energy difference: %f\n", energyDiff);
  } else {
    forceDiff = getForceDiff(Potential::POT_TERSOFF_SI, 11.145002);
    if (abs(forceDiff) > tolerance) {
      printf("WARNING: Tersoff force difference: %f\n", forceDiff);
    } else {
      printf("OK: Tersoff\n");
    }
  }

  energyDiff = getEnergyDiff(Potential::POT_SW_SI, -1449.795645);
  if (abs(energyDiff) > tolerance) {
    printf("WARNING: SW energy difference: %f\n", energyDiff);
  } else {
    forceDiff = getForceDiff(Potential::POT_SW_SI, 2.530904);
    if (abs(forceDiff) > tolerance) {
      printf("WARNING: SW force difference: %f\n", forceDiff);
    } else {
      printf("OK: SW\n");
    }
  }

  energyDiff = getEnergyDiff(Potential::POT_LENOSKY_SI, -1410.679106);
  if (abs(energyDiff) > tolerance) {
    printf("Lenosky energy difference: %f\n", energyDiff);
  } else {
    forceDiff = getForceDiff(Potential::POT_LENOSKY_SI, 2.320168);
    if (abs(forceDiff) > tolerance) {
      printf("Lenosky force difference: %f\n", forceDiff);
    } else {
      printf("OK: Lenosky\n");
    }
  }

  energyDiff = getEnergyDiff(Potential::POT_EAM_AL, -1206.825825);
  if (abs(energyDiff) > tolerance) {
    printf("WARNING: Aluminum energy difference: %f\n", energyDiff);
  } else {
    forceDiff = getForceDiff(Potential::POT_EAM_AL, 0.000246);
    if (abs(forceDiff) > tolerance) {
      printf("WARNING: Aluminum force difference: %f\n", forceDiff);
    } else {
      printf("OK: Aluminum\n");
    }
  }

  energyDiff = getEnergyDiff(Potential::POT_QSC, -1232.806318);
  if (abs(energyDiff) > tolerance) {
    printf("WARNING: QSC energy difference: %f\n", energyDiff);
  } else {
    forceDiff = getForceDiff(Potential::POT_QSC, 0.673444);
    if (abs(forceDiff) > tolerance) {
      printf("WARNING: QSC force difference: %f\n", forceDiff);
    } else {
      printf("OK: QSC\n");
    }
  }

  energyDiff = getEnergyDiff(Potential::POT_TIP4P, 4063.865115);
  if (abs(energyDiff) > tolerance) {
    printf("WARNING: TIP4P energy difference: %f\n", energyDiff);
  } else {
    forceDiff = getForceDiff(Potential::POT_TIP4P, 73.655248);
    if (abs(forceDiff) > tolerance) {
      printf("WARNING: TIP4P force difference: %f\n", forceDiff);
    } else {
      printf("OK: TIP4P\n");
    }
  }
}

double TestJob::getEnergyDiff(string pot, double refEnergy) {
  string posFilename("pos_test.con");
  params->potential = pot;
  Matter *pos = std::make_unique<Matter>(pot, params);
  pos->con2matter(posFilename);
  //    printf("Energy: %f\n", pos->getPotentialEnergy());
  return pos->getPotentialEnergy() - refEnergy;
}

double TestJob::getForceDiff(string pot, double refForce) {
  string posFilename("pos_test.con");
  params->potential = pot;
  Matter *pos = std::make_unique<Matter>(pot, params);
  pos->con2matter(posFilename);
  //    printf("Force: %f\n", pos->maxForce());
  return pos->maxForce() - refForce;
}
