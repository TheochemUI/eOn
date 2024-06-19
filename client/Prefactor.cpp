#include "Prefactor.h"
#include "Hessian.h"

int Prefactor::getPrefactors(Parameters *parameters, Matter *min1,
                             Matter *saddle, Matter *min2, double &pref1,
                             double &pref2) {
  VectorXd min1Freqs, saddleFreqs, min2Freqs;

  // determine which atoms moved in the process
  VectorXi atoms;

  if (parameters->prefactorFilterScheme == Prefactor::FILTER_FRACTION) {
    atoms = movedAtomsPct(parameters, min1, saddle, min2);
  } else {
    atoms = movedAtoms(parameters, min1, saddle, min2);
  }

  int size = 3 * atoms.rows();
  assert(size > 0);

  // calculate min1 frequencies
  Hessian hessian(parameters, min1);
  min1Freqs = hessian.getFreqs(min1, atoms);
  if (min1Freqs.size() == 0) {
    SPDLOG_ERROR("[Prefactor] Bad hessian: min1");
    return -1;
  }
  // remove zero modes
  if (parameters->checkRotation) {
    min1Freqs = hessian.removeZeroFreqs(min1Freqs);
  }

  // calculate saddle frequencies
  saddleFreqs = hessian.getFreqs(saddle, atoms);
  if (saddleFreqs.size() == 0) {
    SPDLOG_ERROR("[Prefactor] Bad hessian: saddle");
    return -1;
  }
  // remove zero modes
  if (parameters->checkRotation) {
    saddleFreqs = hessian.removeZeroFreqs(saddleFreqs);
  }

  // calculate min2 frequencies
  min2Freqs = hessian.getFreqs(min2, atoms);
  if (min2Freqs.size() == 0) {
    if (!parameters->quiet) {
      SPDLOG_ERROR("[Prefactor] Bad hessian: min2");
    }
    return -1;
  }
  // remove zero modes
  if (parameters->checkRotation) {
    min2Freqs = hessian.removeZeroFreqs(min2Freqs);
  }

  // check Hessian sizes
  if ((min1Freqs.size() != saddleFreqs.size()) ||
      (min1Freqs.size() != saddleFreqs.size())) {
    if (!parameters->quiet) {
      SPDLOG_ERROR("[Prefactor] Bad prefactor: Hessian sizes do not match");
    }
    return -1;
  }

  logFreqs(min1Freqs, (char *)"minimum 1");
  logFreqs(saddleFreqs, (char *)"saddle");
  logFreqs(min2Freqs, (char *)"minimum 2");

  // check for correct number of negative modes
  int i, numNegFreq = 0;
  for (i = 0; i < size; i++) {
    if (min1Freqs(i) < 0) {
      SPDLOG_DEBUG("[Prefactor] min1 had negative mode of {}", min1Freqs(i));
      numNegFreq++;
    }
  }
  if (numNegFreq != 0) {
    SPDLOG_DEBUG("[Prefactor] Error: {} negative modes at min1", numNegFreq);
    return -1;
  }

  numNegFreq = 0;
  for (i = 0; i < size; i++) {
    if (saddleFreqs(i) < 0) {
      numNegFreq++;
    }
  }
  if (numNegFreq != 1) {
    SPDLOG_DEBUG("Error: {} negative modes at saddle", numNegFreq);
    return -1;
  }

  numNegFreq = 0;
  for (i = 0; i < size; i++) {
    if (min2Freqs(i) < 0) {
      numNegFreq++;
    }
  }
  if (numNegFreq != 0) {
    SPDLOG_DEBUG("Error: {} negative modes at min2", numNegFreq);
    return -1;
  }

  // calculate the prefactors
  pref1 = 1.0;
  pref2 = 1.0;

  if (parameters->prefactorRate == Prefactor::RATE_HTST) {

    // products are calculated this way in order to avoid overflow
    for (int i = 0; i < saddleFreqs.size(); i++) {
      pref1 *= min1Freqs[i];
      pref2 *= min2Freqs[i];
      if (saddleFreqs[i] > 0) {
        pref1 /= saddleFreqs[i];
        pref2 /= saddleFreqs[i];
      }
    }
    pref1 = sqrt(pref1) / (2 * M_PI * 10.18e-15);
    pref2 = sqrt(pref2) / (2 * M_PI * 10.18e-15);
  } else if (parameters->prefactorRate == Prefactor::RATE_QQHTST) {
    float kB_T = parameters->temperature * 8.617332e-5; // eV
    float h_bar = 6.582119e-16;                         // eV*s
    float h = 4.135667e-15;                             // eV*s
    float temp = (h_bar / (2. * kB_T));

    for (int i = 0; i < min1Freqs.size(); i++) {
      pref1 = pref1 * (sinh(temp * (sqrt(min1Freqs[i]) / 10.18e-15)));
      pref2 = pref2 * (sinh(temp * (sqrt(min2Freqs[i]) / 10.18e-15)));

      if (saddleFreqs[i] > 0) {
        pref1 = pref1 / (sinh(temp * (sqrt(saddleFreqs[i]) / 10.18e-15)));
        pref2 = pref2 / (sinh(temp * (sqrt(saddleFreqs[i]) / 10.18e-15)));
      }
    }
    pref1 = 2. * kB_T / (h)*pref1;
    pref2 = 2. * kB_T / (h)*pref2;
  }
  SPDLOG_DEBUG("[Prefactor] reactant to product prefactor: {:.3e}", pref1);
  SPDLOG_DEBUG("[Prefactor] product to reactant prefactor: {:.3e}", pref2);
  return 0;
}

void Prefactor::logFreqs(VectorXd freqs, char *name) {
  std::shared_ptr<spdlog::logger> fileLogger;
  fileLogger = spdlog::basic_logger_mt("prefactor", "freqs.dat", true);
  fileLogger->set_pattern("%v");
  fileLogger->debug("[Prefactor] Frequencies at {}", name);
  int i;
  for (i = 0; i < freqs.size(); i++) {
    fileLogger->debug("");
    fileLogger->debug("{:10.6f} ", freqs(i));
    if ((i + 1) % 5 == 0) {
      fileLogger->debug("");
    }
  }
  fileLogger->debug("");
  spdlog::drop("prefactor");
  fileLogger.reset();
}

VectorXi Prefactor::movedAtoms(Parameters *parameters, Matter *min1,
                               Matter *saddle, Matter *min2) {
  long nAtoms = saddle->numberOfAtoms();

  VectorXi moved(nAtoms);
  moved.setConstant(-1);

  AtomMatrix diffMin1 =
      saddle->pbc(saddle->getPositions() - min1->getPositions());
  AtomMatrix diffMin2 =
      saddle->pbc(saddle->getPositions() - min2->getPositions());

  diffMin1.array() *= saddle->getFree().array();
  diffMin2.array() *= saddle->getFree().array();

  int nMoved = 0;
  for (int i = 0; i < nAtoms; i++) {
    if ((diffMin1.row(i).norm() > parameters->prefactorMinDisplacement) ||
        (diffMin2.row(i).norm() > parameters->prefactorMinDisplacement)) {
      if (!(moved.array() == i).any()) {
        moved[nMoved] = i;
        nMoved++;
      }
      for (int j = 0; j < nAtoms; j++) {
        double diffRSaddle = saddle->distance(i, j);

        if (diffRSaddle < parameters->prefactorWithinRadius &&
            (!saddle->getFixed(j))) {
          if (!(moved.array() == j).any()) {
            moved[nMoved] = j;
            nMoved++;
          }
        }
      }
    }
  }
  return (VectorXi)moved.block(0, 0, nMoved, 1);
}

VectorXi Prefactor::movedAtomsPct(Parameters *parameters, Matter *min1,
                                  Matter *saddle, Matter *min2) {
  long nAtoms = saddle->numberOfAtoms();
  long nFree = saddle->numberOfFreeAtoms();

  VectorXi moved(nAtoms);
  moved.setConstant(-1);

  AtomMatrix diffMin1 =
      saddle->pbc(saddle->getPositions() - min1->getPositions());
  AtomMatrix diffMin2 =
      saddle->pbc(saddle->getPositions() - min2->getPositions());

  diffMin1.array() *= saddle->getFree().array();
  diffMin2.array() *= saddle->getFree().array();

  VectorXd diff(nAtoms);
  diff.setConstant(0.0);

  SPDLOG_DEBUG(
      "[Prefactor] including all atoms that make up {:.3f}% of the motion",
      parameters->prefactorFilterFraction * 100);
  double sum = 0.0;
  int mini = 0;
  for (int i = 0; i < nAtoms; i++) {
    diff[i] = max(diffMin1.row(i).norm(), diffMin2.row(i).norm());
    sum += diff[i];
    if (diff[i] <= diff[mini]) {
      mini = i;
    }
  }

  SPDLOG_DEBUG("[Prefactor] sum of atom distances moved {:.4f}", sum);
  SPDLOG_DEBUG("[Prefactor] max moved atom distance: {:.4f}", diff.maxCoeff());

  int nMoved = 0;
  double d = 0.0;
  while (d / sum <= parameters->prefactorFilterFraction && nMoved < nFree) {
    int maxi = mini;
    for (int i = 0; i < nAtoms; i++) {
      if (diff[i] >= diff[maxi]) {
        if (!(moved.array() == i).any()) {
          maxi = i;
        }
      }
    }
    moved[nMoved] = maxi;
    nMoved++;
    d += diff[maxi];
  }

  int totalAtoms = nMoved;
  for (int i = 0; i < nMoved; i++) {
    for (int j = 0; j < nAtoms; j++) {
      if (moved[i] == j)
        continue;

      double diffRSaddle = saddle->distance(moved[i], j);

      if (diffRSaddle < parameters->prefactorWithinRadius &&
          (!saddle->getFixed(j))) {

        if (!(moved.array() == j).any()) {
          moved[totalAtoms++] = j;
        }
      }
    }
  }
  SPDLOG_DEBUG("[Prefactor] including {} atoms in the hessian ({} moved + {} "
               "neighbors)",
               totalAtoms, nMoved, totalAtoms - nMoved);
  return (VectorXi)moved.block(0, 0, totalAtoms, 1);
}

VectorXi Prefactor::allFreeAtoms(Matter *matter) {
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
  return moved.head(nMoved);
}
