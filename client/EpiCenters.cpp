#include "EpiCenters.h"
#include "HelperFunctions.h"

#include <cassert>
#include <climits>
#include <vector>

using namespace helper_functions;
using std::vector;

long EpiCenters::cnaEpiCenter(const Matter *matter, double neighborCutoff) {
  long *cnaList;
  long j, nAtoms, indexEpiCenter;
  double tempDouble;
  nAtoms = matter->numberOfAtoms();
  cnaList = new long[nAtoms];
  indexEpiCenter = -2; // initialize to a value that will fail the assert if no
                       // EpiCenter is found
  //----- Initialize end -----
  // std::cout<<"cnaEpiCenter\n";

  cna(cnaList, matter, neighborCutoff);

  // count atoms that are not FCC or HCP and are free to move
  j = 0;
  for (int i = 0; i < nAtoms; i++) {
    if ((cnaList[i] == 2) && !(matter->getFixed(i)))
      j++;
  }
  // pick a random atom being both free and not FCC or HCP coordinated
  tempDouble = randomDouble(j);
  j = (long)tempDouble + 1;
  for (int i = 0; i < nAtoms; i++) {
    if ((cnaList[i] == 2) && !(matter->getFixed(i))) {
      j--;
      if (!j) {
        indexEpiCenter = i;
        break;
      }
    }
  }
  delete[] cnaList;

  // make certain the chosen atom index is within the domain of the atom list
  assert(indexEpiCenter > -1 && indexEpiCenter < nAtoms);
  return (indexEpiCenter);
}

long EpiCenters::minCoordinatedEpiCenter(const Matter *matter,
                                         double neighborCutoff) {
  bool *minCoordinatedList;
  long j, nAtoms, indexEpiCenter, minCoordinationVal;
  double tempDouble;
  nAtoms = matter->numberOfAtoms();
  minCoordinatedList = new bool[nAtoms];
  indexEpiCenter = -2; // initialize to a value that will fail the assert if no
                       // EpiCenter is found
  //----- Initialize end -----
  // std::cout<<"minCoordinatedEpiCenter\n";

  minCoordinationVal = minCoordination(matter, neighborCutoff);
  coordinationLessOrEqual(minCoordinatedList, minCoordinationVal, matter,
                          neighborCutoff);

  // count all atoms that are minimally coordinated and free to move
  j = 0;
  for (int i = 0; i < nAtoms; i++) {
    if ((minCoordinatedList[i]) && !(matter->getFixed(i)))
      j++;
  }
  // pick a random atom that is free and minimally coordinated
  tempDouble = randomDouble(j);
  j = (long)tempDouble;
  for (int i = 0; i < nAtoms; i++) {
    if ((minCoordinatedList[i]) && !(matter->getFixed(i))) {
      if (!j) {
        indexEpiCenter = i;
        break;
      } else
        j--;
    }
  }
  delete[] minCoordinatedList;

  // make certain the chosen atom is in the atom list
  assert(indexEpiCenter > -1 && indexEpiCenter < nAtoms);
  return (indexEpiCenter);
}

long EpiCenters::lastAtom(const Matter *matter) {
  long nAtoms, indexEpiCenter;
  nAtoms = matter->numberOfAtoms();
  indexEpiCenter = nAtoms - 1;

  // make certain the chosen atom is in the atom list
  assert(indexEpiCenter > -1 && indexEpiCenter < nAtoms);
  return (indexEpiCenter);
}

long EpiCenters::randomFreeAtomEpiCenter(const Matter *matter) {
  long j, nAtoms, indexEpiCenter;
  double tempDouble;
  nAtoms = matter->numberOfAtoms();
  indexEpiCenter = -2; // initialize to a value that will fail the assert if no
                       // EpiCenter is found
  j = 0;

  j = matter->numberOfFreeAtoms() - 1;
  tempDouble = randomDouble(j);
  j = (long)tempDouble;

  // pick a random atom that is free
  for (int i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      if (!j) {
        indexEpiCenter = i;
        break;
      } else
        j--;
    }
  }
  // make certain the chosen atom is in the atom list
  assert(indexEpiCenter > -1 && indexEpiCenter < nAtoms);
  return (indexEpiCenter);
}

// long EpiCenters::randomFreeAtomEpiCenter(const Matter *matter)
//{
//     long indexEpiCenter;
//     long nAtoms = matter->numberOfAtoms();
//     do
//     {
//         indexEpiCenter = (long)randomDouble(nAtoms - 1);
//     } while (matter->getFixed(indexEpiCenter));
//     return(indexEpiCenter);
// }

void EpiCenters::cna(long *cna, const Matter *matter, double neighborCutoff) {
  int a1 = 0;
  int a2 = 0;
  int a3 = 0;
  int nAtoms;
  int unsigned n = 0;
  int unsigned n1 = 0;
  int unsigned n2 = 0;
  int unsigned m2 = 0;
  int j1 = 0;
  int j2 = 0;
  double diffR;
  int nBonds = 0;
  int bondsSum = 0;
  nAtoms = (int)matter->numberOfAtoms();
  vector<int> nFCC(nAtoms);
  vector<int> nHCP(nAtoms);
  vector<vector<int>> neighborLists(nAtoms);
  //----- Initialize end -----
  // std::cout<<"cna\n";

  for (int i = 0; i < nAtoms - 1; i++) {
    for (int j = i + 1; j < nAtoms; j++) {
      diffR = matter->distance(i, j);
      if (diffR < neighborCutoff) {
        neighborLists[i].push_back(j);
        neighborLists[j].push_back(i);
      }
    }
  }
  for (a2 = 0; a2 < nAtoms; a2++) {
    vector<int> &nbs2 = neighborLists[a2];
    for (n2 = 0; n2 < nbs2.size(); n2++) {
      a1 = nbs2[n2];
      if (a1 < a2) {
        vector<int> common;
        vector<int> &nbs1 = neighborLists[a1];
        for (n1 = 0; n1 < nbs1.size(); n1++) {
          a3 = nbs1[n1];
          for (m2 = 0; m2 < nbs2.size(); m2++)
            if (a3 == nbs2[m2])
              common.push_back(a3);
        }
        if (common.size() == 4) {
          nBonds = 0;
          bondsSum = 0;
          for (j2 = 1; j2 < 4; j2++) {
            vector<int> &nbs = neighborLists[common[j2]];
            for (j1 = 0; j1 < j2; j1++) {
              for (n = 0; n < nbs.size(); n++) {
                if (common[j1] == nbs[n]) {
                  nBonds++;
                  bondsSum += j1 + j2;
                  break;
                }
              }
            }
          }
          if (nBonds == 2) {
            if (bondsSum == 6) {
              nFCC[a1]++;
              nFCC[a2]++;
            } else {
              nHCP[a1]++;
              nHCP[a2]++;
            }
          }
        }
      }
    }
  }
  // 0: fcc (421), 1: hcp (422), 2: other
  for (int i = 0; i < nAtoms; i++) {
    if (neighborLists[i].size() == 12) {
      if (nFCC[i] == 12)
        cna[i] = 0;
      else if (nFCC[i] == 6 && nHCP[i] == 6)
        cna[i] = 1;
      else
        cna[i] = 2;
    } else
      cna[i] = 2;
  }
  return;
}

void EpiCenters::coordination(long *coordinationVal, const Matter *matter,
                              double neighborCutoff) {
  long nAtoms;
  double diffR;
  nAtoms = matter->numberOfAtoms();
  for (int i = 0; i < nAtoms; i++)
    coordinationVal[i] = 0;
  //----- Initialize end -----
  // std::cout<<"coordination\n";

  for (int i = 0; i < nAtoms - 1; i++) {
    for (int j = i + 1; j < nAtoms; j++) {
      // determine coordination number
      diffR = matter->distance(i, j);
      if (diffR < neighborCutoff) {
        coordinationVal[i] = coordinationVal[i] + 1;
        coordinationVal[j] = coordinationVal[j] + 1;
      }
    }
  }
  return;
}

void EpiCenters::coordinationLessOrEqual(bool *result, long coordinationMaxVal,
                                         const Matter *matter,
                                         double neighborCutoff) {
  long *coordinationVal;
  long nAtoms;
  nAtoms = matter->numberOfAtoms();
  coordinationVal = new long[nAtoms];
  //----- Initialize end -----
  // std::cout<<"coordinationLessOrEqual\n";

  coordination(coordinationVal, matter, neighborCutoff);

  for (int i = 0; i < nAtoms; i++) {
    if ((coordinationVal[i] < coordinationMaxVal + 1) && !(matter->getFixed(i)))
      result[i] = true;
    else
      result[i] = false;
  }
  delete[] coordinationVal;

  return;
}

long EpiCenters::minCoordination(const Matter *matter, double neighborCutoff) {
  long *coordinationVal;
  long nAtoms;
  long minCoordinationVal;
  nAtoms = matter->numberOfAtoms();
  coordinationVal = new long[nAtoms];
  //----- Initialize end -----
  // std::cout<<"minCoordination\n";

  coordination(coordinationVal, matter, neighborCutoff);
  // LONG_MAX is a the maximal value a long can achieve, from library limits
  minCoordinationVal = LONG_MAX;

  for (int i = 0; i < nAtoms; i++) {
    if ((coordinationVal[i] < minCoordinationVal) && !(matter->getFixed(i))) {
      minCoordinationVal = coordinationVal[i];
    }
  }
  delete[] coordinationVal;

  return (minCoordinationVal);
}
