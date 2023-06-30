
#include "BondBoost.h"
#include "HelperFunctions.h"
#include <math.h>
const char Hyperdynamics::NONE[] = "none";
const char Hyperdynamics::BOND_BOOST[] = "bond_boost";

BondBoost::BondBoost(Matter *matt, Parameters *params) {
  matter = matt;
  parameters = params;
  nAtoms = matter->numberOfAtoms();
  BAList = nullptr;
  RAList = nullptr;
  TABAList = nullptr;
  BBAList = nullptr;
  log = spdlog::get("combi");
}

BondBoost::~BondBoost() {
  if (BAList != nullptr) {
    delete[] BAList;
  }
  if (RAList != nullptr) {
    delete[] RAList;
  }
  if (TABAList != nullptr) {
    delete[] TABAList;
  }
  if (BBAList != nullptr) {
    delete[] BBAList;
  }
  return;
}

void BondBoost::initialize() {
  long i, j, k = 0, count = 0;
  bool flag = 1;
  string BALstring; // Boosted Atom List String
  vector<int> atoms;
  int leng_strlist;
  nBBs = 0;
  nReg = 1;
  BALstring = parameters->bondBoostBALS;
  atoms = helper_functions::split_string_int(BALstring, ",");
  leng_strlist = atoms.size();
  if (BALstring.c_str() == string("all") or atoms.size() == 0) {
    SPDLOG_LOGGER_DEBUG(log, "boost all atoms that are set free\n");
    nBAs = matter->numberOfFreeAtoms();
    nRAs = nAtoms - nBAs;    // nRestAtoms
    BAList = new long[nBAs]; // BoostAtomsList
    RAList = new long[nRAs]; // RestAtomsList
    for (i = 0; i < nAtoms; i++) {
      if (!matter->getFixed(i)) {
        BAList[k] = i;
        k++;
      }
    }
    for (i = 0; i < nAtoms; i++) {
      flag = 1;
      for (j = 0; j < nBAs; j++) {
        if (i == BAList[j]) {
          flag = 0;
        }
      }
      if (flag == 1) {
        RAList[count] = i;
        count++;
      }
    }
    if (count != nRAs) {
      SPDLOG_LOGGER_DEBUG(
          log, "Error: nRestAtoms does not equal to counted number!\n");
    }
  } else {
    SPDLOG_LOGGER_DEBUG(log, "boost the following selected atoms:");
    for (i = 0; i < leng_strlist; i++) {
      SPDLOG_LOGGER_DEBUG(log, "{} ", atoms[i]);
    }
    SPDLOG_LOGGER_DEBUG(log, "\n");
    nBAs = leng_strlist;
    nRAs = nAtoms - nBAs;    // nRestAtoms
    BAList = new long[nBAs]; // BoostAtomsList
    RAList = new long[nRAs]; // RestAtomsList
    for (i = 0; i < nBAs; i++) {
      BAList[k] = atoms[i];
      k++;
    }
    for (i = 0; i < nAtoms; i++) {
      flag = 1;
      for (j = 0; j < nBAs; j++) {
        if (i == BAList[j]) {
          flag = 0;
        }
      }
      if (flag == 1) {
        RAList[count] = i;
        count++;
      }
    }
    if (count != nRAs) {
      SPDLOG_LOGGER_DEBUG(
          log, "Error: nRestAtoms does not equal to counted number!\n");
    }
  }

  nTABs = nBAs * (nBAs - 1) / 2 +
          nBAs * nRAs;            // number of Bonds involved with Tagged Atoms
  TABAList = new long[2 * nTABs]; // CorrespondingAtomsList of TABLList
  TABLList.setZero(nTABs, 1);
  SPDLOG_LOGGER_DEBUG(log, "BondBoost Used !\n");

  // TEST PRINT
  /*
      for(i=0 ;i< nBAs;i++){
      SPDLOG_LOGGER_DEBUG(log, "nBoostAtoms {}: {}\n", i, BAList[i]);
      }

      for(i=0; i< nRAs;i++){
      SPDLOG_LOGGER_DEBUG(log, "nRestAtoms {}: {}\n", i, RAList[i]);
      }
  */
}

double BondBoost::boost() {
  long RMDS;
  double biasPot; // , AVE_Boost_Fact;
  Matrix<double, Eigen::Dynamic, 1> TABL_tmp(nTABs, 1);
  bool flag = 0;

  RMDS = int(parameters->bondBoostRMDTime / parameters->mdTimeStep);
  biasPot = 0.0;

  if (nReg <= RMDS) {
    flag = 0;
  } else {
    flag = 1;
  }

  if (flag == 0) {
    TABL_tmp = Rmdsteps();
    TABLList = TABLList + (1.0 / RMDS) * TABL_tmp;
    nReg++;
    // TEST PRINT
    /*  for(long i=nTABs-1;i<nTABs;i++){
        SPDLOG_LOGGER_DEBUG(log, "Distance between Atoms {} and {} is {}\n",
       TABAList[2 * i], TABAList[2 * i + 1], TABLList(i, 0));
            }

            SPDLOG_LOGGER_DEBUG(log, "\n");
    */
    // TEST END
  } else {
    // SPDLOG_LOGGER_DEBUG(log, "nreg = {}; RMDS={}\n", nReg, RMDS);
    if (nReg == RMDS + 1) {
      nBBs = BondSelect();

      //   for (long i=0;i<nBBs;i++){
      // SPDLOG_LOGGER_DEBUG(log, "Equilibrium Distance between Atoms {} and {}
      // is {}\n", BBAList[2 * i], BBAList[2 * i + 1], EBBLList(i, 0));
      //    }
    }

    Epsr_Q = new double[nBBs];
    CBBLList.setZero(nBBs, 1);
    biasPot = Booststeps();
    nReg++;
    delete Epsr_Q;
  }
  return biasPot;
}

double BondBoost::Booststeps() {
  long i, j; //,Mi;
  long AtomI_1, AtomI_2;
  double QRR, PRR, Epsr_MAX, A_EPS_M, Sum_V, Boost_Fact, DVMAX;
  double Dforce, Fact_1, Fact_2; //, Mforce
  double Ri[3] = {0.0}, R = 0.0;

  AtomMatrix Free(nAtoms, 3);
  AtomMatrix AddForces(nBBs, 3);
  AtomMatrix TADF(nAtoms, 3);
  AtomMatrix BiasForces(nAtoms, 3);

  TADF.setZero();
  Free.setZero();
  AddForces.setZero();
  BiasForces.setZero();

  QRR = parameters->bondBoostQRR;
  PRR = parameters->bondBoostPRR;
  DVMAX = parameters->bondBoostDVMAX;
  Epsr_MAX = 0.0;
  A_EPS_M = 0.0;
  Boost_Fact = 0.0;
  Sum_V = 0.0;
  Dforce = 0.0;
  // Mforce = 0.0;
  Fact_1 = 0.0;
  Fact_2 = 0.0;

  for (i = 0; i < nBBs; i++) {
    AtomI_1 = BBAList[2 * i];
    AtomI_2 = BBAList[2 * i + 1];
    CBBLList(i, 0) = matter->distance(AtomI_1, AtomI_2);
  }

  for (i = 0; i < nBBs; i++) {
    Epsr_Q[i] = (CBBLList(i, 0) - EBBLList(i, 0)) / EBBLList(i, 0) / QRR;
    if (abs(Epsr_Q[i]) >= Epsr_MAX) {
      Epsr_MAX = abs(Epsr_Q[i]);
      // Mi = i;
    }
  }

  A_EPS_M = (1.0 - Epsr_MAX * Epsr_MAX) * (1.0 - Epsr_MAX * Epsr_MAX) /
            (1.0 - PRR * PRR * Epsr_MAX * Epsr_MAX);
  /*
      for (i=0;i<nBBs;i++){
      SPDLOG_LOGGER_DEBUG(log, "Boost::Distance between Atoms {} and {} is {},
     EQ= {}; EPSR_Q={}\n", BBAList[2 * i], BBAList[2 * i + 1], CBBLList(i, 0),
     EBBLList(i, 0), Epsr_Q[i]);
      }

      SPDLOG_LOGGER_DEBUG(log, "Boost::Epsr_MAX= {}, Atoms {} and {}\n",
     Epsr_MAX, BBAList[2 * Mi], BBAList[2 * Mi + 1]); SPDLOG_LOGGER_DEBUG(log,
     "Boost::A_EPS_M= {}\n", A_EPS_M);
  */
  if (Epsr_MAX < 1.0) {
    for (i = 0; i < nBBs; i++) {
      Sum_V += DVMAX * (1.0 - Epsr_Q[i] * Epsr_Q[i]) / double(nBBs);
    }
  } else if (Epsr_MAX >= 1.0) {
    Sum_V = 0.0;
    A_EPS_M = 0.0;
  }

  Boost_Fact = A_EPS_M * Sum_V;
  for (i = 0; i < nBBs; i++) {
    if (abs(Epsr_Q[i]) < Epsr_MAX) {
      Fact_1 = 2.0 * A_EPS_M * DVMAX * Epsr_Q[i] / QRR / EBBLList(i, 0) /
               double(nBBs);
      Dforce = Fact_1;
    } else if (abs(Epsr_Q[i]) == Epsr_MAX) {
      Fact_1 = 2.0 * A_EPS_M * DVMAX * Epsr_Q[i] / QRR / EBBLList(i, 0) /
               double(nBBs);
      double Fact_tmp1 = (1.0 - PRR * PRR * Epsr_Q[i] * Epsr_Q[i]);
      double Fact_tmp2 = (1.0 - Epsr_Q[i] * Epsr_Q[i]);
      Fact_2 = 2.0 * Fact_tmp2 * Epsr_Q[i] *
               (2.0 * Fact_tmp1 - PRR * PRR * Fact_tmp2) / QRR /
               EBBLList(i, 0) / Fact_tmp1 / Fact_tmp1;
      Dforce = Fact_1 + Sum_V * Fact_2;
    }

    AtomI_1 = BBAList[2 * i];
    AtomI_2 = BBAList[2 * i + 1];

    R = CBBLList(i, 0);
    //    matter->distance(AtomI_1,AtomI_2);
    // SPDLOG_LOGGER_DEBUG(log, "R={}, CLList={}\n", R, CBBLList(i, 0));

    for (j = 0; j < 3; j++) {
      Ri[j] = matter->pdistance(AtomI_1, AtomI_2, j);
      AddForces(i, j) = Ri[j] / R * Dforce;
      TADF(AtomI_1, j) = TADF(AtomI_1, j) + AddForces(i, j);
      TADF(AtomI_2, j) = TADF(AtomI_2, j) - AddForces(i, j);
    }
    // SPDLOG_LOGGER_DEBUG(log, "{} Bond:: Rx={}; Ry={}; Rz={}; Rsum={};
    // R={}\n", i, Ri[0], Ri[1], Ri[2], std::sqrt(Ri[0] * Ri[0] + Ri[1] * Ri[1]
    // + Ri[2] * Ri[2]), R);
  }

  // if (i != Mi) {
  //   SPDLOG_LOGGER_DEBUG(
  //       log, "DeltaF(Epsr_Q = {}) = {}, {}, {}, {}; Atom {} and Atom {}\n",
  //       Epsr_Q[i], AddForces(i, 0), AddForces(i, 1), AddForces(i, 2),
  //       BBAList[2 * i], BBAList[2 * i + 1]);
  // } else {
  //   SPDLOG_LOGGER_DEBUG(
  //       log, "DeltaF(Epsr_MAX = {}) = {}, {}, {}, {}; Atom {} and Atom {}\n",
  //       Epsr_MAX, Mforce, AddForces(i, 0), AddForces(i, 1), AddForces(i, 2),
  //       BBAList[2 * i], BBAList[2 * i + 1]);
  // }

  // SPDLOG_LOGGER_DEBUG(log, "TADF::");
  // for (i = 0; i < nAtoms; i++) {
  //   SPDLOG_LOGGER_DEBUG(log, "{}   {}   {}", TADF(i, 0), TADF(i, 1),
  //                       TADF(i, 2));
  // }

  // SPDLOG_LOGGER_DEBUG(log, "OldForces::");
  // for (i = 0; i < nAtoms; i++) {
  //   SPDLOG_LOGGER_DEBUG(log, "{}   {}   {}", OldForce(i, 0), OldForce(i, 1),
  //                       OldForce(i, 2));
  // }

  // SPDLOG_LOGGER_DEBUG(log, "NewForces::");
  // for (i = 0; i < nAtoms; i++) {
  //   SPDLOG_LOGGER_DEBUG(log, "{}   {}   {}", NewForce(i, 0), NewForce(i, 1),
  //                       NewForce(i, 2));
  // }

  // OldForce = matter->getForces();
  // SPDLOG_LOGGER_DEBUG(log, "NewSettedForces::");
  // for (i = 0; i < nAtoms; i++) {
  //   SPDLOG_LOGGER_DEBUG(log, "{}   {}   {}", OldForce(i, 0), OldForce(i, 1),
  //                       OldForce(i, 2));
  // }

  // SPDLOG_LOGGER_DEBUG(log, "boost_fact= {}, totE= {}\n", Boost_Fact,
  // Boost_Fact + matter->getKineticEnergy() + matter->getPotentialEnergy());
  BiasForces = TADF;
  Free = matter->getFree();
  BiasForces = BiasForces.array() * Free.array();
  matter->setBiasForces(BiasForces);
  return Boost_Fact;
}

Matrix<double, Eigen::Dynamic, 1> BondBoost::Rmdsteps() {
  Matrix<double, Eigen::Dynamic, 1> TABL_t(nTABs, 1);
  long count = 0, i, j;

  for (i = 0; i < nBAs; i++) {
    for (j = i + 1; j < nBAs; j++) {
      // SPDLOG_LOGGER_DEBUG(log, "distance = {}\n", matter->distance(BAList[i],
      // BAList[j]));
      TABL_t(count, 0) = matter->distance(BAList[i], BAList[j]);
      TABAList[2 * count] = BAList[i];
      TABAList[2 * count + 1] = BAList[j];
      count++;
    }
  }

  for (i = 0; i < nBAs; i++) {
    for (j = 0; j < nRAs; j++) {
      // SPDLOG_LOGGER_DEBUG(log, "distance = {}\n", matter->distance(BAList[i],
      // RAList[j]));
      TABL_t(count, 0) = matter->distance(BAList[i], RAList[j]);
      TABAList[2 * count] = BAList[i];
      TABAList[2 * count + 1] = RAList[j];
      count++;
    }
  }
  // SPDLOG_LOGGER_DEBUG(log, "count={}\n", count);
  if (count != nTABs) {
    SPDLOG_LOGGER_DEBUG(
        log, "Total Involved Bond number does not equal counted number\n");
  }

  // SPDLOG_LOGGER_DEBUG(log, "test here");
  // for (i = count - 1; i < count; i++) {
  //   SPDLOG_LOGGER_DEBUG(log, "test ::Distance between Atoms {} and {} is
  //   {}\n",
  //                       TABAList[2 * i], TABAList[2 * i + 1], TABL_t(i, 0));
  // }
  return TABL_t;
}

long BondBoost::BondSelect() {
  long count = 0, i, nBBs_tmp = 0;
  double Qcutoff = parameters->bondBoostQcut;

  for (i = 0; i < nTABs; i++) {
    if (TABLList(i, 0) <= Qcutoff) {
      count++;
    }
  }
  nBBs_tmp = count;
  EBBLList.setZero(nBBs_tmp, 1);
  BBAList = new long[2 * nBBs_tmp];
  count = 0;
  for (i = 0; i < nTABs; i++) {
    if (TABLList(i, 0) <= Qcutoff) {
      EBBLList(count, 0) = TABLList(i, 0);
      BBAList[2 * count] = TABAList[2 * i];
      BBAList[2 * count + 1] = TABAList[2 * i + 1];
      count++;
    }
  }
  return nBBs_tmp;
}
