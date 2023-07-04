#include "GlobalOptimizationJob.h"
#include "GlobalOptimization.h"
// #include "MinimizationJob.h"
#include "Dynamics.h"
#include "HelperFunctions.h"
#include <spdlog/spdlog.h>

using namespace helper_functions;

std::vector<std::string> GlobalOptimizationJob::run(void) {
  // int status;
  GlobalOptimization globopt = GlobalOptimization(params.get());
  string reactant_passed =
      helper_functions::getRelevantFile(params->conFilename);
  std::vector<std::string> returnFiles;
  // returnFiles.push_back(reactant_output);
  Matter *matter_cur = new Matter(pot, params);
  Matter *matter_hop = new Matter(pot, params);
  matter_cur->con2matter(reactant_passed);
  bool converged;
  long nstep = params->globalOptimizationSteps;
  AtomMatrix rat_t(matter_cur->numberOfAtoms(), 3);
  // double epot_hop;
  // std::vector<double> earr;
  // std::vector<Matter> allmatter;
  SPDLOG_LOGGER_DEBUG(log, "\nBeginning minima hopping of %s",
                      reactant_passed.c_str());
  // long fcalls;
  // SPDLOG_LOGGER_DEBUG(log, "fcalls= {}", matter_cur->getForceCalls());
  // SPDLOG_LOGGER_DEBUG(log, "epot= {:24.15E}",
  // matter_cur->getPotentialEnergy()); SPDLOG_LOGGER_DEBUG(log, "fcalls= {}",
  // matter_cur->getForceCalls()); SPDLOG_LOGGER_DEBUG(log, "epot= {:24.15E}",
  // matter_cur->getPotentialEnergy()); SPDLOG_LOGGER_DEBUG(log, "fcalls= {}",
  // matter_cur->getForceCalls());
  converged = matter_cur->relax(false, params->writeMovies, params->checkpoint,
                                "min", "matter_cur");
  SPDLOG_LOGGER_DEBUG(log, "converged {}", (converged) ? "TRUE" : "FALSE");
  // nlmin=0;
  // if(nlmin==0)
  earr.push_back(matter_cur->getPotentialEnergy());
  matter_hop[0] = matter_cur[0];
  GlobalOptimizationJob::report(matter_hop);
  for (long istep = 1; istep <= nstep; istep++) {
    // hoppingStep attempts to hop into a new state followed by a minimization
    hoppingStep(istep, matter_cur, matter_hop);
    // decisionStep decides to accept or reject this step if it escaped
    decisionStep(matter_cur, matter_hop);
    // analyzing what happened in this step if it escaped
    analyze(matter_cur, matter_hop);
    // reporting useful information about this hop
    report(matter_hop);
    if (matter_cur->getPotentialEnergy() <
        params->globalOptimizationTargetEnergy)
      break;
  }
  for (size_t i = 0; i < earr.size(); i++) {
    double earrim1;
    if (i == 0) {
      earrim1 = earr[0];
    } else {
      earrim1 = earr[i - 1];
    }
    fprintf(earrfile, "%5lu  %15.5f  %15.5f  %15.5f  ", i + 1, earr[i],
            earr[i] - earr[0], earr[i] - earrim1);
  }
  // globopt.run();
  return returnFiles;
} // end of GlobalOptimizationJob::run

void GlobalOptimizationJob::analyze(Matter *matter_cur, Matter *matter_hop) {
  if (escapeResult == "failure") {
    // matter_hop[0]=matter_cur[0];
    hoppingResult = "same";
    return;
  }
  double epot = matter_hop->getPotentialEnergy();
  size_t jlo = hunt(epot);
  // SPDLOG_LOGGER_DEBUG(log, "REZA: {}", jlo);
  if (abs(epot - earr[jlo]) < params->energyDifference) {
    hoppingResult = "already_visited";
  } else {
    hoppingResult = "new";
  }
  if (decisionResult == "accepted") {
    double epot_hop = matter_hop->getPotentialEnergy();
    matter_cur[0] = matter_hop[0];
    size_t jlo = hunt(matter_hop->getPotentialEnergy());
    if (hoppingResult == "new" && jlo == 0 && epot_hop < earr[0]) {
      SPDLOG_LOGGER_DEBUG(
          log, "new lowest: nlmin, epot_hop, dE {:7d} {:15.5f} {:10.5f}", 1,
          epot_hop, epot_hop - earr[0]);
    }
    // log = spdlog::get("_traceback");
    // SPDLOG_LOGGER_CRITICAL(log, "{:15.5f}  {:15.5f}  {:15.5f}  ", epot_hop,
    // earr[0], earr[1]); std::exit(1);
    insert(matter_cur);
  } else if (decisionResult == "rejected") {
    // matter_hop[0]=matter_cur[0];
  } else {
    log = spdlog::get("_traceback");
    SPDLOG_LOGGER_CRITICAL(
        log,
        "ERROR: new minimum is neither accepted nor rejected: client stops.");
    std::exit(1);
  }
}

void GlobalOptimizationJob::examineEscape(Matter *matter_cur,
                                          Matter *matter_hop) {
  // fSPDLOG_LOGGER_DEBUG(log, monfile,"%15.5f  %15.5f  %15.5f
  // ",matter_hop->getPotentialEnergy(),
  // matter_cur->getPotentialEnergy(),matter_hop->getPotentialEnergy()-matter_cur->getPotentialEnergy());
  double epot, epot_hop;
  epot = matter_cur->getPotentialEnergy();
  epot_hop = matter_hop->getPotentialEnergy();
  if (abs(epot_hop - epot) < params->energyDifference) {
    escapeResult = "failure";
  } else {
    escapeResult = "success";
    // fprintf(monfile,"%15.5f  %15.5f  %15.5f           ", epot_hop, ediff,
    // ekin);
  }
}

void GlobalOptimizationJob::applyMoveFeedbackMD(void) {
  if (firstStep) {
    firstStep = false;
    return;
  }
  if (hoppingResult == "same") {
    ekin *= beta1;
  } else if (hoppingResult == "already_visited") {
    ekin *= beta2;
  } else if (hoppingResult == "new") {
    ekin *= beta3;
  } else {
    log = spdlog::get("_traceback");
    SPDLOG_LOGGER_CRITICAL(log,
                           "ERROR: client does not know what to do with ekin.");
    SPDLOG_LOGGER_CRITICAL(log, "ERROR: client stops in applyMoveFeedbackMD.");
    std::exit(1);
  }
}

void GlobalOptimizationJob::applyDecisionFeedback(void) {
  if (decisionResult == "accepted") {
    ediff *= alpha1;
  } else {
    ediff *= alpha2;
  }
}

void GlobalOptimizationJob::report(Matter *matter_hop) {
  char C1, C2;
  double ekin_p;
  if (hoppingResult == "same") {
    C1 = 'S';
    C2 = '-';
    ekin_p = ekin * beta1;
  } else if (hoppingResult == "already_visited") {
    C1 = 'O';
    ekin_p = ekin * beta2;
  } else if (hoppingResult == "new") {
    C1 = 'N';
    ekin_p = ekin * beta3;
  } else {
    C1 = '-';
    ekin_p = ekin;
  }
  if (decisionResult == "accepted") {
    C2 = 'A';
  } else if (decisionResult == "rejected") {
    C2 = 'R';
  } else {
    C2 = '-';
  }
  double epot_hop = matter_hop->getPotentialEnergy();
  double temp = (2.0 * ekin_p / params->kB);
  double dt = params->mdTimeStep;
  fprintf(monfile, "%15.5f  %15.5f  %11lu  %12.2f         %c%c  %5ld  %5ld",
          epot_hop, ediff, (size_t)temp, dt, C1, C2, fcallsMove, fcallsRelax);
}

void GlobalOptimizationJob::decisionStep(Matter *matter_cur,
                                         Matter *matter_hop) {
  decisionResult = "unknown";
  examineEscape(matter_cur, matter_hop);
  if (escapeResult == "failure") {
    // matter_hop[0] = matter_cur[0];
    return;
  }
  if (params->globalOptimizationDecisionMethod == "npew") {
    acceptRejectNPEW(matter_cur, matter_hop);
    // GlobalOptimizationJob::update_minhop_param(matter_hop);
  } else if (params->globalOptimizationDecisionMethod == "boltzmann") {
    acceptRejectBoltzmann(matter_cur, matter_hop);
  } else {
    log = spdlog::get("_traceback");
    SPDLOG_LOGGER_CRITICAL(
        log, "ERROR: accept/reject method not specified. client stops.");
    std::exit(1);
  }
  applyDecisionFeedback();
}

// NPEW: non-probablistic energy window (used in minima hopping method)
void GlobalOptimizationJob::acceptRejectNPEW(Matter *matter_cur,
                                             Matter *matter_hop) {
  if (matter_hop->getPotentialEnergy() <
      matter_cur->getPotentialEnergy() + ediff) {
    decisionResult = "accepted";
  } else {
    decisionResult = "rejected";
  }
}

void GlobalOptimizationJob::acceptRejectBoltzmann(Matter *matter_cur,
                                                  Matter *matter_hop) {
  double eTrial = matter_hop->getPotentialEnergy();
  double eCurrent = matter_hop->getPotentialEnergy();

  double deltaE = eTrial - eCurrent;
  double kB = 8.6173324e-5;

  double p;
  if (deltaE <= 0.0) {
    p = 1.0;
  } else {
    p = exp(-deltaE / params->temperature * kB);
  }

  if (randomDouble(1.0) < p) {
    decisionResult = "accepted";
  } else {
    decisionResult = "rejected";
  }
}

void GlobalOptimizationJob::hoppingStep(long istep, Matter *matter_cur,
                                        Matter *matter_hop) {
  bool converged;
  matter_hop[0] = matter_cur[0];
  long fcalls1 = matter_hop->getForceCalls();
  if (params->globalOptimizationMoveMethod == "md") {
    applyMoveFeedbackMD();
    mdescape(matter_hop);
  } else if (params->globalOptimizationMoveMethod == "random") {
    randomMove(matter_hop);
  }
  long fcalls2 = matter_hop->getForceCalls();
  hoppingResult = "unknown";
  converged = matter_hop->relax(true, params->writeMovies, params->checkpoint,
                                "min", "matter_hop");
  SPDLOG_LOGGER_DEBUG(log, "converged {}", (converged) ? "TRUE" : "FALSE");
  long fcalls3 = matter_hop->getForceCalls();
  fcallsMove = fcalls2 - fcalls1;
  fcallsRelax = fcalls3 - fcalls2;
}

void GlobalOptimizationJob::randomMove(Matter *matter) {
  // create a random displacement
  AtomMatrix displacement;
  displacement.resize(matter->numberOfAtoms(), 3);
  displacement.setZero();
  int num = matter->numberOfAtoms();

  for (int i = 0; i < num; i++) {
    double disp = params->basinHoppingDisplacement;
    if (!matter->getFixed(i)) {
      for (int j = 0; j < 3; j++) {
        if (params->basinHoppingDisplacementDistribution == "uniform") {
          displacement(i, j) = randomDouble(2 * disp) - disp;
        } else if (params->basinHoppingDisplacementDistribution == "gaussian") {
          displacement(i, j) = gaussRandom(0.0, disp);
        } else {
          log = spdlog::get("_traceback");
          SPDLOG_LOGGER_CRITICAL(log, "Unknown displacement_distribution");
          std::exit(1);
        }
      }
    }
  }
  matter->setPositions(matter->getPositions() + displacement);
}

void GlobalOptimizationJob::mdescape(Matter *matter) {
  int nmd;
  double ekinc, epot, etot, epot0, etot0;
  Dynamics *dyn = new Dynamics(matter, params.get());
  velopt(matter);
  epot = matter->getPotentialEnergy();
  ekinc = matter->getKineticEnergy();
  etot = ekinc + epot;
  epot0 = epot;
  etot0 = etot;
  size_t nummax = 0, nummin = 0;
  double enmin1 = 0.0, enmin2 = 0.0, en0000 = 0.0;
  double econs_max = -1.E100, econs_min = 1.E100, devcon;
  bool md_presumably_escaped = false;
  SPDLOG_LOGGER_DEBUG(log, "MD  {:5d}  {:20.10E}  {:15.5E}  {:15.5E}  ", 0,
                      epot - epot0, ekinc, etot - etot0);
  nmd = 1000;
  for (int imd = 1; imd <= nmd; imd++) {
    enmin2 = enmin1;
    enmin1 = en0000;
    dyn->velocityVerlet();
    epot = matter->getPotentialEnergy();
    ekinc = matter->getKineticEnergy();
    etot = ekinc + epot;
    en0000 = epot - epot0;
    if (enmin1 > enmin2 && enmin1 > en0000)
      nummax += 1;
    if (enmin1 < enmin2 && enmin1 < en0000)
      nummin += 1;
    // SPDLOG_LOGGER_DEBUG(log, "MD  {:5d}  {:15.5f}  {:15.5f}  {:12.2E}  {:4lu}
    // {:4lu}", imd, epot - epot0, ekinc, etot - etot0, nummax, nummin);
    econs_max = max(econs_max, ekinc + epot);
    econs_min = min(econs_min, ekinc + epot);
    if (nummin >= (size_t)mdmin) {
      if (nummax != nummin)
        SPDLOG_LOGGER_WARN(log, "WARNING: iproc,nummin,nummax {:4lu} {:4lu}",
                           nummin, nummax);
      md_presumably_escaped = true;
      break;
    }
  } // end of loop over imd
  devcon = econs_max - econs_min;
  if (md_presumably_escaped) {
    devcon = devcon / (double)(matter->numberOfFreeAtoms() * 3);
    if (devcon / ekin < 2.E-3) {
      params->mdTimeStep *= 1.1;
    } else {
      params->mdTimeStep /= 1.1;
    }
  } else {
    SPDLOG_LOGGER_DEBUG(log, "TOO MANY MD STEPS  ");
    params->mdTimeStep *= 2.0;
  }
}

void GlobalOptimizationJob::velopt(Matter *matter) {
  AtomMatrix vat(matter->numberOfAtoms(), 3);
  double tt1, tt2, tt3, vtot[3]; //, ekin_t;
  int iat;
  // matter->numberOfAtoms();
  vtot[0] = 0.0;
  vtot[1] = 0.0;
  vtot[2] = 0.0;
  for (iat = 0; iat < matter->numberOfAtoms(); iat++) {
    tt1 = randomDouble();
    tt2 = randomDouble();
    tt3 = randomDouble();
    vat(iat, 0) = (tt1 - 0.5) * 2.0;
    vat(iat, 1) = (tt2 - 0.5) * 2.0;
    vat(iat, 2) = (tt3 - 0.5) * 2.0;
    vtot[0] += vat(iat, 0);
    vtot[1] += vat(iat, 1);
    vtot[2] += vat(iat, 2);
  }
  SPDLOG_LOGGER_DEBUG(log, "Linear momentum  {:15.5E}  {:15.5E}  {:15.5E}  ",
                      vtot[0], vtot[1], vtot[2]);
  vtot[0] /= matter->numberOfAtoms();
  vtot[1] /= matter->numberOfAtoms();
  vtot[2] /= matter->numberOfAtoms();
  for (iat = 0; iat < matter->numberOfAtoms(); iat++) {
    vat(iat, 0) -= vtot[0];
    vat(iat, 1) -= vtot[1];
    vat(iat, 2) -= vtot[2];
  }
  matter->setVelocities(vat);
  long nFreeCoords = matter->numberOfFreeAtoms() * 3;
  double kinE = matter->getKineticEnergy();
  double kB = params->kB;
  double kinT = (2.0 * kinE / nFreeCoords / kB);
  double temperature = (2.0 * ekin / kB);
  matter->setVelocities(vat * sqrt(temperature / kinT));
}

/*
void ::rescaleVelocity()
{
    AtomMatrix velocity = matter->getVelocities();
    double kinE = matter->getKineticEnergy();
    double kinT = (2.0*kinE/nFreeCoords/kb);
    matter->setVelocities(velocity*sqrt(temperature/kinT));
}
*/

void GlobalOptimizationJob::insert(Matter *matter) {
  double epot;
  // vector<double> epot_hop;
  size_t jlo, jlo_insert;
  // vector<double>::iterator it;
  epot = matter->getPotentialEnergy();
  jlo = hunt(epot);
  SPDLOG_LOGGER_DEBUG(log, "JLO= {:5lu}  {:10.5f}  ", jlo,
                      std::abs(epot - earr[jlo]));
  // it=earr.begin()+jlo;
  // epot_hop.push_back(epot);
  if (!(abs(epot - earr[jlo]) < params->energyDifference)) {
    // earr.insert(it,epot_hop.begin(),epot_hop.end());
    jlo_insert = jlo;
    if (epot > earr[jlo])
      jlo_insert++;
    earr.insert(earr.begin() + jlo_insert, 1, epot);
  }
}

size_t GlobalOptimizationJob::hunt(double epot) {
  // epot is in interval [earr(jlo),earr(jlo+1)[ ; earr(0)=-Infinity ; earr(n+1)
  // = Infinity
  size_t jlo;
  double de;
  for (jlo = 0; jlo < earr.size(); jlo++)
    if (epot < earr[jlo])
      break;
  if (jlo == earr.size())
    jlo--;
  de = abs(epot - earr[jlo]);
  if (jlo > 0)
    if (abs(epot - earr[jlo - 1]) < de)
      jlo--; //{jlo--;de=abs(epot-earr[jlo]);}
  // if(jlo!=earr.size()-1) if(abs(epot-earr[jlo+1])<params->energyDifference)
  // jlo++;
  return jlo;
}
