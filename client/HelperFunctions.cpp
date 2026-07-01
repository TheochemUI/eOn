/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "HelperFunctions.h"
#include "EonLogger.h"
#include "GeometryAnalysis.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "SafeMath.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#ifndef _WIN32
#include <sys/resource.h>
#include <sys/time.h>
#endif
using std::ifstream;
using std::string;

// Vector functions.
// Make v1 orthogonal to v2
AtomMatrix eonc::helpers::makeOrthogonal(const AtomMatrix v1,
                                         const AtomMatrix v2) {
  return v1 - matDot(v1, v2) * eonc::safemath::safe_normalized(v2);
}

void eonc::helpers::getTime(double *real, double *user, double *sys) {
  // Wall-clock time via C++11 chrono (portable)
  using namespace std::chrono;
  auto now = steady_clock::now();
  *real = duration<double>(now.time_since_epoch()).count();

#ifdef _WIN32
  if (user)
    *user = 0.0;
  if (sys)
    *sys = 0.0;
#else
  struct rusage r_usage;
  if (getrusage(RUSAGE_SELF, &r_usage) != 0) {
    EONC_LOG_WARNING("problem getting usage info: {}", strerror(errno));
  }
  if (user) {
    *user = static_cast<double>(r_usage.ru_utime.tv_sec) +
            static_cast<double>(r_usage.ru_utime.tv_usec) / 1e6;
  }
  if (sys) {
    *sys = static_cast<double>(r_usage.ru_stime.tv_sec) +
           static_cast<double>(r_usage.ru_stime.tv_usec) / 1e6;
  }
#endif
}

bool eonc::helpers::existsFile(string filename) {
  return std::filesystem::exists(filename);
}

string eonc::helpers::getRelevantFile(string filename) {
  string filenameRelevant;
  string filenamePrefix;
  string filenamePostfix;

  // check if the _cp version of the file is present
  int i = filename.rfind(".");
  filenamePrefix.assign(filename, 0, i);
  filenamePostfix.assign(filename, i, filename.size());
  filenameRelevant = filenamePrefix + "_cp" + filenamePostfix;
  if (existsFile(filenameRelevant)) {
    return filenameRelevant;
  }
  // check if the _in version of the file is present
  filenameRelevant = filenamePrefix + "_in" + filenamePostfix;
  if (existsFile(filenameRelevant)) {
    return filenameRelevant;
  }
  // otherwise return original filename
  return filename;
}

VectorXd eonc::helpers::loadMasses(string filename, int nAtoms) {
  ifstream massFile(filename.c_str());
  if (!massFile.is_open()) {
    EONC_LOG_CRITICAL("File {} was not found", filename);
    std::exit(1);
  }

  VectorXd masses(nAtoms);
  for (int i = 0; i < nAtoms; i++) {
    double mass;
    if (!(massFile >> mass)) {
      EONC_LOG_CRITICAL("Error reading {}", filename);
      std::exit(1);
    }
    masses(i) = mass;
  }

  massFile.close();

  return masses;
}

AtomMatrix eonc::helpers::loadMode(FILE *modeFile, int nAtoms) {
  AtomMatrix mode;
  mode.resize(nAtoms, 3);
  mode.setZero();
  for (int i = 0; i < nAtoms; i++) {
    fscanf(modeFile, "%lf %lf %lf", &mode(i, 0), &mode(i, 1), &mode(i, 2));
  }
  return mode;
}

AtomMatrix eonc::helpers::loadMode(string filename, int nAtoms) {
  // Unique FILE* with RAII cleanup
  auto closer = [](FILE *f) {
    if (f)
      std::fclose(f);
  };
  std::unique_ptr<FILE, decltype(closer)> modeFile(
      std::fopen(filename.c_str(), "rb"), closer);
  if (!modeFile) {
    EONC_LOG_CRITICAL("File {} was not found\n Stopping", filename);
    std::exit(1);
  }
  return loadMode(modeFile.get(), nAtoms);
}

bool eonc::helpers::loadOrSynthesizeDisplacement(
    Matter &target, const Matter &initial, const std::string &displacementPath,
    const std::string &modePath, double scale) {
  if (eonc::io::io_ok(target.con2matter(displacementPath))) {
    // displacement.con may carry stale fixed-atom coordinates from a prior run.
    const AtomMatrix &initPos = initial.getPositions();
    AtomMatrix pos = target.getPositionsCopy();
    const long n = initial.numberOfAtoms();
    for (long i = 0; i < n; i++) {
      if (initial.getFixed(i)) {
        pos.row(i) = initPos.row(i);
      }
    }
    target.setPositions(pos);
    return true;
  }
  if (!existsFile(modePath)) {
    return false;
  }
  AtomMatrix mode =
      loadMode(modePath, static_cast<int>(initial.numberOfAtoms()));
  const double norm = mode.norm();
  if (!(norm > 0.0)) {
    return false;
  }
  mode *= (scale / norm);
  target = initial;
  AtomMatrix pos = initial.getPositionsCopy();
  pos += mode;
  const AtomMatrix &initPos = initial.getPositions();
  const long n = initial.numberOfAtoms();
  for (long i = 0; i < n; i++) {
    if (initial.getFixed(i)) {
      pos.row(i) = initPos.row(i);
    }
  }
  target.setPositions(pos);
  EONC_LOG_INFO("Synthesized displacement from pos.con + scale {:.6g} * unit "
                "mode in {} (missing {})",
                scale, modePath, displacementPath);
  return true;
}

void eonc::helpers::saveMode(FILE *modeFile, std::shared_ptr<Matter> matter,
                             AtomMatrix mode) {
  long const nAtoms = matter->numberOfAtoms();
  for (long i = 0; i < nAtoms; ++i) {
    if (matter->getFixed(i)) {
      fprintf(modeFile, "0 0 0\n");
    } else {
      fprintf(modeFile, "%lf\t%lf \t%lf\n", mode(i, 0), mode(i, 1), mode(i, 2));
    }
  }
  return;
}

void eonc::helpers::saveMode(const std::string &filename,
                             std::shared_ptr<Matter> matter, AtomMatrix mode) {
  std::ofstream out(filename);
  if (!out)
    return;
  long const nAtoms = matter->numberOfAtoms();
  for (long i = 0; i < nAtoms; ++i) {
    if (matter->getFixed(i)) {
      out << "0 0 0\n";
    } else {
      out << std::format("{:f}\t{:f} \t{:f}\n", mode(i, 0), mode(i, 1),
                         mode(i, 2));
    }
  }
}

std::vector<int> eonc::helpers::split_string_int(std::string s,
                                                 std::string delim) {
  std::vector<int> list;
  if (s.empty())
    return list;

  size_t start = 0;
  size_t end = s.find_first_of(delim);
  while (start < s.size()) {
    auto token = s.substr(start, end - start);
    if (!token.empty()) {
      try {
        list.push_back(std::stoi(token));
      } catch (const std::exception &) {
        return {}; // Parse error
      }
    }
    if (end == std::string::npos)
      break;
    start = end + 1;
    end = s.find_first_of(delim, start);
  }
  return list;
}

namespace {
class MatterObjectiveFunction : public ObjectiveFunction {
  Matter &m_matter; // non-owning reference, avoids copy
public:
  MatterObjectiveFunction(Matter &mat, const Parameters &parametersPassed)
      : ObjectiveFunction(parametersPassed),
        m_matter{mat} {}
  ~MatterObjectiveFunction() = default;
  double getEnergy() { return m_matter.getPotentialEnergy(); }
  VectorXd getGradient(bool fdstep = false) {
    return -m_matter.getForcesFreeV();
  }
  void setPositions(const VectorXd &x) { m_matter.setPositionsFreeV(x); }
  VectorXd getPositions() { return m_matter.getPositionsFreeV(); }
  int degreesOfFreedom() { return 3 * m_matter.numberOfFreeAtoms(); }
  bool isConverged() {
    return getConvergence() < params.optimizer_options.converged_force;
  }
  double getConvergence() {
    if (params.optimizer_options.convergence_metric == "norm") {
      return m_matter.getForcesFreeV().norm();
    } else if (params.optimizer_options.convergence_metric == "max_atom") {
      return m_matter.maxForce();
    } else if (params.optimizer_options.convergence_metric == "max_component") {
      return m_matter.getForces().maxCoeff();
    } else {
      EONC_LOG_CRITICAL("{} Unknown opt_convergence_metric: {}", "[Matter]",
                        params.optimizer_options.convergence_metric);
      std::exit(1);
    }
  }
  VectorXd difference(const VectorXd &a, const VectorXd &b) {
    return m_matter.pbcV(a - b);
  }
};
} // namespace

bool eonc::helpers::relaxMatter(Matter &matter, const Parameters &params,
                                bool quiet, bool writeMovie, bool checkpoint,
                                std::string prefixMovie,
                                std::string prefixCheckpoint) {
  eonc::log::Scoped m_log;
  auto objf = std::make_shared<MatterObjectiveFunction>(matter, params);
  auto optim = eonc::helpers::create::mkOptim(
      objf, params.optimizer_options.method, params);

  std::ostringstream min;
  min << prefixMovie;
  std::string minDatFilename = prefixMovie + ".dat";
  auto write_movie_frame = [&](uint64_t frameIndex, bool append,
                               double stepSize) {
    eonc::io::ConFrameMetadata metadata;
    metadata.frame_index = frameIndex;
    metadata.energy = matter.getPotentialEnergy();
    metadata.scalars.push_back({"step_size", stepSize});
    metadata.scalars.push_back({"convergence", objf->getConvergence()});
    matter.matter2con(min.str(), append, &metadata);

    if (params.debug_options.write_deprecated_outs) {
      std::ofstream minDat(minDatFilename,
                           append ? (std::ios::binary | std::ios::app)
                                  : std::ios::binary);
      if (minDat) {
        if (!append) {
          minDat << "iteration\tstep_size\tconvergence\tenergy\n";
        }
        minDat << std::format("{}\t{:.5e}\t{:.5e}\t{:.6f}\n", frameIndex,
                              stepSize, objf->getConvergence(),
                              matter.getPotentialEnergy());
      }
    }
  };
  if (writeMovie) {
    write_movie_frame(0, false, 0.0);
  }

  int iteration = 0;
  if (!quiet) {
    QUILL_LOG_DEBUG(m_log, "{} {:10s}  {:14s}  {:18s}  {:13s}\n", "[Matter]",
                    "Iter", "Step size",
                    params.optimizer_options.convergence_metric_label,
                    "Energy");
    QUILL_LOG_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}\n",
                    "[Matter]", iteration, 0.0, objf->getConvergence(),
                    matter.getPotentialEnergy());
  }

  while (!objf->isConverged() &&
         iteration < params.optimizer_options.max_iterations) {

    AtomMatrix pos = matter.getPositions();

    optim->step(params.optimizer_options.max_move);
    iteration++;

    double stepSize =
        eonc::geometry::maxAtomMotion(matter.pbc(matter.getPositions() - pos));

    if (!quiet) {
      QUILL_LOG_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                      "[Matter]", iteration, stepSize, objf->getConvergence(),
                      matter.getPotentialEnergy());
    }

    if (writeMovie) {
      write_movie_frame(static_cast<uint64_t>(iteration), true, stepSize);
    }

    if (checkpoint) {
      std::ostringstream chk;
      chk << prefixCheckpoint << "_cp";
      matter.matter2con(chk.str(), false);
    }
  }

  if (iteration == 0) {
    if (!quiet) {
      QUILL_LOG_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                      "[Matter]", iteration, 0.0, objf->getConvergence(),
                      matter.getPotentialEnergy());
    }
  }
  return objf->isConverged();
}
