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
#include "eon/HelperFunctions.h"
#include "eon/EonLogger.h"
#include "eon/EpiCenters.h"
#include "eon/GeometryAnalysis.h"
#include "eon/ObjectiveFunction.h"
#include "eon/Optimizer.h"
#include "eon/Parameters.h"
#include "eon/SafeMath.h"

#include <cassert>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#ifndef _WIN32
#include <sys/resource.h>
#include <sys/time.h>
#endif
// Vector functions.
// Make v1 orthogonal to v2
AtomMatrix eonc::helpers::makeOrthogonal(const AtomMatrix v1,
                                         const AtomMatrix v2) {
  return v1 - matDot(v1, v2) * eonc::safemath::safe_normalized(v2);
}

void eonc::helpers::getTime(double *real, double *user, double *sys) {
  using namespace std::chrono;
  auto now = steady_clock::now();
  if (real) {
    *real = duration<double>(now.time_since_epoch()).count();
  }

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

bool eonc::helpers::existsFile(std::string filename) {
  return std::filesystem::exists(filename);
}

std::string eonc::helpers::getRelevantFile(std::string filename) {
  const auto dot = filename.rfind('.');
  const std::string prefix =
      (dot == std::string::npos) ? filename : filename.substr(0, dot);
  const std::string postfix =
      (dot == std::string::npos) ? std::string{} : filename.substr(dot);
  std::string filenameRelevant = prefix + "_cp" + postfix;
  if (existsFile(filenameRelevant)) {
    return filenameRelevant;
  }
  filenameRelevant = prefix + "_in" + postfix;
  if (existsFile(filenameRelevant)) {
    return filenameRelevant;
  }
  return filename;
}

VectorXd eonc::helpers::loadMasses(std::string filename, int nAtoms) {
  std::ifstream massFile(filename.c_str());
  if (!massFile.is_open()) {
    EONC_LOG_CRITICAL("File {} was not found", filename);
    throw std::runtime_error(std::format("cannot open {}", filename));
  }

  VectorXd masses(nAtoms);
  for (int i = 0; i < nAtoms; i++) {
    double mass;
    if (!(massFile >> mass)) {
      EONC_LOG_CRITICAL("Error reading {}", filename);
      throw std::runtime_error(
          std::format("{} ended after {} of {} masses", filename, i, nAtoms));
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
    if (fscanf(modeFile, "%lf %lf %lf", &mode(i, 0), &mode(i, 1),
               &mode(i, 2)) != 3) {
      EONC_LOG_CRITICAL("Mode file ended after {} of {} atoms", i, nAtoms);
      throw std::runtime_error(
          std::format("mode file ended after {} of {} atoms", i, nAtoms));
    }
  }
  return mode;
}

AtomMatrix eonc::helpers::loadMode(std::string filename, int nAtoms) {
  // Unique FILE* with RAII cleanup
  auto closer = [](FILE *f) {
    if (f)
      std::fclose(f);
  };
  std::unique_ptr<FILE, decltype(closer)> modeFile(
      std::fopen(filename.c_str(), "rb"), closer);
  if (!modeFile) {
    EONC_LOG_CRITICAL("File {} was not found", filename);
    throw std::runtime_error(std::format("cannot open {}", filename));
  }
  return loadMode(modeFile.get(), nAtoms);
}

bool eonc::helpers::loadOrSynthesizeDisplacement(
    Matter &target, const Matter &initial, const std::string &displacementPath,
    const std::string &modePath, double scale) {
  if (eonc::io::io_ok(target.con2matter(displacementPath))) {
    if (target.numberOfAtoms() != initial.numberOfAtoms()) {
      EONC_LOG_ERROR("{} holds {} atoms, the initial structure has {}",
                     displacementPath, target.numberOfAtoms(),
                     initial.numberOfAtoms());
      return false;
    }
    // displacement.con may carry stale fixed-atom coordinates from a prior run.
    // It also usually has sequential column-5 ids; keep the reactant's.
    const AtomMatrix &initPos = initial.getPositions();
    AtomMatrix pos = target.getPositionsCopy();
    const long n = initial.numberOfAtoms();
    std::vector<long> fileMap(static_cast<size_t>(n));
    for (long i = 0; i < n; i++) {
      if (initial.getFixed(i)) {
        pos.row(i) = initPos.row(i);
      }
      target.setAtomIndex(i, initial.getAtomIndex(i));
      fileMap[static_cast<size_t>(i)] = initial.mapFileRow(i);
    }
    target.setFileToMatter(std::move(fileMap));
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

bool eonc::helpers::applyClientDisplacement(Matter &target,
                                            const Matter &initial,
                                            const Parameters &params,
                                            AtomMatrix *modeOut) {
  using namespace eonc::EpiCenters;
  const auto &opt = params.saddle_search_options();
  const std::string &dtype = opt.displace_type;
  if (dtype == DISP_LOAD) {
    return false;
  }

  long epicenter = -1;
  const double cutoff = params.structure_comparison_options().neighbor_cutoff;
  if (dtype == DISP_LISTED_ATOMS) {
    epicenter = listedAtomEpiCenter(&initial, opt.displace_atom_list);
  } else if (dtype == DISP_RANDOM) {
    epicenter = randomFreeAtomEpiCenter(&initial);
  } else if (dtype == DISP_LAST_ATOM) {
    epicenter = lastAtom(&initial);
  } else if (dtype == DISP_MIN_COORDINATED) {
    epicenter = minCoordinatedEpiCenter(&initial, cutoff);
  } else if (dtype == DISP_NOT_FCC_OR_HCP) {
    epicenter = cnaEpiCenter(&initial, cutoff);
  } else {
    return false;
  }

  target = initial;
  const long n = initial.numberOfAtoms();
  const double radius = opt.displace_radius;
  const double mag = opt.displace_magnitude;
  AtomMatrix pos = initial.getPositionsCopy();
  AtomMatrix mode = AtomMatrix::Zero(n, 3);
  for (long i = 0; i < n; ++i) {
    if (initial.getFixed(i)) {
      continue;
    }
    const double dist = (i == epicenter) ? 0.0 : initial.distance(epicenter, i);
    if (dist <= radius) {
      for (int a = 0; a < 3; ++a) {
        mode(i, a) = eonc::rng::gaussRandom(0.0, mag);
      }
    }
  }
  const double norm = mode.norm();
  if (norm > 0.0) {
    pos += mode;
    mode /= norm;
  } else if (epicenter >= 0 && epicenter < n && !initial.getFixed(epicenter)) {
    mode(epicenter, 0) = 1.0;
    pos(epicenter, 0) += mag;
  }
  target.setPositions(pos);
  if (modeOut != nullptr) {
    *modeOut = std::move(mode);
  }
  return true;
}

void eonc::helpers::saveMode(FILE *modeFile, std::shared_ptr<Matter> matter,
                             AtomMatrix mode) {
  const AtomMatrix free = matter->getFree();
  long const nAtoms = matter->numberOfAtoms();
  for (long i = 0; i < nAtoms; ++i) {
    fprintf(modeFile, "%.17g\t%.17g\t%.17g\n", free(i, 0) * mode(i, 0),
            free(i, 1) * mode(i, 1), free(i, 2) * mode(i, 2));
  }
  return;
}

void eonc::helpers::saveMode(const std::string &filename,
                             std::shared_ptr<Matter> matter, AtomMatrix mode) {
  std::ofstream out(filename);
  if (!out)
    return;
  const AtomMatrix free = matter->getFree();
  long const nAtoms = matter->numberOfAtoms();
  for (long i = 0; i < nAtoms; ++i) {
    out << std::format("{:.17g}\t{:.17g}\t{:.17g}\n", free(i, 0) * mode(i, 0),
                       free(i, 1) * mode(i, 1), free(i, 2) * mode(i, 2));
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

std::optional<std::string_view>
eonc::helpers::convergenceMetricLabel(std::string_view metric) {
  if (metric == "max_atom") {
    return "Max atom force";
  }
  if (metric == "max_component") {
    return "Max force comp";
  }
  if (metric == "norm") {
    return "||Force||";
  }
  if (metric == "rms") {
    return "RMS force";
  }
  return std::nullopt;
}

void eonc::helpers::requireKnownConvergenceMetric(std::string_view metric,
                                                  std::string_view context) {
  if (convergenceMetricLabel(metric)) {
    return;
  }
  throw std::invalid_argument(
      std::format("{} unknown convergence_metric: {}", context, metric));
}

namespace {
class MatterObjectiveFunction : public eonc::ObjectiveFunction {
  eonc::Matter &m_matter; // non-owning reference, avoids copy
public:
  MatterObjectiveFunction(eonc::Matter &mat,
                          const eonc::Parameters &parametersPassed)
      : eonc::ObjectiveFunction(parametersPassed),
        m_matter{mat} {
    eonc::helpers::requireKnownConvergenceMetric(
        params.optimizer_options().convergence_metric, "[Matter]");
  }
  ~MatterObjectiveFunction() = default;
  double getEnergy() { return m_matter.getPotentialEnergy(); }
  VectorXd getGradient(bool fdstep = false) {
    return -m_matter.getForcesFreeV();
  }
  void setPositions(const VectorXd &x) { m_matter.setPositionsFreeV(x); }
  VectorXd getPositions() { return m_matter.getPositionsFreeV(); }
  int degreesOfFreedom() { return 3 * m_matter.numberOfFreeAtoms(); }
  bool isConverged() {
    return getConvergence() < params.optimizer_options().converged_force;
  }
  double getConvergence() {
    if (params.optimizer_options().convergence_metric == "norm") {
      return m_matter.getForcesFreeV().norm();
    } else if (params.optimizer_options().convergence_metric == "rms") {
      const VectorXd f = m_matter.getForcesFreeV();
      const auto n = f.size();
      return n > 0 ? f.norm() / std::sqrt(static_cast<double>(n)) : 0.0;
    } else if (params.optimizer_options().convergence_metric == "max_atom") {
      return m_matter.maxForce();
    } else if (params.optimizer_options().convergence_metric ==
               "max_component") {
      return m_matter.getForces().cwiseAbs().maxCoeff();
    } else {
      EONC_LOG_CRITICAL("{} Unknown opt_convergence_metric: {}", "[Matter]",
                        params.optimizer_options().convergence_metric);
      throw std::invalid_argument(
          std::format("[Matter] unknown convergence_metric: {}",
                      params.optimizer_options().convergence_metric));
    }
  }
  VectorXd difference(const VectorXd &a, const VectorXd &b) {
    return m_matter.pbcV(a - b);
  }
  VectorXd getMasses() const override {
    const auto all = m_matter.getMasses();
    const auto mask = m_matter.getFree();
    VectorXd out(m_matter.numberOfFreeAtoms());
    long k = 0;
    for (long i = 0; i < m_matter.numberOfAtoms(); ++i) {
      if (mask.row(i).sum() > 0.5) {
        out[k++] = all[i];
      }
    }
    return out;
  }
  bool getPeriodic() const override { return m_matter.getPeriodic(); }
  void minimumImage(Eigen::Ref<Eigen::Vector3d> dr) const override {
    AtomMatrix m(1, 3);
    m.row(0) = dr.transpose();
    m = m_matter.pbc(m);
    dr = m.row(0).transpose();
  }
};
} // namespace

bool eonc::helpers::relaxMatter(Matter &matter, const Parameters &params,
                                bool quiet, bool writeMovie, bool checkpoint,
                                std::string prefixMovie,
                                std::string prefixCheckpoint,
                                std::vector<readcon::ConFrame> *outFrames) {
  eonc::log::Scoped m_log;
  auto objf = std::make_shared<MatterObjectiveFunction>(matter, params);
  auto optim = eonc::helpers::create::mkOptim(
      objf, params.optimizer_options().method, params);

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
    if (outFrames) {
      outFrames->push_back(eonc::io::matterToConFrame(matter, &metadata));
    }
    if (writeMovie) {
      if (!eonc::io::io_ok(matter.matter2con(min.str(), append, &metadata))) {
        QUILL_LOG_WARNING(m_log, "Failed to write movie frame {}", min.str());
      }
    }

    if (params.debug_options().write_deprecated_outs) {
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
  if (writeMovie || outFrames) {
    write_movie_frame(0, false, 0.0);
  }

  int iteration = 0;
  if (!quiet) {
    QUILL_LOG_DEBUG(m_log, "{} {:10s}  {:14s}  {:18s}  {:13s}\n", "[Matter]",
                    "Iter", "Step size",
                    params.optimizer_options().convergence_metric_label,
                    "Energy");
    QUILL_LOG_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}\n",
                    "[Matter]", iteration, 0.0, objf->getConvergence(),
                    matter.getPotentialEnergy());
  }

  while (!objf->isConverged() &&
         iteration < params.optimizer_options().max_iterations) {

    AtomMatrix pos = matter.getPositions();

    optim->step(params.optimizer_options().max_move);
    iteration++;

    double stepSize =
        eonc::geometry::maxAtomMotion(matter.pbc(matter.getPositions() - pos));

    if (!quiet) {
      QUILL_LOG_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                      "[Matter]", iteration, stepSize, objf->getConvergence(),
                      matter.getPotentialEnergy());
    }

    if (writeMovie || outFrames) {
      write_movie_frame(static_cast<uint64_t>(iteration), true, stepSize);
    }

    if (checkpoint) {
      std::ostringstream chk;
      chk << prefixCheckpoint << "_cp";
      if (!eonc::io::io_ok(matter.matter2con(chk.str(), false))) {
        QUILL_LOG_WARNING(m_log, "Failed to write checkpoint {}", chk.str());
      }
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
