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
#include "NEBSplineExtrema.h"
#include "ConFileIO.h"
#include <cmath>
#include <format>
#include <fstream>
#include <limits>

namespace fs = std::filesystem;

namespace eonc::neb {

namespace {

eonc::io::ConFrameMetadata neb_frame_metadata(
    const std::vector<std::shared_ptr<Matter>> &path,
    const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
    const std::vector<std::shared_ptr<EigenmodeStrategy>> &eigenmode_solvers,
    long numImages, bool estimateEigenvalues, long imageIndex,
    double reactionCoordinate, std::optional<size_t> bandIndex) {
  AtomMatrix tang;
  if (imageIndex == 0) {
    tang = path[0]->pbc(path[1]->getPositions() - path[0]->getPositions());
  } else if (imageIndex == numImages + 1) {
    tang = path[numImages]->pbc(path[numImages + 1]->getPositions() -
                                path[numImages]->getPositions());
  } else {
    tang = *tangent[imageIndex];
  }
  tang.normalize();

  const double reference_energy = path[0]->getPotentialEnergy();
  const double absolute_energy = path[imageIndex]->getPotentialEnergy();
  const double relative_energy = absolute_energy - reference_energy;
  const double parallel_force = matDot(path[imageIndex]->getForces(), tang);

  eonc::io::ConFrameMetadata metadata;
  metadata.frame_index = static_cast<uint64_t>(imageIndex);
  metadata.energy = absolute_energy;
  metadata.neb_bead = static_cast<uint64_t>(imageIndex);
  if (bandIndex) {
    metadata.neb_band = static_cast<uint64_t>(*bandIndex);
  }
  metadata.scalars.push_back({"reaction_coordinate", reactionCoordinate});
  metadata.scalars.push_back({"relative_energy", relative_energy});
  metadata.scalars.push_back({"parallel_force", parallel_force});

  if (estimateEigenvalues && imageIndex >= 0 &&
      imageIndex < static_cast<long>(eigenmode_solvers.size()) &&
      eigenmode_solvers[imageIndex]) {
    eonc::eigenmodeCompute(*eigenmode_solvers[imageIndex], path[imageIndex],
                           tang);
    metadata.scalars.push_back(
        {"lowest_eigenvalue",
         eonc::eigenmodeGetEigenvalue(*eigenmode_solvers[imageIndex])});
  }

  return metadata;
}

} // namespace

ExtremaResult
findSplineExtrema(const std::vector<std::shared_ptr<Matter>> &path,
                  const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
                  long numImages) {

  auto *log = eonc::log::get();

  // Calculate cubic parameters for each interval
  AtomMatrix tangentEndpoint;
  std::vector<double> a(numImages + 1), b(numImages + 1), c(numImages + 1),
      d(numImages + 1);
  double F1, F2, U1, U2, dist;

  for (long i = 0; i <= numImages; i++) {
    dist = path[i]->distanceTo(*path[i + 1]);
    if (i == 0) {
      tangentEndpoint =
          path[i]->pbc(path[1]->getPositions() - path[0]->getPositions());
      tangentEndpoint.normalize();
      F1 = matDot(path[i]->getForces(), tangentEndpoint) * dist;
    } else {
      F1 = matDot(path[i]->getForces(), *tangent[i]) * dist;
    }
    if (i == numImages) {
      tangentEndpoint = path[i + 1]->pbc(path[numImages + 1]->getPositions() -
                                         path[numImages]->getPositions());
      tangentEndpoint.normalize();
      F2 = matDot(path[i + 1]->getForces(), tangentEndpoint) * dist;
    } else {
      F2 = matDot(path[i + 1]->getForces(), *tangent[i + 1]) * dist;
    }
    U1 = path[i]->getPotentialEnergy();
    U2 = path[i + 1]->getPotentialEnergy();
    a[i] = U1;
    b[i] = -F1;
    c[i] = 3. * (U2 - U1) + 2. * F1 + F2;
    d[i] = -2. * (U2 - U1) - (F1 + F2);
  }

  ExtremaResult result;
  result.positions.resize(2 * (numImages + 1));
  result.energies.resize(2 * (numImages + 1));
  result.curvatures.resize(2 * (numImages + 1));

  double discriminant, f;

  for (long i = 0; i <= numImages; i++) {
    discriminant = c[i] * c[i] - 3.0 * b[i] * d[i];
    if (discriminant >= 0) {
      f = -1;

      // Quadratic case
      if ((d[i] == 0) && (c[i] != 0)) {
        f = (-b[i] / (2. * c[i]));
      }
      // Cubic case 1
      else if (d[i] != 0) {
        f = -(c[i] + std::sqrt(discriminant)) / (3. * d[i]);
      }
      if ((f >= 0) && (f <= 1)) {
        result.positions[result.numExtrema] = i + f;
        result.energies[result.numExtrema] =
            ((d[i] * f + c[i]) * f + b[i]) * f + a[i]; // Horner's method
        result.curvatures[result.numExtrema] = 6.0 * d[i] * f + 2 * c[i];
        result.numExtrema++;
      }
      // Cubic case 2
      if (d[i] != 0) {
        f = (-(c[i] - std::sqrt(discriminant)) / (3. * d[i]));
      }
      if ((f >= 0) && (f <= 1)) {
        result.positions[result.numExtrema] = i + f;
        result.energies[result.numExtrema] =
            ((d[i] * f + c[i]) * f + b[i]) * f + a[i]; // Horner's method
        result.curvatures[result.numExtrema] = 6 * d[i] * f + 2 * c[i];
        result.numExtrema++;
      }
    }
  }

  QUILL_LOG_DEBUG(log, "Found {} extrema", result.numExtrema);
  QUILL_LOG_DEBUG(log, "Energy reference: {}", path[0]->getPotentialEnergy());
  for (long i = 0; i < result.numExtrema; i++) {
    QUILL_LOG_DEBUG(
        log, "extrema #{} at image position {} with energy {} and curvature {}",
        i + 1, result.positions[i],
        result.energies[i] - path[0]->getPotentialEnergy(),
        result.curvatures[i]);
  }

  return result;
}

void printImageData(
    const std::vector<std::shared_ptr<Matter>> &path,
    const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
    const std::vector<std::shared_ptr<EigenmodeStrategy>> &eigenmode_solvers,
    long numImages, bool estimateEigenvalues, bool writeToFile, size_t idx,
    eonc::log::Scoped log) {

  double dist, distTotal = 0;
  AtomMatrix tangentStart =
      path[0]->pbc(path[1]->getPositions() - path[0]->getPositions());
  AtomMatrix tangentEnd = path[numImages]->pbc(
      path[numImages + 1]->getPositions() - path[numImages]->getPositions());
  AtomMatrix tang;
  std::string header;
  if (estimateEigenvalues) {
    header = std::format("{:>3s} {:>12s} {:>12s} {:>12s} {:>12s}", "img",
                         "rxn_coord", "energy", "f_para", "eigval");
  } else {
    header = std::format("{:>3s} {:>12s} {:>12s} {:>12s}", "img", "rxn_coord",
                         "energy", "f_para");
  }

  tangentStart.normalize();
  tangentEnd.normalize();

  std::ofstream fileLogger;
  if (writeToFile) {
    std::string neb_dat_fs;
    if (idx == std::numeric_limits<size_t>::max()) {
      neb_dat_fs = "neb.dat";
    } else {
      neb_dat_fs = std::format("neb_{:03}.dat", idx);
    }
    if (fs::exists(neb_dat_fs)) {
      fs::remove(neb_dat_fs);
    }
    fileLogger.open(neb_dat_fs);
    if (fileLogger.is_open()) {
      fileLogger << header << "\n";
    }
  }
  const double energy_reactant = path[0]->getPotentialEnergy();

  for (long i = 0; i <= numImages + 1; i++) {
    if (i == 0) {
      tang = tangentStart;
    } else if (i == numImages + 1) {
      tang = tangentEnd;
    } else {
      tang = *tangent[i];
    }

    if (i > 0) {
      dist = path[i]->distanceTo(*path[i - 1]);
      distTotal += dist;
    }

    double relative_energy = path[i]->getPotentialEnergy() - energy_reactant;
    double parallel_force = matDot(path[i]->getForces(), tang);

    if (estimateEigenvalues) {
      eonc::eigenmodeCompute(*eigenmode_solvers[i], path[i], tang);
      double lowest_eigenvalue =
          eonc::eigenmodeGetEigenvalue(*eigenmode_solvers[i]);
      if (fileLogger.is_open()) {
        fileLogger << std::format(
            "{:>3} {:>12.6f} {:>12.6f} {:>12.6f} {:>12.6f}\n", i, distTotal,
            relative_energy, parallel_force, lowest_eigenvalue);
      } else {
        QUILL_LOG_DEBUG(log, "{:>3} {:>12.6f} {:>12.6f} {:>12.6f} {:>12.6f}", i,
                        distTotal, relative_energy, parallel_force,
                        lowest_eigenvalue);
      }
    } else {
      if (fileLogger.is_open()) {
        fileLogger << std::format("{:>3} {:>12.6f} {:>12.6f} {:>12.6f}\n", i,
                                  distTotal, relative_energy, parallel_force);
      } else {
        QUILL_LOG_DEBUG(log, "{:>3} {:>12.6f} {:>12.6f} {:>12.6f}", i,
                        distTotal, relative_energy, parallel_force);
      }
    }
  }
}

bool writePathCon(
    const std::vector<std::shared_ptr<Matter>> &path,
    const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
    const std::vector<std::shared_ptr<EigenmodeStrategy>> &eigenmode_solvers,
    long numImages, bool estimateEigenvalues, std::string filename,
    std::optional<size_t> bandIndex) {
  double distTotal = 0.0;

  for (long i = 0; i <= numImages + 1; i++) {
    if (i > 0) {
      distTotal += path[i]->distanceTo(*path[i - 1]);
    }

    auto metadata = neb_frame_metadata(path, tangent, eigenmode_solvers,
                                       numImages, estimateEigenvalues, i,
                                       distTotal, bandIndex);
    if (!path[i]->matter2con(filename, /*append=*/i > 0, &metadata)) {
      return false;
    }
  }

  return true;
}

} // namespace eonc::neb
