#include "NEBInitialPaths.hpp"
#include "BaseStructures.h"
#include "IDPPObjectiveFunction.hpp"
#include "Optimizer.h"
#include "Parameters.h"
#include <filesystem>
#include <fstream>
#include <memory>
#include <span>
#include <vector>

#include "EonLogger.h"
namespace fs = std::filesystem;

namespace eonc::helpers::neb_paths {

// Forward declaration of ZBL setup helper to keep code clean
std::shared_ptr<Potential> createZBLPotential() {
  auto zbl_params = Parameters{};
  zbl_params.potential_options.potential = PotType::ZBL;
  // Strong short-range repulsion
  zbl_params.zbl_options.cut_inner = 0.5;
  // Cutoff sufficient to push overlapping atoms apart
  zbl_params.zbl_options.cut_global = 3.0;
  return eonc::helpers::makePotential(PotType::ZBL, zbl_params);
}

std::vector<Matter> linearPath(const Matter &initImg, const Matter &finalImg,
                               const size_t nimgs) {
  std::vector<Matter> all_images_on_path(nimgs + 2, initImg);
  all_images_on_path.front() = Matter(initImg);
  all_images_on_path.back() = Matter(finalImg);
  AtomMatrix posInitial = all_images_on_path.front().getPositions();
  AtomMatrix posFinal = all_images_on_path.back().getPositions();
  AtomMatrix imageSep = initImg.pbc(posFinal - posInitial) / (nimgs + 1);
  // Only the ones which are not the front and back
  for (auto it{std::next(all_images_on_path.begin())};
       it != std::prev(all_images_on_path.end()); ++it) {
    *it = Matter(initImg);
    (*it).setPositions(posInitial +
                       imageSep *
                           int(std::distance(all_images_on_path.begin(), it)));
  }
  return all_images_on_path;
}

std::vector<Matter> filePathInit(const std::vector<fs::path> &fsrcs,
                                 const Matter &refImg, const size_t nimgs) {
  std::vector<Matter> all_images_on_path;
  if (nimgs + 2 != fsrcs.size()) {
    throw std::runtime_error("Error in filePathInit: Expected " +
                             std::to_string(nimgs + 2) + " files, but got " +
                             std::to_string(fsrcs.size()) + ".");
  }
  all_images_on_path.reserve(nimgs + 2);
  // For all images
  for (const auto &filePath : fsrcs) {
    Matter img(refImg);
    img.con2matter(filePath.string());
    all_images_on_path.push_back(img);
  }
  return all_images_on_path;
}

std::vector<fs::path> readFilePaths(const std::string &listFilePath) {
  std::vector<fs::path> paths;
  std::ifstream inputFile(listFilePath);

  if (!inputFile.is_open()) {
    throw std::runtime_error("Error: Could not open path list file: " +
                             listFilePath);
  }

  std::string line;
  while (std::getline(inputFile, line)) {
    // Skip any empty lines in the input file
    if (!line.empty()) {
      paths.emplace_back(line);
    }
  }

  return paths;
}

MatrixXd getDistanceMatrix(const Matter &m) {
  int natoms = m.numberOfAtoms();
  MatrixXd d(natoms, natoms);
  AtomMatrix pos = m.getPositions();
  for (int i = 0; i < natoms; ++i) {
    for (int j = 0; j < natoms; ++j) {
      if (i == j) {
        d(i, j) = 0.0;
      } else {
        d(i, j) = m.pbc(pos.row(i) - pos.row(j)).norm();
      }
    }
  }
  return d;
}

std::vector<Matter> idppPath(const Matter &initImg, const Matter &finalImg,
                             const size_t nimgs, const Parameters &params,
                             bool use_zbl) {

  auto log = eonc::log::get();
  QUILL_LOG_INFO(log, "Generating initial path using IDPP...");
  if (use_zbl) {
    QUILL_LOG_WARNING(
        log, "ZBL Repulsion not implemented for iterative IDPP (idppPath). "
             "Using standard IDPP.");
  }

  // Start with a linear interpolation to get initial Cartesian coordinates
  std::vector<Matter> path = linearPath(initImg, finalImg, nimgs);

  // Pre-calculate endpoint distance matrices
  MatrixXd dInit = getDistanceMatrix(initImg);
  MatrixXd dFinal = getDistanceMatrix(finalImg);

  // Optimize intermediate images
  // Note: path[0] and path[nimgs+1] are fixed endpoints
  for (size_t i = 1; i <= nimgs; ++i) {

    // Calculate the interpolation factor (Reaction Coordinate)
    double xi = static_cast<double>(i) / (nimgs + 1);

    // Linear interpolation of the distance matrix (The "Image Dependent" part)
    MatrixXd dTarget = (1.0 - xi) * dInit + xi * dFinal;

    // Create the IDPP Objective Function
    auto idpp_objf = std::make_shared<IDPPObjectiveFunction>(
        std::make_shared<Matter>(path[i]), params, dTarget);

    // Create an Optimizer
    // Defaults to taking the same one as optimizer
    auto idpp_optim = eonc::helpers::create::mkOptim(
        idpp_objf, params.neb_options.opt_method, params);

    // Run the optimization
    int status =
        idpp_optim->run(params.neb_options.initialization.max_iterations,
                        params.neb_options.initialization.max_move);

    // Log progress
    double residual = idpp_objf->getConvergence();
    QUILL_LOG_DEBUG(log,
                    "IDPP Image {:2d}/{:2d} | xi: {:.2f} | Residual: {:.4e}", i,
                    nimgs, xi, residual);

    // Explicitly sync positions back to the path vector just to be safe
    path[i].setPositions(AtomMatrix::Map(idpp_objf->getPositions().data(),
                                         path[i].numberOfAtoms(), 3));
  }

  QUILL_LOG_INFO(log, "IDPP path generation complete.");
  return path;
}

std::vector<Matter> idppCollectivePath(const Matter &initImg,
                                       const Matter &finalImg, size_t nimgs,
                                       const Parameters &params, bool use_zbl) {
  auto log = eonc::log::get();
  QUILL_LOG_INFO(log, "Generating initial path using Collective IDPP-NEB...");

  std::vector<Matter> path = linearPath(initImg, finalImg, nimgs);

  // 1. Base Objective
  std::shared_ptr<ObjectiveFunction> idpp_objf =
      std::make_shared<CollectiveIDPPObjectiveFunction>(path, params);

  // 2. Optional ZBL Wrapper
  if (use_zbl) {
    QUILL_LOG_INFO(log, "Enabling ZBL repulsive penalty for IDPP...");
    auto zbl_pot = createZBLPotential();
    // Wrap the IDPP objective with ZBL repulsion (weight = 1.0)
    idpp_objf = std::make_shared<ZBLRepulsiveIDPPObjective>(idpp_objf, zbl_pot,
                                                            path, params, 1.0);
  }

  auto optim = eonc::helpers::create::mkOptim(
      idpp_objf, params.neb_options.initialization.opt_method, params);

  int maxSteps = params.neb_options.initialization.max_iterations;
  int currentStep = 0;
  int checkInterval = 40;

  while (currentStep < maxSteps) {
    optim->run(checkInterval, params.optimizer_options.max_move);
    currentStep += checkInterval;

    if (idpp_objf->isConverged()) {
      QUILL_LOG_INFO(log,
                     "IDPP-NEB converged after {} steps. Max Residual: {:.4f}",
                     currentStep, idpp_objf->getConvergence());
      return path;
    }
  }

  QUILL_LOG_WARNING(log,
                    "IDPP-NEB reached max_iterations ({}) without full "
                    "convergence. Residual: {:.4f}",
                    maxSteps, idpp_objf->getConvergence());
  return path;
}

// Helper to insert an image linearly between two others
Matter interpolateImage(const Matter &A, const Matter &B, double fraction) {
  Matter newImg(A);
  AtomMatrix posA = A.getPositions();
  AtomMatrix posB = B.getPositions();
  // Use PBC-aware interpolation
  AtomMatrix diff = A.pbc(posB - posA);
  newImg.setPositions(posA + fraction * diff);
  return newImg;
}

std::vector<Matter> sidppPath(const Matter &initImg, const Matter &finalImg,
                              size_t target_nimgs, const Parameters &params,
                              bool use_zbl) {

  auto log = eonc::log::get();
  const auto &init = params.neb_options.initialization;
  QUILL_LOG_INFO(log, "Generating initial path using S-IDPP{} ({} images, "
                      "alpha={:.2f}, frontier_tol={:.4f})...",
                 use_zbl ? "-ZBL" : "", target_nimgs, init.sidpp_alpha,
                 init.sidpp_frontier_tol);

  // 1. Start with endpoints [Reactant, Product]
  std::vector<Matter> path;
  path.push_back(initImg);
  path.push_back(finalImg);

  std::shared_ptr<Potential> zbl_pot = nullptr;
  if (use_zbl) {
    zbl_pot = createZBLPotential();
  }

  // Track frontier counts: nLeft images from reactant, nRight from product
  int nLeft = 0;
  int nRight = 0;
  int nIntermediate = 0;
  bool addToLeft = true; // Alternate L/R, starting with left

  // Helper: create IDPP objective with optional ZBL wrapping
  auto makeIDPP = [&]() -> std::shared_ptr<ObjectiveFunction> {
    auto objf =
        std::make_shared<CollectiveIDPPObjectiveFunction>(path, params);
    if (use_zbl && zbl_pot) {
      return std::make_shared<ZBLRepulsiveIDPPObjective>(objf, zbl_pot, path,
                                                         params, 1.0);
    }
    return objf;
  };

  // Helper: relax current path on IDPP surface
  auto relaxPath = [&](int maxSteps) -> double {
    auto objf = makeIDPP();
    auto optim =
        eonc::helpers::create::mkOptim(objf, init.opt_method, params);
    int step = 0;
    while (step < maxSteps) {
      optim->run(5, init.max_move);
      step += 5;
      if (objf->isConverged())
        break;
    }
    return objf->getConvergence();
  };

  // 2. Sequential growth loop: alternate adding images from L and R
  while (nIntermediate < static_cast<int>(target_nimgs)) {

    if (addToLeft && nIntermediate < static_cast<int>(target_nimgs)) {
      // Add to left (reactant) frontier
      Matter &frontier = path[nLeft];
      Matter &next = path[nLeft + 1];
      Matter newImg = interpolateImage(frontier, next, init.sidpp_alpha);
      path.insert(path.begin() + nLeft + 1, newImg);
      nLeft++;
      nIntermediate++;
      QUILL_LOG_DEBUG(log, "S-IDPP: +L frontier (nL={}, nR={}, total={})",
                      nLeft, nRight, nIntermediate);
    } else if (nIntermediate < static_cast<int>(target_nimgs)) {
      // Add to right (product) frontier
      int rightIdx = static_cast<int>(path.size()) - 1 - nRight;
      Matter &frontier = path[rightIdx];
      Matter &prev = path[rightIdx - 1];
      Matter newImg = interpolateImage(frontier, prev, init.sidpp_alpha);
      path.insert(path.begin() + rightIdx, newImg);
      nRight++;
      nIntermediate++;
      QUILL_LOG_DEBUG(log, "S-IDPP: +R frontier (nL={}, nR={}, total={})",
                      nLeft, nRight, nIntermediate);
    }
    addToLeft = !addToLeft; // Alternate sides

    // Optimize current path on IDPP surface
    double residual = relaxPath(init.nsteps);

    // Frontier convergence gating: if not converged, keep relaxing
    // before adding more images (up to max_iterations total)
    if (residual > init.sidpp_frontier_tol) {
      double residual2 = relaxPath(init.max_iterations - init.nsteps);
      QUILL_LOG_DEBUG(log,
                      "S-IDPP: Extended relaxation {:.4f} -> {:.4f}",
                      residual, residual2);
      residual = residual2;
    }

    QUILL_LOG_DEBUG(log,
                    "S-IDPP: {} images | Residual: {:.4f}",
                    nIntermediate, residual);
  }

  // 3. Reparameterize: redistribute images evenly along arc length
  if (init.sidpp_reparam && path.size() > 3) {
    QUILL_LOG_INFO(log, "S-IDPP: Reparameterizing {} images along arc length",
                   path.size() - 2);
    path = resamplePath(path, target_nimgs);
  }

  // 4. Final relaxation of the complete path
  QUILL_LOG_INFO(log, "S-IDPP: Final relaxation of full path...");
  double finalResidual = relaxPath(init.max_iterations);
  QUILL_LOG_INFO(log, "S-IDPP: Final residual: {:.4f}", finalResidual);

  return path;
}

AtomMatrix cubicInterpolate(const AtomMatrix &P0, const AtomMatrix &T0,
                            const AtomMatrix &P1, const AtomMatrix &T1,
                            double f) {
  double f2 = f * f;
  double f3 = f2 * f;

  // Hermite basis functions
  double h00 = 2 * f3 - 3 * f2 + 1;
  double h10 = f3 - 2 * f2 + f;
  double h01 = -2 * f3 + 3 * f2;
  double h11 = f3 - f2;

  return h00 * P0 + h10 * T0 + h01 * P1 + h11 * T1;
}

std::vector<Matter> resamplePath(const std::vector<Matter> &densePath,
                                 size_t targetCount) {
  if (densePath.size() == targetCount + 2)
    return densePath;

  size_t n = densePath.size();

  // Calculate cumulative arc length along the path
  std::vector<double> arcLength(n, 0.0);
  for (size_t i = 1; i < n; ++i) {
    AtomMatrix diff = densePath[i].pbc(densePath[i].getPositions() -
                                       densePath[i - 1].getPositions());
    arcLength[i] = arcLength[i - 1] + diff.norm();
  }
  double totalLength = arcLength.back();

  // Calculate tangents for cubic interpolation
  std::vector<AtomMatrix> tangents(n);
  for (size_t i = 0; i < n; ++i) {
    AtomMatrix T;
    if (i == 0) {
      T = densePath[i].pbc(densePath[i + 1].getPositions() -
                           densePath[i].getPositions());
    } else if (i == n - 1) {
      T = densePath[i].pbc(densePath[i].getPositions() -
                           densePath[i - 1].getPositions());
    } else {
      AtomMatrix dNext = densePath[i].pbc(densePath[i + 1].getPositions() -
                                          densePath[i].getPositions());
      AtomMatrix dPrev = densePath[i].pbc(densePath[i].getPositions() -
                                          densePath[i - 1].getPositions());
      T = 0.5 * (dNext + dPrev);
    }
    // Scale tangent by local segment length for proper spline parameterization
    if (i > 0 && i < n - 1) {
      double localScale = (arcLength[i + 1] - arcLength[i - 1]) / 2.0;
      T = T.normalized() * localScale;
    }
    tangents[i] = T;
  }

  std::vector<Matter> resampled;
  resampled.reserve(targetCount + 2);
  resampled.push_back(densePath.front());

  // Place new images at equal arc-length intervals
  double segmentLength = totalLength / (targetCount + 1);

  for (size_t i = 1; i <= targetCount; ++i) {
    double targetArc = i * segmentLength;

    // Find the segment containing this arc length
    size_t lowIdx = 0;
    for (size_t j = 1; j < n; ++j) {
      if (arcLength[j] >= targetArc) {
        lowIdx = j - 1;
        break;
      }
    }
    size_t highIdx = lowIdx + 1;

    // Interpolation parameter within this segment
    double segmentArc = arcLength[highIdx] - arcLength[lowIdx];
    double f = (segmentArc > 1e-10)
                   ? (targetArc - arcLength[lowIdx]) / segmentArc
                   : 0.0;

    Matter newImg(densePath[0]);
    newImg.setPositions(cubicInterpolate(
        densePath[lowIdx].getPositions(), tangents[lowIdx],
        densePath[highIdx].getPositions(), tangents[highIdx], f));
    resampled.push_back(newImg);
  }

  resampled.push_back(densePath.back());
  return resampled;
}

void resamplePathInPlace(std::span<std::shared_ptr<Matter>> path) {
  size_t n = path.size();
  if (n < 3)
    return;

  // Calculate cumulative arc length
  std::vector<double> arcLength(n, 0.0);
  for (size_t i = 1; i < n; ++i) {
    AtomMatrix diff =
        path[i]->pbc(path[i]->getPositions() - path[i - 1]->getPositions());
    arcLength[i] = arcLength[i - 1] + diff.norm();
  }
  double totalLength = arcLength.back();
  if (totalLength < 1e-12)
    return;

  // Tangents for cubic interpolation
  std::vector<AtomMatrix> tangents(n);
  for (size_t i = 0; i < n; ++i) {
    AtomMatrix T;
    if (i == 0) {
      T = path[i]->pbc(path[i + 1]->getPositions() -
                        path[i]->getPositions());
    } else if (i == n - 1) {
      T = path[i]->pbc(path[i]->getPositions() -
                        path[i - 1]->getPositions());
    } else {
      AtomMatrix dN = path[i]->pbc(path[i + 1]->getPositions() -
                                    path[i]->getPositions());
      AtomMatrix dP = path[i]->pbc(path[i]->getPositions() -
                                    path[i - 1]->getPositions());
      T = 0.5 * (dN + dP);
    }
    if (i > 0 && i < n - 1) {
      double localScale = (arcLength[i + 1] - arcLength[i - 1]) / 2.0;
      T = T.normalized() * localScale;
    }
    tangents[i] = T;
  }

  // Store original positions for interpolation source
  std::vector<AtomMatrix> origPos(n);
  for (size_t i = 0; i < n; ++i)
    origPos[i] = path[i]->getPositions();

  // Redistribute interior images at equal arc-length intervals (in-place)
  size_t nInterior = n - 2;
  double segLen = totalLength / (nInterior + 1);

  for (size_t i = 1; i <= nInterior; ++i) {
    double targetArc = i * segLen;
    size_t lo = 0;
    for (size_t j = 1; j < n; ++j) {
      if (arcLength[j] >= targetArc) {
        lo = j - 1;
        break;
      }
    }
    size_t hi = lo + 1;
    double sArc = arcLength[hi] - arcLength[lo];
    double f = (sArc > 1e-10) ? (targetArc - arcLength[lo]) / sArc : 0.0;

    path[i]->setPositions(
        cubicInterpolate(origPos[lo], tangents[lo], origPos[hi], tangents[hi], f));
  }
}

} // namespace eonc::helpers::neb_paths
