#include "NEBInitialPaths.hpp"
#include "BaseStructures.h"
#include "IDPPObjectiveFunction.hpp"
#include "Optimizer.h"
#include "Parameters.h"
#include <filesystem>
#include <fstream>
#include <vector>

namespace fs = std::filesystem;

namespace helper_functions::neb_paths {

// Forward declaration of ZBL setup helper to keep code clean
std::shared_ptr<Potential> createZBLPotential() {
  auto zbl_params = std::make_shared<Parameters>();
  zbl_params->potential = PotType::ZBL;
  // Strong short-range repulsion
  zbl_params->zbl_options.cut_inner = 0.5;
  // Cutoff sufficient to push overlapping atoms apart
  zbl_params->zbl_options.cut_global = 3.0;
  return helper_functions::makePotential(PotType::ZBL, zbl_params);
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

Eigen::MatrixXd getDistanceMatrix(const Matter &m) {
  int natoms = m.numberOfAtoms();
  Eigen::MatrixXd d(natoms, natoms);
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
                             const size_t nimgs,
                             std::shared_ptr<Parameters> params, bool use_zbl) {

  auto log = spdlog::get("combi");
  SPDLOG_LOGGER_INFO(log, "Generating initial path using IDPP...");
  if (use_zbl) {
    SPDLOG_LOGGER_WARN(
        log, "ZBL Repulsion not implemented for iterative IDPP (idppPath). "
             "Using standard IDPP.");
  }

  // Start with a linear interpolation to get initial Cartesian coordinates
  std::vector<Matter> path = linearPath(initImg, finalImg, nimgs);

  // Pre-calculate endpoint distance matrices
  Eigen::MatrixXd dInit = getDistanceMatrix(initImg);
  Eigen::MatrixXd dFinal = getDistanceMatrix(finalImg);

  // Optimize intermediate images
  // Note: path[0] and path[nimgs+1] are fixed endpoints
  for (size_t i = 1; i <= nimgs; ++i) {

    // Calculate the interpolation factor (Reaction Coordinate)
    double xi = static_cast<double>(i) / (nimgs + 1);

    // Linear interpolation of the distance matrix (The "Image Dependent" part)
    Eigen::MatrixXd dTarget = (1.0 - xi) * dInit + xi * dFinal;

    // Create the IDPP Objective Function
    auto idpp_objf = std::make_shared<IDPPObjectiveFunction>(
        std::make_shared<Matter>(path[i]), params, dTarget);

    // Create an Optimizer
    // Defaults to taking the same one as optimizer
    auto idpp_optim = helpers::create::mkOptim(
        idpp_objf, params->neb_options.opt_method, params);

    // Run the optimization
    int status =
        idpp_optim->run(params->neb_options.initialization.max_iterations,
                        params->neb_options.initialization.max_move);

    // Log progress
    double residual = idpp_objf->getConvergence();
    SPDLOG_LOGGER_DEBUG(
        log, "IDPP Image {:2d}/{:2d} | xi: {:.2f} | Residual: {:.4e}", i, nimgs,
        xi, residual);

    // Explicitly sync positions back to the path vector just to be safe
    path[i].setPositions(MatrixXd::Map(idpp_objf->getPositions().data(),
                                       path[i].numberOfAtoms(), 3));
  }

  SPDLOG_LOGGER_INFO(log, "IDPP path generation complete.");
  return path;
}

std::vector<Matter> idppCollectivePath(const Matter &initImg,
                                       const Matter &finalImg, size_t nimgs,
                                       std::shared_ptr<Parameters> params,
                                       bool use_zbl) {
  auto log = spdlog::get("combi");
  SPDLOG_LOGGER_INFO(log,
                     "Generating initial path using Collective IDPP-NEB...");

  std::vector<Matter> path = linearPath(initImg, finalImg, nimgs);

  // 1. Base Objective
  std::shared_ptr<ObjectiveFunction> idpp_objf =
      std::make_shared<CollectiveIDPPObjectiveFunction>(path, params);

  // 2. Optional ZBL Wrapper
  if (use_zbl) {
    SPDLOG_LOGGER_INFO(log, "Enabling ZBL repulsive penalty for IDPP...");
    auto zbl_pot = createZBLPotential();
    // Wrap the IDPP objective with ZBL repulsion (weight = 1.0)
    idpp_objf = std::make_shared<ZBLRepulsiveIDPPObjective>(idpp_objf, zbl_pot,
                                                            path, params, 1.0);
  }

  auto optim = helpers::create::mkOptim(
      idpp_objf, params->neb_options.initialization.opt_method, params);

  int maxSteps = params->neb_options.initialization.max_iterations;
  int currentStep = 0;
  int checkInterval = 40;

  while (currentStep < maxSteps) {
    optim->run(checkInterval, params->optMaxMove);
    currentStep += checkInterval;

    if (idpp_objf->isConverged()) {
      SPDLOG_LOGGER_INFO(
          log, "IDPP-NEB converged after {} steps. Max Residual: {:.4f}",
          currentStep, idpp_objf->getConvergence());
      return path;
    }
  }

  SPDLOG_LOGGER_WARN(log,
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
                              size_t target_nimgs,
                              std::shared_ptr<Parameters> params,
                              bool use_zbl) {

  auto log = spdlog::get("combi");
  if (use_zbl) {
    SPDLOG_LOGGER_INFO(
        log, "Generating initial path using Sequential IDPP-ZBL (S-IDPP)...");
  } else {
    SPDLOG_LOGGER_INFO(
        log, "Generating initial path using Sequential IDPP (S-IDPP)...");
  }

  // 1. Start with [Init, Final]
  std::vector<Matter> path;
  path.push_back(initImg);
  path.push_back(finalImg);

  // Initialize ZBL potential if requested (created once to be efficient)
  std::shared_ptr<Potential> zbl_pot = nullptr;
  if (use_zbl) {
    zbl_pot = createZBLPotential();
  }

  // Track how many images we have added to the Reactant (left) and Product
  // (right) sides In ORCA this is nR and nP. We start with 0 intermediate
  // images.
  int nLeft = 0;
  int nRight = 0;
  int nIntermediate = 0;

  // 2. Growth Loop
  while (nIntermediate < target_nimgs) {

    // --- STEP A: ADD IMAGES ---
    // We add images if we haven't reached the target yet.
    // Strategy: Add one to left, one to right (if space permits)

    // Add to Left (Reaction side)
    if (nIntermediate < target_nimgs) {
      // Insert after the last "Left" image (index nLeft)
      // Interpolate between path[nLeft] and path[nLeft+1] (which is the first
      // "Right" image) For the very first step, this interpolates between Init
      // and Final.

      // We place it close to the frontier (e.g., 20% of the way to the next
      // image) effectively "growing" slowly.
      Matter frontier = path[nLeft];
      Matter next = path[nLeft + 1];
      Matter newImg = interpolateImage(
          frontier, next, params->neb_options.initialization.sidpp_alpha);

      path.insert(path.begin() + nLeft + 1, newImg);
      nLeft++;
      nIntermediate++;
      SPDLOG_LOGGER_DEBUG(log, "S-IDPP: Added Left Frontier. Total: {}",
                          nIntermediate);
    }

    // Add to Right (Product side)
    if (nIntermediate < target_nimgs) {
      // Insert before the first "Right" image.
      // Current indices: [0 ... nLeft] [NEW] [nLeft+1 ... END]
      // We want to insert before the last element (Final).

      int rightFrontierIdx = path.size() - 1 - nRight;
      Matter frontier = path[rightFrontierIdx];
      Matter prev = path[rightFrontierIdx - 1];

      // Grow backwards from product
      Matter newImg = interpolateImage(
          frontier, prev, params->neb_options.initialization.sidpp_alpha);

      path.insert(path.begin() + rightFrontierIdx, newImg);
      nRight++;
      nIntermediate++;
      SPDLOG_LOGGER_DEBUG(log, "S-IDPP: Added Right Frontier. Total: {}",
                          nIntermediate);
    }

    // --- STEP B: OPTIMIZE CURRENT SET ---
    int steps = params->neb_options.initialization.nsteps;

    // Create Base IDPP Objective
    std::shared_ptr<ObjectiveFunction> idpp_objf =
        std::make_shared<CollectiveIDPPObjectiveFunction>(path, params);

    // Apply ZBL Wrapper if enabled
    if (use_zbl && zbl_pot) {
      idpp_objf = std::make_shared<ZBLRepulsiveIDPPObjective>(
          idpp_objf, zbl_pot, path, params, 1.0);
    }

    // TODO(rg): this is a headache, since it uses the optimizer stanza but with
    // the NEB OptType
    auto optim = helpers::create::mkOptim(
        idpp_objf, params->neb_options.initialization.opt_method, params);

    // Relax the current intermediate path until it meets the tolerance
    int growthRelaxSteps = params->neb_options.initialization.max_iterations;
    int step = 0;
    while (step < growthRelaxSteps) {
      optim->run(5, params->optMaxMove); // Run small batches
      step += 5;
      if (idpp_objf->isConverged())
        break;
    }

    SPDLOG_LOGGER_DEBUG(
        log,
        "S-IDPP Frontier Relaxed: {} images | Steps: {} | Residual: {:.4f}",
        nIntermediate, step, idpp_objf->getConvergence());

    SPDLOG_LOGGER_DEBUG(log, "S-IDPP: Relaxed with {} images. Residual: {:.4f}",
                        nIntermediate, idpp_objf->getConvergence());
  }

  // 3. Final Interpolation / Alignment
  // The path now has size target_nimgs + 2.
  // However, the "growth" heuristic might have left the middle gap uneven.
  // It is good practice to run one final IDPP on the FULL path to evenly space
  // everything.

  SPDLOG_LOGGER_INFO(log, "S-IDPP: Final relaxation of full path...");

  std::shared_ptr<ObjectiveFunction> final_objf =
      std::make_shared<CollectiveIDPPObjectiveFunction>(path, params);

  // Apply ZBL Wrapper for final relaxation as well
  if (use_zbl && zbl_pot) {
    final_objf = std::make_shared<ZBLRepulsiveIDPPObjective>(
        final_objf, zbl_pot, path, params, 1.0);
  }

  auto final_optim =
      helpers::create::mkOptim(final_objf, OptType::LBFGS, params);
  final_optim->run(500, params->optMaxMove);

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

} // namespace helper_functions::neb_paths
