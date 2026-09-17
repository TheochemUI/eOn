#pragma once
#include "Matter.h"
#include <filesystem>
#include <format>
#include <memory>
#include <span>
#include <stdexcept>
#include <string_view>

namespace eonc {

namespace helpers::neb_paths {
namespace fs = std::filesystem;

/// Abort before Eigen subtracts two position matrices of different size.
inline void requireSameAtomCount(const Matter &a, const Matter &b,
                                 std::string_view what) {
  if (a.numberOfAtoms() == b.numberOfAtoms()) {
    return;
  }
  throw std::invalid_argument(std::format(
      "NEB: {} do not have the same number of atoms ({} vs {})", what,
      a.numberOfAtoms(), b.numberOfAtoms()));
}

std::vector<Matter> linearPath(const Matter &initImg, const Matter &finalImg,
                               const size_t nimgs);

std::vector<Matter> filePathInit(const std::vector<fs::path> &fsrcs,
                                 const Matter &refImg, const size_t nimgs);

/**
 * @brief Interpolates positions using a cubic Hermite spline.
 * @param P0 Starting positions
 * @param T0 Tangent at P0 (scaled by segment length)
 * @param P1 Ending positions
 * @param T1 Tangent at P1 (scaled by segment length)
 * @param f Fraction between 0 and 1
 */
AtomMatrix cubicInterpolate(const AtomMatrix &P0, const AtomMatrix &T0,
                            const AtomMatrix &P1, const AtomMatrix &T1,
                            double f);

std::vector<Matter> resamplePath(const std::vector<Matter> &densePath,
                                 size_t targetCount);

/// In-place path reparameterization for NEB shared_ptr paths.
/// Redistributes interior images at equal arc-length intervals using
/// cubic Hermite interpolation, without allocating new Matter objects.
///
/// @param path The FULL path including endpoints. path[0] and path[n-1]
///             are treated as fixed endpoints and are never modified.
///             Interior images path[1] through path[n-2] are repositioned.
///             Pass the entire NEB path vector, not a sub-span of
///             intermediate images.
void resamplePathInPlace(std::span<std::shared_ptr<Matter>> path);

/**
 * @brief Reads a file where each line contains a path to another file.
 *
 * @param listFilePath The path to the file containing the list of file paths.
 * @return A vector of filesystem paths. Returns an empty vector if the
 * file cannot be opened.
 */
std::vector<std::filesystem::path>
readFilePaths(const std::string &listFilePath);

MatrixXd getDistanceMatrix(const Matter &m);

std::vector<Matter> idppPath(const Matter &initImg, const Matter &finalImg,
                             size_t nimgs, const Parameters &params,
                             bool use_zbl = false);

std::vector<Matter> idppCollectivePath(const Matter &initImg,
                                       const Matter &finalImg, size_t nimgs,
                                       const Parameters &params,
                                       bool use_zbl = false);

std::vector<Matter> sidppPath(const Matter &initImg, const Matter &finalImg,
                              size_t target_nimgs, const Parameters &params,
                              bool use_zbl = false);

/// Adjacent images closer than min_sep (RMSD, PBC) are a collapsed path.
/// SIDPP + resample can emit bit-identical intermediates; those starve
/// climbing_image_converged_only (eOn-bghy).
void ensureDistinctAdjacentImages(const std::vector<Matter> &path,
                                  double min_sep);

// Helper to insert an image linearly between two others
Matter interpolateImage(const Matter &A, const Matter &B, double fraction);
// Helper to construct ZBL potentials
std::shared_ptr<Potential> createZBLPotential();

} // namespace helpers::neb_paths

} // namespace eonc
