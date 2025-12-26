#pragma once
#include "Matter.h"
#include <filesystem>

namespace helper_functions::neb_paths {
namespace fs = std::filesystem;
std::vector<Matter> linearPath(const Matter &initImg, const Matter &finalImg,
                               const size_t nimgs);

std::vector<Matter> filePathInit(const std::vector<fs::path> &fsrcs,
                                 const Matter &refImg, const size_t nimgs);

/**
 * @brief Reads a file where each line contains a path to another file.
 *
 * @param listFilePath The path to the file containing the list of file paths.
 * @return A vector of filesystem paths. Returns an empty vector if the
 * file cannot be opened.
 */
std::vector<std::filesystem::path>
readFilePaths(const std::string &listFilePath);

Eigen::MatrixXd getDistanceMatrix(const Matter &m);

std::vector<Matter> idppPath(const Matter &initImg, const Matter &finalImg,
                             size_t nimgs, std::shared_ptr<Parameters> params,
                             bool use_zbl = false);

std::vector<Matter> idppCollectivePath(const Matter &initImg,
                                       const Matter &finalImg, size_t nimgs,
                                       std::shared_ptr<Parameters> params,
                                       bool use_zbl = false);

std::vector<Matter> sidppPath(const Matter &initImg, const Matter &finalImg,
                              size_t target_nimgs,
                              std::shared_ptr<Parameters> params,
                              bool use_zbl = false);

// Helper to insert an image linearly between two others
Matter interpolateImage(const Matter &A, const Matter &B, double fraction);

} // namespace helper_functions::neb_paths
