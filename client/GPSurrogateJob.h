#pragma once

#include "HelperFunctions.h"
#include "Job.h"
#include "Parameters.h"

#ifdef WITH_CATLEARN
#include "potentials/CatLearnPot/CatLearnPot.h"
#endif
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>

#include "NudgedElasticBand.h"

class GPSurrogateJob : public Job {
public:
  GPSurrogateJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {
    // debugging
#ifndef NDEBUG
    py::module_ sys_mod = py::module_::import("sys");
    py::module_ ipdb_mod = py::module_::import("ipdb");
    sys_mod.attr("breakpointhook") = ipdb_mod.attr("set_trace");
#endif // NDEBUG
  }
  ~GPSurrogateJob(void) = default;
  std::vector<std::string> run(void) override;

private:
  void saveData(NudgedElasticBand::NEBStatus status,
                std::unique_ptr<NudgedElasticBand> neb);
  std::vector<std::string> returnFiles;
#ifndef WITH_ASE_ORCA
  // When ASE_ORCA is used as a potential ClientEON.cpp holds lock in main
  pybind11::scoped_interpreter guard{};
#endif
};

namespace helper_functions::surrogate {
Eigen::MatrixXd get_features(const std::vector<Matter> &matobjs);
Eigen::MatrixXd
get_features(const std::vector<std::shared_ptr<Matter>> &matobjs);
Eigen::MatrixXd get_targets(std::vector<std::shared_ptr<Matter>> &matobjs,
                            std::shared_ptr<Potential> true_pot);
Eigen::MatrixXd get_targets(std::vector<Matter> &matobjs,
                            std::shared_ptr<Potential> true_pot);
Eigen::VectorXd make_target(Matter &m1, std::shared_ptr<Potential> true_pot);
std::pair<Eigen::VectorXd, Eigen::VectorXd>
getNewDataPoint(const std::vector<std::shared_ptr<Matter>> &matobjs,
                std::shared_ptr<Potential> true_pot);
std::vector<Matter> getMidSlice(const std::vector<Matter> &matobjs);
bool accuratePES(std::vector<std::shared_ptr<Matter>> &matobjs,
                 std::shared_ptr<Potential> true_pot);
std::pair<double, Eigen::VectorXd::Index>
getMaxUncertainty(const std::vector<std::shared_ptr<Matter>> &matobjs);
} // namespace helper_functions::surrogate

namespace helper_functions::eigen {
Eigen::MatrixXd vertCat(const Eigen::MatrixXd &m1, const Eigen::MatrixXd &m2);
void addVectorRow(Eigen::MatrixXd &data, const Eigen::VectorXd &newrow);
// Modifies data
} // namespace helper_functions::eigen
