#pragma once

#include "Job.h"
#include "Parameters.h"
#include "HelperFunctions.h"

#include <pybind11/embed.h>
#include <pybind11/eigen.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include "potentials/PySurrogate/PySurrogate.h"

#include "NudgedElasticBand.h"

class GPSurrogateJob : public Job {
public:
  GPSurrogateJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~GPSurrogateJob(void) = default;
  std::vector<std::string> run(void) override;

private:
    pybind11::scoped_interpreter guard{};
};

namespace helper_functions::surrogate {
  Eigen::MatrixXd get_features(const std::vector<Matter>& matobjs);
  Eigen::MatrixXd get_features(const std::vector<std::shared_ptr<Matter>>& matobjs);
  Eigen::MatrixXd get_targets(std::vector<std::shared_ptr<Matter>>& matobjs);
  Eigen::MatrixXd get_targets(std::vector<Matter>& matobjs);
  std::vector<Matter> getMidSlice(const std::vector<Matter>& matobjs);
}

namespace helper_functions::eigen {
  Eigen::MatrixXd vertCat(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2);
}
