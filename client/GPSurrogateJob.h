#pragma once

#include "Job.h"
#include "Parameters.h"
#include "HelperFunctions.h"

class GPSurrogateJob : public Job {
public:
  GPSurrogateJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~GPSurrogateJob(void) = default;
  std::vector<std::string> run(void) override;
};

namespace helper_functions::surrogate {
  Eigen::MatrixXd get_features(const std::vector<Matter>& matobjs);
  Eigen::MatrixXd get_targets(std::vector<Matter>& matobjs);
  std::vector<Matter> getMidSlice(const std::vector<Matter>& matobjs);
}
