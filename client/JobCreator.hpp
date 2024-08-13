#pragma once

#include "BaseStructures.h"
#include "Job.h"
#include "Parser.hpp"
#include "PointJob.h"
#include "thirdparty/toml.hpp"

namespace eonc {

/*
** Do the makers of Jobs need Matter objects? Perhaps not..
** In fact, the Job should just take a job, other children (single, double ended
*etc. can have more / less methods of inputs..)
** Unlike Potentials, generated Jobs are to be unique
** Put another way, Job arity in terms of how many matter objects are required
*is a design constraint..
*/
// Forwarding Matter objects
template <typename... Args>
std::unique_ptr<JobBase> makeJob(const toml::table &config, Args &&...args) {
  config_section(config, "Main");
  auto jtype = get_enum_toml<JobType>(config["Main"]["job"]);

  switch (jtype) {
  case JobType::Point: {
    return std::make_unique<PointJob>(std::forward<Args>(args)...);
  }
  // Add other cases for different JobTypes here
  default: {
    throw std::runtime_error("No known job could be constructed");
  }
  }
}

} // namespace eonc
