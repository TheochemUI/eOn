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
#include "client/MinimizationJob.hpp"
#include "ConjugateGradients.h"
#include "OptimCreator.hpp"
#include "client/HelperFunctions.h"
#include "client/objectives/MatterObjf.hpp"
#include "parsers/ParseOptim.hpp"

namespace eonc {
bool MinimizationJob::runImpl(Matter &mat) {
  // TODO(rg): params are from toml
  const auto optparams = std::visit(opt::OptBaseVisitor{}, m_p.oparams);
  eonc::objf::MatterObjectiveFunction objf(
      {optparams.optCM, optparams.optConvergedForce}, mat);

  auto optim = std::visit(
      [&](const OptParams &optParam) -> decltype(auto) {
        return eonc::mkOptim{}(objf,
                               std::get<ConjugateGradients::Params>(optParam));
      },
      m_p.oparams);

  // const auto config = toml::table{{"Optimizer", toml::table{{"method",
  // "cg"}}}}; auto optim = helpers::create::mkOptim(objf, config);

  // std::ostringstream min;
  // min << prefixMovie;
  // if (writeMovie) {
  //   matter2con(min.str(), false);
  // }
  bool quiet = false;
  const size_t maxIter = optparams.optMaxIter;
  const double maxMove = optparams.optMaxMove;
  int iteration = 0;
  if (!quiet) {
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10s}  {:14s}  {:18s}  {:13s}", "[Matter]",
                        "Iter", "Step size",
                        magic_enum::enum_name(optparams.optCM), "Energy");
    SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                        "[Matter]", iteration, 0.0, objf.getConvergence(),
                        objf.getEnergy());
  }

  // TODO(rg): Refactor further
  while (!objf.isConverged() && iteration < maxIter) {

    VectorType pos = objf.getPositions();

    optim.step(maxMove);
    iteration++;
    mat.setPositionsFreeV(objf.getPositions());

    double stepSize = helper_functions::maxAtomMotionV(
        mat.pbcV(mat.getPositionsFreeV() - pos));

    if (!quiet) {
      SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                          "[Matter]", iteration, stepSize,
                          objf.getConvergence(), objf.getEnergy());
    }

    // if (writeMovie) {
    //   matter2con(min.str(), true);
    // }

    // if (checkpoint) {
    //   ostringstream chk;
    //   chk << prefixCheckpoint << "_cp";
    //   matter2con(chk.str(), false);
    // }
  }

  if (iteration == 0) {
    if (!quiet) {
      SPDLOG_LOGGER_DEBUG(m_log, "{} {:10}  {:14.5e}  {:18.5e}  {:13.5f}",
                          "[Matter]", iteration, 0.0, objf.getConvergence(),
                          objf.getEnergy());
    }
  }
  //    bool converged = optimizer->run(parameters->optMaxIterations,
  //    parameters->optMaxMove);
  return objf.isConverged();
}

} // namespace eonc
