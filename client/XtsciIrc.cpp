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
#include "eon/XtsciIrc.h"

#include <stdexcept>

#ifdef WITH_RGSADDLE
#include <rgsaddle.h>
#if RGSADDLE_ABI_MINOR >= 14
#define EONC_RGSADDLE_IRC 1
#endif
#endif

#ifdef EONC_RGSADDLE_IRC
#include <algorithm>
#include <deque>
#include <memory>
#include <string>
#endif

namespace eonc {

#ifdef EONC_RGSADDLE_IRC
namespace {

// rgsaddle evaluates the accepted point of a step again when the next step
// starts, and a rejected trial leaves that point a few ulps away.
constexpr double kSamePoint = 1.0e-12;

struct Evaluation {
  Eigen::VectorXd x;
  double energy;
  Eigen::VectorXd gradient;
};

struct IrcUser {
  Matter *probe;
  std::deque<Evaluation> seen;
};

const Evaluation &evaluate(IrcUser &ctx, const double *x, Eigen::Index n) {
  const Eigen::Map<const Eigen::VectorXd> at(x, n);
  for (const auto &e : ctx.seen) {
    if ((e.x - at).cwiseAbs().maxCoeff() <= kSamePoint) {
      return e;
    }
  }
  ctx.probe->setPositionsFreeV(at);
  Evaluation e{at, ctx.probe->getPotentialEnergy(),
               -ctx.probe->getForcesFreeV()};
  if (ctx.seen.size() >= 4) {
    ctx.seen.pop_front();
  }
  ctx.seen.push_back(std::move(e));
  return ctx.seen.back();
}

rgsaddle_status_t surfaceCallback(void *user, rgsaddle_surface_request_t *req) {
  auto *ctx = static_cast<IrcUser *>(user);
  if (ctx == nullptr || req == nullptr || req->positions == nullptr ||
      req->energies == nullptr || req->gradients == nullptr ||
      req->n_images != 1) {
    return RGSADDLE_SURFACE_FAILED;
  }
  if (req->n_atoms != ctx->probe->numberOfFreeAtoms()) {
    return RGSADDLE_SHAPE;
  }
  try {
    const Evaluation &e = evaluate(*ctx, req->positions, 3 * req->n_atoms);
    req->energies[0] = e.energy;
    std::copy(e.gradient.begin(), e.gradient.end(), req->gradients);
    return RGSADDLE_OK;
  } catch (...) {
    return RGSADDLE_SURFACE_FAILED;
  }
}

void checkStatus(int rc, const char *what) {
  if (rc != RGSADDLE_OK) {
    throw std::runtime_error(
        std::string(what) + ": " +
        rgsaddle_status_name(static_cast<rgsaddle_status_t>(rc)));
  }
}

IrcBranch followBranch(RgsaddleIrc *session, IrcUser &ctx, long maxSteps,
                       Eigen::Index n) {
  IrcBranch branch;
  Eigen::VectorXd x(n);
  Eigen::VectorXd last;
  for (long k = 0; k < maxSteps && !branch.reachedMinimum; ++k) {
    rgsaddle_report_t report{};
    checkStatus(rgsaddle_irc_step(session, surfaceCallback, &ctx, &report),
                "rgsaddle_irc_step");
    checkStatus(rgsaddle_irc_position(session, x.data()),
                "rgsaddle_irc_position");
    ++branch.steps;
    // rgsaddle reports a minimum only after the first step and at a largest
    // atomic force at or below config.force_tol, which followIrc sets.
    branch.reachedMinimum = report.status == RGSADDLE_STATUS_CONVERGED;
    if (last.size() != n || (x - last).cwiseAbs().maxCoeff() > kSamePoint) {
      const Evaluation &e = evaluate(ctx, x.data(), n);
      branch.energies.push_back(e.energy);
      branch.maxForces.push_back(AtomMatrix::Map(e.gradient.data(), n / 3, 3)
                                     .rowwise()
                                     .norm()
                                     .maxCoeff());
      ctx.probe->setPositionsFreeV(x);
      branch.positions.push_back(ctx.probe->getPositions());
      last = x;
    }
  }
  return branch;
}

} // namespace
#endif

bool ircAvailable() {
#ifdef EONC_RGSADDLE_IRC
  rgsaddle_version_t abi{};
  return rgsaddle_abi_stamp(&abi) == RGSADDLE_OK &&
         abi.major == RGSADDLE_ABI_MAJOR && abi.minor >= 14;
#else
  return false;
#endif
}

IrcPath followIrc(const Matter &saddle, const AtomMatrix &mode,
                  const Parameters &params) {
  if (!ircAvailable()) {
    throw std::runtime_error("rgsaddle IRC is not available in this build");
  }
#ifdef EONC_RGSADDLE_IRC
  const auto &opts = params.saddle_search_options().irc;
  std::vector<double> masses;
  std::vector<long> rows;
  for (long i = 0; i < saddle.numberOfAtoms(); ++i) {
    if (!saddle.getFixed(i)) {
      masses.push_back(saddle.getMass(i));
      rows.push_back(i);
    }
  }
  if (rows.empty() || mode.rows() != saddle.numberOfAtoms()) {
    throw std::runtime_error("IRC needs free atoms and one mode row per atom");
  }
  const auto nFree = static_cast<Eigen::Index>(rows.size());
  AtomMatrix freeMode(nFree, 3);
  for (Eigen::Index j = 0; j < nFree; ++j) {
    freeMode.row(j) = mode.row(rows[static_cast<size_t>(j)]);
  }
  const Eigen::VectorXd x0 = saddle.getPositionsFreeV();

  Matter probe(saddle);
  IrcUser ctx{&probe, {}};
  // One session per branch: a GS2 session halves dx on every rejected trial
  // and keeps the shorter step through rgsaddle_irc_set_direction.
  auto branch = [&](int32_t direction) {
    rgsaddle_irc_config_t config{};
    config.version = RGSADDLE_VERSION_INIT;
    config.dx = opts.step_size;
    config.force_tol = opts.converged_force;
    config.force_gate = RGSADDLE_FORCE_MAX_ATOM;
    config.max_move = params.optimizer_options().max_move;
    config.max_inner = 10;
    config.kind = RGSADDLE_IRC_GS2;
    config.direction = direction;
    std::unique_ptr<RgsaddleIrc, decltype(&rgsaddle_irc_free)> session(
        rgsaddle_irc_create(&config, nFree, x0.data(), masses.data(),
                            freeMode.data()),
        &rgsaddle_irc_free);
    if (!session) {
      throw std::runtime_error("rgsaddle_irc_create failed");
    }
    return followBranch(session.get(), ctx, opts.max_steps, 3 * nFree);
  };
  IrcBranch forward = branch(RGSADDLE_IRC_FORWARD);
  IrcBranch reverse = branch(RGSADDLE_IRC_REVERSE);

  // The forward kick takes its sign from rgsaddle.
  const bool forwardAlongMode =
      forward.positions.empty() ||
      saddle.pbc(forward.positions.front() - saddle.getPositions())
              .cwiseProduct(mode)
              .sum() > 0.0;
  IrcPath path;
  path.plus = forwardAlongMode ? std::move(forward) : std::move(reverse);
  path.minus = forwardAlongMode ? std::move(reverse) : std::move(forward);
  return path;
#else
  static_cast<void>(saddle);
  static_cast<void>(mode);
  static_cast<void>(params);
  return {};
#endif
}

} // namespace eonc
