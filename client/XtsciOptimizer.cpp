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
#include "eon/XtsciOptimizer.h"
#include "eon/PairHessian.h"
#include "eon/XtsciEindir.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <rgmin.h>

#if RGMIN_ABI_VERSION_MINOR < 10
#error "rgmin ABI minor 10 or newer is required for set_periodic"
#endif

namespace {

// Packwood et al., J. Chem. Phys. 144, 164109 (2016), assembled and
// solved inside rgmin.
constexpr const char *kRgminExp = "rgmin_exp";

// rgmin keeps the pair list it is given and measures a missing distance
// without periodic images, so the host sends minimum-image distances.
struct ExpPairs {
  std::vector<size_t> i;
  std::vector<size_t> j;
  std::vector<double> r;
  double r_nn{0.0};
};

ExpPairs exp_pairs(const eonc::ObjectiveFunction &objf,
                   const Eigen::VectorXd &x, double rcut) {
  const auto nat = static_cast<size_t>(x.size() / 3);
  auto dist = [&](size_t a, size_t b) {
    Eigen::Vector3d d = x.segment<3>(3 * a) - x.segment<3>(3 * b);
    objf.minimumImage(d);
    return d.norm();
  };
  std::vector<double> nn(nat, std::numeric_limits<double>::infinity());
  for (size_t a = 0; a < nat; ++a) {
    for (size_t b = a + 1; b < nat; ++b) {
      const double r = dist(a, b);
      nn[a] = std::min(nn[a], r);
      nn[b] = std::min(nn[b], r);
    }
  }
  ExpPairs out;
  for (const double v : nn) {
    if (std::isfinite(v)) {
      out.r_nn = std::max(out.r_nn, v);
    }
  }
  // ASE and Packwood et al. keep pairs within 2 r_nn.
  const double cut = rcut > 0.0 ? rcut : 2.0 * out.r_nn;
  for (size_t a = 0; a < nat; ++a) {
    for (size_t b = a + 1; b < nat; ++b) {
      const double r = dist(a, b);
      if (r > 0.0 && r < cut) {
        out.i.push_back(a);
        out.j.push_back(b);
        out.r.push_back(r);
      }
    }
  }
  return out;
}

struct XtsObjectiveContext {
  eonc::ObjectiveFunction *objective;
  eonc::PotType pot;
  std::string precon;
  double precon_A;
  double precon_mu;
  double precon_rcut;
  Eigen::VectorXd *cached_x;
  eonc::xtsci_eindir::State *eindir;
};

Eigen::Map<const Eigen::VectorXd>
map_input(const DLManagedTensorVersioned *tensor) {
  const auto &dl = tensor->dl_tensor;
  if (dl.ndim != 1 || dl.dtype.code != kDLFloat || dl.dtype.bits != 64 ||
      dl.dtype.lanes != 1 || dl.device.device_type != kDLCPU ||
      dl.shape == nullptr || dl.data == nullptr) {
    throw std::runtime_error("rgmin objective requires a CPU f64 vector");
  }
  return {static_cast<const double *>(dl.data) +
              dl.byte_offset / sizeof(double),
          static_cast<Eigen::Index>(dl.shape[0])};
}

Eigen::Map<Eigen::VectorXd> map_output(DLManagedTensorVersioned *tensor) {
  const auto &dl = tensor->dl_tensor;
  if (dl.ndim != 1 || dl.dtype.code != kDLFloat || dl.dtype.bits != 64 ||
      dl.dtype.lanes != 1 || dl.device.device_type != kDLCPU ||
      dl.shape == nullptr || dl.data == nullptr) {
    throw std::runtime_error("rgmin tensor requires a CPU f64 vector");
  }
  return {static_cast<double *>(dl.data) + dl.byte_offset / sizeof(double),
          static_cast<Eigen::Index>(dl.shape[0])};
}

// Matter::setPositions always dirties the PES cache. Energy and
// forces come from one potential->force() call, so eval then grad at
// the same x must not setPositions twice.
void set_positions_if_changed(eonc::ObjectiveFunction *obj,
                              const Eigen::VectorXd &x,
                              Eigen::VectorXd *cached) {
  if (cached != nullptr && cached->size() == x.size() &&
      (*cached - x).isZero(0.0)) {
    return;
  }
  // First host iteration: relaxMatter already called isConverged()
  // at this geometry. Do not dirty that cache.
  const auto cur = obj->getPositions();
  if (cur.size() == x.size() && (cur - x).isZero(0.0)) {
    if (cached != nullptr) {
      *cached = x;
    }
    return;
  }
  obj->setPositions(x);
  if (cached != nullptr) {
    *cached = x;
  }
}

// One Matter::computePotential. Energy and forces share that call.
rgmin_status_t evaluate_gradient(void *user, const DLManagedTensorVersioned *x,
                                 double *value_out,
                                 DLManagedTensorVersioned *gradient_out) {
  try {
    auto *context = static_cast<XtsObjectiveContext *>(user);
    if (context->eindir != nullptr) {
      return static_cast<rgmin_status_t>(eonc::xtsci_eindir::eval_grad(
          context->eindir, x, value_out, gradient_out));
    }
    const auto positions = map_input(x);
    set_positions_if_changed(context->objective, positions, context->cached_x);
    *value_out = context->objective->getEnergy();
    const auto gradient_value = context->objective->getGradient();
    auto output = map_output(gradient_out);
    if (output.size() != gradient_value.size()) {
      return RGMIN_INVALID_PARAMETER;
    }
    output = gradient_value;
    return RGMIN_SUCCESS;
  } catch (...) {
    return RGMIN_INTERNAL_ERROR;
  }
}

rgmin_status_t hessian(void *user, const DLManagedTensorVersioned *x,
                       DLManagedTensorVersioned *hess_out) {
  try {
    auto *context = static_cast<XtsObjectiveContext *>(user);
    const auto positions = map_input(x);
    const Eigen::MatrixXd H = eonc::pairhess::build(
        positions, context->precon, context->pot, context->precon_A,
        context->precon_mu, context->precon_rcut, *context->objective);
    auto output = map_output(hess_out);
    const Eigen::Index n = H.rows();
    if (output.size() != n * n) {
      return RGMIN_INVALID_PARAMETER;
    }
    for (Eigen::Index i = 0; i < n; ++i) {
      for (Eigen::Index j = 0; j < n; ++j) {
        output[i * n + j] = H(i, j);
      }
    }
    return RGMIN_SUCCESS;
  } catch (...) {
    return RGMIN_INTERNAL_ERROR;
  }
}

rgmin_method_t method_from_name(const std::string &name) {
  if (name == "lbfgs") {
    return RGMIN_LBFGS;
  }
  if (name == "bfgs") {
    return RGMIN_BFGS;
  }
  if (name == "sr1") {
    return RGMIN_SR1;
  }
  if (name == "sr2") {
    return RGMIN_SR2;
  }
  if (name == "newton") {
    return RGMIN_NEWTON;
  }
  if (name == "rfo") {
    return RGMIN_RFO;
  }
  if (name == "steepest") {
    return RGMIN_STEEPEST;
  }
  if (name == "adam") {
    return RGMIN_ADAM;
  }
  if (name == "pso") {
    return RGMIN_PSO;
  }
  if (name == "polak_ribiere" || name == "nlcg" || name == "pr") {
    return RGMIN_POLAK_RIBIERE;
  }
  if (name == "fletcher_reeves") {
    return RGMIN_FLETCHER_REEVES;
  }
  if (name == "hestenes_stiefel") {
    return RGMIN_HESTENES_STIEFEL;
  }
  if (name == "dai_yuan") {
    return RGMIN_DAI_YUAN;
  }
  if (name == "conjugate_descent") {
    return RGMIN_CONJUGATE_DESCENT;
  }
  if (name == "hager_zhang") {
    return RGMIN_HAGER_ZHANG;
  }
  if (name == "liu_storey") {
    return RGMIN_LIU_STOREY;
  }
  if (name == "fr_pr") {
    return RGMIN_FR_PR;
  }
  if (name == "fire") {
    return RGMIN_FIRE;
  }
  if (name == "bb" || name == "barzilai_borwein") {
    return RGMIN_BB;
  }
  if (name == "dogleg") {
    return RGMIN_DOGLEG;
  }
  if (name == "fire2") {
    return RGMIN_FIRE2;
  }
#if RGMIN_ABI_VERSION_MINOR >= 29
  if (name == "quickmin") {
    return RGMIN_QUICKMIN;
  }
#endif
  throw std::invalid_argument(
      "unknown Xtsci.method '" + name +
      "' (lbfgs, bfgs, sr1, sr2, newton, rfo, steepest, adam, pso, "
      "polak_ribiere, fletcher_reeves, hestenes_stiefel, dai_yuan, "
      "conjugate_descent, hager_zhang, liu_storey, fr_pr, fire, bb, "
      "dogleg, fire2, quickmin)");
}

bool is_host_precon(const std::string &name) {
  return name == "pair" || name == "pair_abs" || name == "pair_full" ||
         name == "exp" || name == "c1" || name == "lindh" ||
         name == "lindh_full" || name == "fischer" || name == "schlegel" ||
         name == "swart";
}

// [Xtsci] qn_step / precon are the natural knobs. [LBFGS] lbfgs_step
// and lbfgs_precon stay on the native optimizer and fill in when the
// Xtsci fields are still at their defaults.
std::string
resolved_qn_step(const eonc::Parameters::optimizer_options_t &opts) {
  if (opts.xtsci.qn_step != "lbfgs") {
    return opts.xtsci.qn_step;
  }
  if (opts.lbfgs.step == "newton" || opts.lbfgs.step == "rfo") {
    return opts.lbfgs.step;
  }
  return opts.xtsci.qn_step;
}

std::string resolved_precon(const eonc::Parameters::optimizer_options_t &opts) {
  if (opts.xtsci.precon != "none") {
    return opts.xtsci.precon;
  }
  return opts.lbfgs.precon;
}

std::string resolved_accept(const eonc::Parameters::optimizer_options_t &opts) {
  if (opts.xtsci.accept != "none") {
    return opts.xtsci.accept;
  }
  if (opts.lbfgs.accept == "energy" || opts.lbfgs.accept == "nonmonotone") {
    return opts.lbfgs.accept;
  }
  return opts.xtsci.accept;
}

} // namespace

namespace eonc {

XtsciOptimizer::~XtsciOptimizer() {
  rgmin_solver_free(m_solver);
  xtsci_eindir::release(m_eindir);
  m_eindir = nullptr;
}

void XtsciOptimizer::ensureSolver(double a_maxMove) {
  if (m_solver != nullptr) {
    return;
  }
  const auto dim = static_cast<size_t>(m_objf->degreesOfFreedom());
  if (dim == 0) {
    return;
  }
  const auto stamp = rgmin_abi_stamp();
  if (rgmin_abi_compatible(&stamp) == 0) {
    throw std::runtime_error("incompatible rgmin ABI");
  }
  const auto &lbfgs = m_optConfig.opts.lbfgs;
  const rgmin_method_t method = method_from_name(m_method);
  double istep = std::max(a_maxMove, 1.0e-12);
  bool dynamics = method == RGMIN_FIRE || method == RGMIN_FIRE2;
#if RGMIN_ABI_VERSION_MINOR >= 29
  dynamics = dynamics || method == RGMIN_QUICKMIN;
#endif
  if (dynamics) {
    istep = m_optConfig.opts.time_step;
    if (istep <= 0.0) {
      istep = m_optConfig.opts.time_step_input;
    }
    if (istep <= 0.0) {
      istep = 0.1;
    }
  }
  rgmin_control_t control{
      static_cast<size_t>(std::max<long>(1, m_optConfig.opts.max_iterations)),
      m_optConfig.opts.converged_force,
      istep,
      static_cast<size_t>(std::max<long>(0, lbfgs.memory)),
      0.0,
  };
  m_solver = rgmin_solver_create(method, &control, dim);
  if (m_solver == nullptr) {
    throw std::runtime_error(rgmin_last_error());
  }
  const auto step = resolved_qn_step(m_optConfig.opts);
  if (step == "newton") {
    rgmin_solver_set_qn_step(m_solver, RGMIN_QN_NEWTON);
  } else if (step == "rfo") {
    rgmin_solver_set_qn_step(m_solver, RGMIN_QN_RFO);
  } else {
    rgmin_solver_set_qn_step(m_solver, RGMIN_QN_LBFGS);
  }
  const auto accept = m_objf->gradientMatchesEnergy()
                          ? resolved_accept(m_optConfig.opts)
                          : std::string("none");
  if (accept == "energy") {
    rgmin_solver_set_accept(m_solver, RGMIN_ACCEPT_ENERGY);
  } else if (accept == "nonmonotone") {
    rgmin_solver_set_accept(m_solver, RGMIN_ACCEPT_NONMONOTONE);
  } else {
    rgmin_solver_set_accept(m_solver, RGMIN_ACCEPT_NONE);
  }
  rgmin_solver_set_project_rigid(m_solver, lbfgs.project_rigid ? 1 : 0);
  rgmin_solver_set_extra_updates(
      m_solver, static_cast<size_t>(std::max<long>(0, lbfgs.extra_updates)));
  if (lbfgs.curvature == "cautious") {
    rgmin_solver_set_cautious(m_solver, lbfgs.cautious_eps,
                              lbfgs.cautious_alpha);
  } else {
    rgmin_solver_set_cautious(m_solver, 0.0, lbfgs.cautious_alpha);
  }
  if (m_optConfig.opts.xtsci.highs) {
    if (rgmin_solver_set_highs(m_solver, 1) != 0) {
      throw std::runtime_error(
          "Xtsci.highs needs rgmin built with --features highs");
    }
  }
  const auto &mani = m_optConfig.opts.xtsci.manifold;
  if (mani == "so3" && dim != 9) {
    throw std::runtime_error(
        "Xtsci.manifold=so3 needs length 9; a 3N cluster is "
        "rigid_quotient (Sella R^{3N}/SE(3)), not a packed rotation");
  }
  if (mani == "se3" && dim != 12) {
    throw std::runtime_error(
        "Xtsci.manifold=se3 needs length 12; a 3N cluster is "
        "rigid_quotient (Sella R^{3N}/SE(3)), not an SE(3) prefix");
  }
  if ((mani == "rigid_quotient" || mani == "mw_rigid" || mani == "sella" ||
       mani == "eckart" || mani == "irc") &&
      (dim < 6 || dim % 3 != 0)) {
    throw std::runtime_error("Xtsci.manifold=" + mani +
                             " needs 3N Cartesians with N >= 2");
  }
  if (mani == "sphere") {
    rgmin_solver_set_manifold(m_solver, RGMIN_MANIFOLD_SPHERE);
  } else if (mani == "so3") {
    rgmin_solver_set_manifold(m_solver, RGMIN_MANIFOLD_SO3);
  } else if (mani == "stiefel") {
    rgmin_solver_set_manifold(m_solver, RGMIN_MANIFOLD_STIEFEL);
  } else if (mani == "se3") {
    rgmin_solver_set_manifold(m_solver, RGMIN_MANIFOLD_SE3);
  } else if (mani == "rigid_quotient" || mani == "sella") {
    rgmin_solver_set_manifold(m_solver, RGMIN_MANIFOLD_RIGID_QUOTIENT);
  } else if (mani == "mw_rigid" || mani == "eckart" || mani == "irc") {
    rgmin_solver_set_manifold(m_solver, RGMIN_MANIFOLD_MW_RIGID);
    const auto masses = m_objf->getMasses();
    if (masses.size() > 0) {
      rgmin_solver_set_masses(m_solver, masses.data(),
                              static_cast<size_t>(masses.size()));
    }
  } else {
    rgmin_solver_set_manifold(m_solver, RGMIN_MANIFOLD_EUCLIDEAN);
  }
  if (mani == "rigid_quotient" || mani == "mw_rigid" || mani == "sella" ||
      mani == "eckart" || mani == "irc") {
    rgmin_solver_set_periodic(m_solver, m_objf->getPeriodic() ? 1 : 0);
  }
  if (resolved_precon(m_optConfig.opts) == kRgminExp &&
      (method == RGMIN_NEWTON || method == RGMIN_RFO ||
       method == RGMIN_DOGLEG || step != "lbfgs")) {
    throw std::invalid_argument(
        "Xtsci.precon = rgmin_exp is a first-order metric. Newton, RFO "
        "and dogleg steps read a host Hessian (pair, exp, lindh, ...)");
  }
  if (m_eindir == nullptr) {
    m_eindir = xtsci_eindir::bind(m_objf.get(), &m_cached_x);
  }
}

// A new pair list on the first step and whenever an atom has moved half
// of r_nn since the last one, the rebuild rule of Packwood et al. The
// estimated mu is carried over, so a rebuild costs no force call.
void XtsciOptimizer::refreshExpPrecon() {
#if RGMIN_ABI_VERSION_MINOR >= 32
  if (m_preconRnn > 0.0) {
    const Eigen::VectorXd d = m_objf->difference(m_x, m_preconAt);
    double moved = 0.0;
    for (Eigen::Index a = 0; a + 2 < d.size(); a += 3) {
      moved = std::max(moved, d.segment<3>(a).norm());
    }
    if (moved <= 0.5 * m_preconRnn) {
      return;
    }
  }
  const auto &lbfgs = m_optConfig.opts.lbfgs;
  const ExpPairs pairs = exp_pairs(*m_objf, m_x, lbfgs.precon_rcut);
  if (pairs.i.empty()) {
    return;
  }
  rgmin_precon_t p;
  rgmin_precon_init(&p, RGMIN_PRECON_EXP);
  p.exp_a = lbfgs.precon_A;
  p.mu = rgmin_solver_precon_scale(m_solver);
  p.nbr_i = pairs.i.data();
  p.nbr_j = pairs.j.data();
  p.nbr_r = pairs.r.data();
  p.npairs = pairs.i.size();
  if (rgmin_solver_set_precon(m_solver, &p) != RGMIN_SUCCESS) {
    throw std::runtime_error(rgmin_last_error());
  }
  m_preconAt = m_x;
  m_preconRnn = pairs.r_nn;
#else
  throw std::runtime_error(
      "Xtsci.precon = rgmin_exp needs rgmin ABI minor 32 or newer");
#endif
}

int XtsciOptimizer::step(double a_maxMove) {
  if (m_objf->degreesOfFreedom() <= 0) {
    return m_objf->isConverged() ? 1 : 0;
  }
  ensureSolver(a_maxMove);
  if (m_solver == nullptr) {
    return m_objf->isConverged() ? 1 : 0;
  }
  // Native LBFGS clips by the largest per-atom move, not ||d||_2.
  rgmin_solver_set_maxmove(m_solver, 0.0);
  rgmin_solver_set_atom_maxmove(m_solver, std::max(a_maxMove, 0.0));

  if (m_x.size() == 0) {
    m_x = m_objf->getPositions();
  }
  if (m_x.size() != m_objf->degreesOfFreedom()) {
    throw std::runtime_error("rgmin objective position dimension mismatch");
  }
  if (resolved_precon(m_optConfig.opts) == kRgminExp) {
    refreshExpPrecon();
  }
  auto *tensor =
      rgmin_tensor_borrow_cpu_f64(m_x.data(), static_cast<size_t>(m_x.size()));
  if (tensor == nullptr) {
    throw std::runtime_error("could not allocate rgmin objective tensor");
  }
  const auto &lbfgs = m_optConfig.opts.lbfgs;
  std::string precon = resolved_precon(m_optConfig.opts);
  const rgmin_method_t method = method_from_name(m_method);
  if (method == RGMIN_DOGLEG && !is_host_precon(precon)) {
    precon = "pair";
  }
  XtsObjectiveContext context{
      m_objf.get(),    m_optConfig.potential, precon,      lbfgs.precon_A,
      lbfgs.precon_mu, lbfgs.precon_rcut,     &m_cached_x, m_eindir,
  };
  rgmin_report_t report{};
  const bool want_hess = method == RGMIN_NEWTON || method == RGMIN_RFO ||
                         method == RGMIN_DOGLEG || is_host_precon(precon);
  rgmin_status_t status;
  if (want_hess) {
    status = rgmin_solver_step_hess_fg(m_solver, evaluate_gradient, hessian,
                                       &context, tensor, &report);
  } else {
    status = rgmin_solver_step_fg(m_solver, evaluate_gradient, &context, tensor,
                                  &report);
  }
  rgmin_tensor_free(tensor);
  if (status != RGMIN_SUCCESS) {
    throw std::runtime_error(rgmin_last_error());
  }
  return m_objf->isConverged() ? 1 : 0;
}

int XtsciOptimizer::run(size_t a_maxIterations, double a_maxMove) {
  if (m_objf->degreesOfFreedom() <= 0 || a_maxIterations == 0) {
    return m_objf->isConverged() ? 1 : 0;
  }
  ensureSolver(a_maxMove);
  const auto &lbfgs = m_optConfig.opts.lbfgs;
  std::string precon = resolved_precon(m_optConfig.opts);
  const rgmin_method_t method = method_from_name(m_method);
  if (method == RGMIN_DOGLEG && !is_host_precon(precon)) {
    precon = "pair";
  }
  const bool want_hess = method == RGMIN_NEWTON || method == RGMIN_RFO ||
                         method == RGMIN_DOGLEG || is_host_precon(precon);
  // A full minimize without a host Hessian goes through the eindir entry.
  // step() stays one session iteration so the host loop keeps its memory.
  if (m_eindir != nullptr && !want_hess && precon != kRgminExp) {
    if (m_x.size() == 0) {
      m_x = m_objf->getPositions();
    }
    double value = 0.0;
    const int status = xtsci_eindir::minimize(
        m_eindir, m_x.data(), static_cast<size_t>(m_x.size()), a_maxIterations,
        m_optConfig.opts.converged_force, std::max(a_maxMove, 1.0e-12),
        static_cast<size_t>(std::max<long>(0, lbfgs.memory)),
        static_cast<int>(method), &value);
    if (status != RGMIN_SUCCESS) {
      throw std::runtime_error(rgmin_last_error());
    }
    m_objf->setPositions(m_x);
    m_cached_x = m_x;
    static_cast<void>(value);
    return m_objf->isConverged() ? 1 : 0;
  }
  for (size_t i = 0; i < a_maxIterations; ++i) {
    if (step(a_maxMove) != 0) {
      return 1;
    }
  }
  return m_objf->isConverged() ? 1 : 0;
}

} // namespace eonc
