/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/

#include "eon/VesinNeighbors.h"

#include <cstdlib>
#include <cstring>

namespace eonc {

VesinNeighbors::VesinNeighbors(VesinNeighbors &&other) noexcept
    : list_(other.list_),
      owns_(other.owns_) {
  other.list_ = VesinNeighborList{};
  other.owns_ = false;
}

VesinNeighbors &VesinNeighbors::operator=(VesinNeighbors &&other) noexcept {
  if (this != &other) {
    free_list();
    list_ = other.list_;
    owns_ = other.owns_;
    other.list_ = VesinNeighborList{};
    other.owns_ = false;
  }
  return *this;
}

VesinNeighbors::~VesinNeighbors() { free_list(); }

void VesinNeighbors::free_list() {
  if (owns_) {
    vesin_free(&list_);
    owns_ = false;
  }
  list_ = VesinNeighborList{};
}

void VesinNeighbors::compute(const double *R, std::size_t n, const double *box,
                             const Options &opt) {
  if (R == nullptr || (n > 0 && box == nullptr)) {
    throw std::invalid_argument("VesinNeighbors::compute: null R or box");
  }
  if (opt.cutoff <= 0.0) {
    free_list();
    return;
  }

  // Pass the existing VesinNeighborList through so vesin can re-use pair /
  // distance / vector buffers (upstream contract). Do not vesin_free first —
  // that forces a full realloc every call and is the #386 ASV regression.
  VesinOptions vopt{};
  vopt.cutoff = opt.cutoff;
  vopt.full = opt.full;
  vopt.sorted = opt.sorted;
  vopt.algorithm = opt.algorithm;
  vopt.skin = opt.skin;
  vopt.n_threads = opt.n_threads;
  vopt.return_shifts = opt.return_shifts;
  vopt.return_distances = opt.return_distances;
  vopt.return_vectors = opt.return_vectors;

  bool periodic[3] = {opt.periodic[0], opt.periodic[1], opt.periodic[2]};
  double box33[3][3];
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      box33[a][b] = box[3 * a + b];
    }
  }

  VesinDevice cpu{VesinCPU, 0};
  const char *error_message = nullptr;
  int status =
      vesin_neighbors(reinterpret_cast<const double (*)[3]>(R), n, box33,
                      periodic, cpu, vopt, &list_, &error_message);
  if (status != EXIT_SUCCESS) {
    std::string err = "vesin_neighbors failed";
    if (error_message != nullptr) {
      err += ": ";
      err += error_message;
    }
    free_list();
    throw std::runtime_error(err);
  }
  owns_ = true;
}

VesinNeighborList *VesinNeighbors::release() {
  if (!owns_) {
    return nullptr;
  }
  auto *heap = new VesinNeighborList(list_);
  list_ = VesinNeighborList{};
  owns_ = false;
  return heap;
}

bool CachedPairList::valid(const double *R, std::size_t n, const double *box,
                           const Options &opt) const {
  if (!built_ || n != n_ || !(opt == opt_)) {
    return false;
  }
  // Any box change invalidates the precomputed shift offsets.
  for (int k = 0; k < 9; ++k) {
    if (box[k] != boxref_[k]) {
      return false;
    }
  }
  if (complete_) {
    return true; // all pairs are candidates; motion cannot invalidate
  }
  // Verlet criterion: every atom within skin/2 of its build position. Early
  // exit on the first violator (a different NEB image rejects on atom ~1).
  const double thr2 = 0.25 * opt_.skin * opt_.skin;
  const double *ref = Rref_.data();
  for (std::size_t a = 0; a < n; ++a) {
    const double dx = R[3 * a] - ref[3 * a];
    const double dy = R[3 * a + 1] - ref[3 * a + 1];
    const double dz = R[3 * a + 2] - ref[3 * a + 2];
    if (dx * dx + dy * dy + dz * dz > thr2) {
      return false;
    }
  }
  return true;
}

void CachedPairList::rebuild(const double *R, std::size_t n, const double *box,
                             const Options &opt) {
  setup(n, box, opt);
  if (mic_) {
    double w[3];
    double inv[3];
    for (int k = 0; k < 3; ++k) {
      w[k] = box[4 * k];
      inv[k] = (opt.periodic[static_cast<std::size_t>(k)] && w[k] != 0.0)
                   ? 1.0 / w[k]
                   : 0.0;
    }
    const double bc = opt.cutoff + opt.skin;
    vesin::cpu::brute_force_visit(
        R, n, w, inv, bc * bc, -1.0, pairsIJ_,
        [](int32_t, int32_t, double, double, double, double) {});
    finishRebuild(R, n, box, opt);
    return;
  }
  rebuildCell(R, n, box, opt);
}

void CachedPairList::setup(std::size_t n, const double *box,
                           const Options &opt) {
  // MIC mode: orthorhombic box with the true cutoff inside half the
  // smallest periodic width. Brute-force nearest-image candidates at
  // cutoff+skin form a valid superset (MIC distance is 1-Lipschitz in the
  // displacements), and forEach re-folds per call, so image identity is
  // resolved at evaluation time exactly like the historical MIC loops. The
  // cell list degenerates to a couple of cells in this regime and costs
  // ~20x more per build.
  const bool orthorhombic = box[1] == 0.0 && box[2] == 0.0 && box[3] == 0.0 &&
                            box[5] == 0.0 && box[6] == 0.0 && box[7] == 0.0;
  mic_ = orthorhombic && n <= 20000;
  for (int k = 0; mic_ && k < 3; ++k) {
    if (opt.periodic[static_cast<std::size_t>(k)] &&
        opt.cutoff > 0.5 * box[4 * k]) {
      mic_ = false;
    }
  }
}

void CachedPairList::rebuildCell(const double *R, std::size_t n,
                                 const double *box, const Options &opt) {
  VesinNeighbors::Options vopt;
  vopt.cutoff = opt.cutoff + opt.skin;
  vopt.full = false;
  vopt.sorted = false; // vesin's permutation sort costs more than it saves
  vopt.algorithm = VesinAutoAlgorithm;
  vopt.return_shifts = true;
  vopt.return_distances = false;
  vopt.return_vectors = false;
  vopt.periodic = opt.periodic;
  nl_.compute(R, n, box, vopt);
  pairsIJ_.clear();
  finishRebuild(R, n, box, opt);
}

void CachedPairList::finishRebuild(const double *R, std::size_t n,
                                   const double *box, const Options &opt) {
  for (int k = 0; k < 3; ++k) {
    micInv_[static_cast<std::size_t>(k)] =
        (mic_ && opt.periodic[static_cast<std::size_t>(k)]) ? 1.0 / box[4 * k]
                                                            : 0.0;
  }

  Rref_.assign(R, R + 3 * n);
  for (int k = 0; k < 9; ++k) {
    boxref_[k] = box[k];
  }
  opt_ = opt;
  n_ = n;
  // Complete graph (small cluster inside the cutoff): with MIC evaluation
  // the candidate set can never lose a pair, so the slot stays valid for
  // arbitrary motion and rebuilds stop entirely.
  complete_ = mic_ && pairsIJ_.size() / 2 == n * (n - 1) / 2;
  built_ = true;
}

std::shared_ptr<const CachedPairList>
PairListCache::ensure(const double *R, std::size_t n, const double *box,
                      const CachedPairList::Options &opt) {
  {
    std::lock_guard<std::mutex> lock(mu_);
    for (std::size_t s = 0; s < slots_.size(); ++s) {
      if (slots_[s]->valid(R, n, box, opt)) {
        if (s != 0) {
          auto slot = std::move(slots_[s]);
          slots_.erase(slots_.begin() + static_cast<std::ptrdiff_t>(s));
          slots_.insert(slots_.begin(), std::move(slot));
        }
        return slots_.front();
      }
    }
  }

  // Build outside the pool lock; concurrent misses on the same geometry cost
  // one duplicate build, never a wrong result.
  auto fresh = std::make_shared<CachedPairList>();
  fresh->rebuild(R, n, box, opt);

  std::lock_guard<std::mutex> lock(mu_);
  slots_.insert(slots_.begin(), fresh);
  if (slots_.size() > kMaxSlots) {
    slots_.pop_back(); // readers holding the shared_ptr keep it alive
  }
  return fresh;
}

PairListCache &PairListCache::global() {
  static PairListCache cache;
  return cache;
}

} // namespace eonc
