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

void CachedPairList::rebuild(VesinNeighbors &scratch, const double *R,
                             std::size_t n, const double *box,
                             const Options &opt) {
  VesinNeighbors::Options vopt;
  vopt.cutoff = opt.cutoff + opt.skin;
  vopt.full = false;
  vopt.sorted = true; // index-sorted pairs: better locality in forEach
  vopt.return_shifts = true;
  vopt.return_distances = false;
  vopt.return_vectors = false;
  vopt.periodic = opt.periodic;
  scratch.compute(R, n, box, vopt);

  const VesinNeighborList &vl = scratch.raw();
  plain_.clear();
  shifted_.clear();
  plain_.reserve(vl.length);
  for (std::size_t p = 0; p < vl.length; ++p) {
    const auto i = static_cast<int32_t>(vl.pairs[p][0]);
    const auto j = static_cast<int32_t>(vl.pairs[p][1]);
    const int32_t sa = vl.shifts[p][0];
    const int32_t sb = vl.shifts[p][1];
    const int32_t sc = vl.shifts[p][2];
    if (sa == 0 && sb == 0 && sc == 0) {
      if (i != j) { // self pairs only arise as periodic images (shifted)
        plain_.push_back({i, j});
      }
    } else {
      // Offset S @ H with row-major cell rows a,b,c.
      const double ox = sa * box[0] + sb * box[3] + sc * box[6];
      const double oy = sa * box[1] + sb * box[4] + sc * box[7];
      const double oz = sa * box[2] + sb * box[5] + sc * box[8];
      shifted_.push_back({i, j, ox, oy, oz});
    }
  }

  Rref_.assign(R, R + 3 * n);
  for (int k = 0; k < 9; ++k) {
    boxref_[k] = box[k];
  }
  opt_ = opt;
  n_ = n;
  built_ = true;
}

const CachedPairList &PairListCache::ensure(const double *R, std::size_t n,
                                            const double *box,
                                            const CachedPairList::Options &opt) {
  for (std::size_t s = 0; s < slots_.size(); ++s) {
    if (slots_[s]->valid(R, n, box, opt)) {
      if (s != 0) {
        auto slot = std::move(slots_[s]);
        slots_.erase(slots_.begin() + static_cast<std::ptrdiff_t>(s));
        slots_.insert(slots_.begin(), std::move(slot));
      }
      return *slots_.front();
    }
  }

  std::unique_ptr<CachedPairList> slot;
  if (slots_.size() < kMaxSlots) {
    slot = std::make_unique<CachedPairList>();
  } else {
    slot = std::move(slots_.back()); // recycle LRU (keeps its allocations)
    slots_.pop_back();
  }
  slot->rebuild(scratch_, R, n, box, opt);
  slots_.insert(slots_.begin(), std::move(slot));
  return *slots_.front();
}

PairListCache &PairListCache::local() {
  thread_local PairListCache cache;
  return cache;
}

} // namespace eonc
