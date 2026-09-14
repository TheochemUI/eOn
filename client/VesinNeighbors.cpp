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

} // namespace eonc
