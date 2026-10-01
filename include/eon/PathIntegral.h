/*
** This file is part of eOn.
**
** i-PI Copyright (C) 2014-2015 i-PI developers
**
** Permission is hereby granted, free of charge, to any person obtaining
** a copy of this software and associated documentation files (the
** "Software"), to deal in the Software without restriction, including
** without limitation the rights to use, copy, modify, merge, publish,
** distribute, sublicense, and/or sell copies of the Software, and to
** permit persons to whom the Software is furnished to do so, subject to
** the following conditions:
**
** The above copyright notice and this permission notice shall be
** included in all copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
** MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
** NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
** BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
** ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
** CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
** SOFTWARE.
**
** SPDX-License-Identifier: MIT
*/
#pragma once

#include "eon/Eigen.h"
#include "eon/Potential.h"

#include <cstdint>
#include <string>
#include <vector>

namespace eonc::pathintegral {

/// Trotter springs, or economised springs fitted to harmonic radii of
/// gyration up to a maximum frequency (Zeng and Manolopoulos,
/// arXiv:2607.06414).
enum class Springs { Trotter, Eco };

/// PILE-L, or a normal-mode GLE on the internal modes with a separate
/// Langevin thermostat on the centroid.
enum class Thermostat { Pile, Piglet };

struct Options {
  long beads{8};
  double temperature{1.0};
  double kB{1.0};
  double hbar{1.0};
  double dt{0.005};
  Springs springs{Springs::Trotter};
  Thermostat thermostat{Thermostat::Pile};
  /// Centroid Langevin time, in the same time unit as dt.
  double pileTau{0.2};
  /// Scales the critical PILE damping of the internal modes.
  double pileScale{1.0};
  /// Highest physical frequency the economised springs reproduce.
  double ecoOmegaMax{0.0};
  /// Normal-mode GLE matrices. One block per internal mode, or one
  /// block per bead with the centroid block ignored.
  std::string gleFile;
  std::uint64_t seed{1};
};

struct Sample {
  double kineticCv{0.0};
  /// Average of n · f_centroid. The derivative of the potential of
  /// mean force along the unit normal is the negative of this value.
  double meanForce{0.0};
  long batches{0};
};

/// Economised springs and a normal-mode GLE are refused. The instanton
/// assumes Trotter springs and is refused the same way.
void requireTrotterSprings(const std::string &springs, const char *use);

/// Dimensionless free-ring eigenvalues, mode 0 equal to 0. Physical
/// frequencies are omegan times these values, with
/// omegan = beads * kB * T / hbar.
VectorXd trotterEigenvalues(long nBeads);
VectorXd ecoEigenvalues(long nBeads, double xmax);

/// Orthogonal bead-to-normal-mode matrix. Row k is mode k.
MatrixXd normalModeMatrix(long nBeads);

/// Ring-polymer NVT step. Physical forces on every bead are one
/// forceBatch call, so a calculator group carries the beads together.
/// The classical velocity Verlet step is not used.
class RingPolymer {
public:
  RingPolymer(long nAtoms, std::vector<double> masses,
              std::vector<int> atomicNumbers, std::vector<char> free,
              Options opt);

  void setAllBeads(const double *q);
  /// One position vector of length 3 * nAtoms per bead. The centroid is
  /// projected onto the hyperplane when one is set.
  void setBeads(const std::vector<VectorXd> &beads);
  [[nodiscard]] const std::vector<VectorXd> &beads() const { return q_; }
  void thermalMomenta();

  /// Hold n · (q_centroid - origin) = 0. n is normalised on the free
  /// coordinates. origin has length 3 * nAtoms.
  void setHyperplane(const VectorXd &normal, const VectorXd &origin);

  void step(Potential &pot, const double *box, bool record);

  [[nodiscard]] double kineticCv() const;
  [[nodiscard]] double meanForce() const;
  [[nodiscard]] long batches() const { return batches_; }
  [[nodiscard]] VectorXd centroid() const;
  [[nodiscard]] VectorXd centroidVelocity() const;

  void resetAverages();

  Sample sample(Potential &pot, const double *box, long equilibration,
                long production);

private:
  void forces(Potential &pot, const double *box);
  void thermostat(double h);
  void kick(double h, bool dropParallel);
  void propagate(double h);
  void projectPosition();
  void projectMomentum();
  void toNormal(const std::vector<VectorXd> &src,
                std::vector<VectorXd> &dst) const;
  void fromNormal(const std::vector<VectorXd> &src,
                  std::vector<VectorXd> &dst) const;
  double gauss();
  void initGle();

  Options opt_;
  long nAtoms_{0};
  long nDof_{0};
  long nBeads_{0};
  long nFree_{0};
  std::vector<double> mass_;
  std::vector<int> atomicNumbers_;
  std::vector<char> free_;
  std::vector<long> freeIndex_;
  MatrixXd modes_;
  VectorXd omegaK_;
  std::vector<VectorXd> q_;
  std::vector<VectorXd> p_;
  std::vector<VectorXd> f_;
  std::vector<VectorXd> qnm_;
  std::vector<VectorXd> pnm_;
  VectorXd planeNormal_;
  VectorXd planeOrigin_;
  bool constrain_{false};
  bool haveForces_{false};

  struct ModeGle {
    MatrixXd drift;
    MatrixXd covariance;
    MatrixXd propagate;
    MatrixXd noise;
    MatrixXd extended;
  };
  std::vector<ModeGle> gle_;

  std::uint64_t rng_;
  long batches_{0};
  long recorded_{0};
  double kineticSum_{0.0};
  double forceSum_{0.0};
};

} // namespace eonc::pathintegral
