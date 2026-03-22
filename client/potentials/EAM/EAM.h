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
#pragma once

#include <array>
#include <cmath>
#include <vector>

#ifndef EAM_STANDALONE
#include "../../Potential.h"
#endif

/// EAM (Embedded Atom Method) potential with cell list neighbor finding.
class EAM
#ifndef EAM_STANDALONE
    : public Potential
#endif
{
public:
  explicit EAM(const Parameters &params)
      : Potential(PotType::EAM_AL, params), rc_{6.0, 6.0, 6.0} {}

  ~EAM() override = default;

  void cleanMemory();
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *fullbox) override;

private:
  struct element_parameters {
    const int Z;                // Atomic number
    const double Dm;            // Morse potential well depth
    const double alphaM;        // Curvature at Morse minimum
    const double Rm;            // Position of Morse minimum
    const double beta1;         // Density parameter 1
    const double beta2;         // Density parameter 2
    const double r_cut;         // Cutoff distance
    const double func_coeff[9]; // 8th order poly for embedding function
  };
  static const element_parameters el_params[];

  std::vector<long> celllist_old_;
  std::vector<long> celllist_new_;
  std::vector<long> neigh_list_;
  bool initialized_{false};
  std::array<double, 3> rc_;

  void calc_force(long N, double *R, const int *atomicNrs, double *F, double *U,
                  const double *box);
  void new_celllist(long N, const double *box, long *num_axis,
                    long *cell_length, long *celllist_new, long num_cells,
                    double *Rnew);
  void cell_to_neighbor(long N, long num_of_cells, long *num_axis,
                        long *cell_length, long *celllist_new,
                        long *neigh_list);
  int update_cell_list(long N, long num_cells, long *num_axis,
                       long *cell_length, long *celllist_old, double *Rnew);
  [[nodiscard]] static double embedding_function(const double *func_coeff,
                                                 double rho);
  [[nodiscard]] static double embedding_force(const double *func_coeff,
                                              double rho);
  [[nodiscard]] element_parameters get_element_parameters(int atomic_number);
};
