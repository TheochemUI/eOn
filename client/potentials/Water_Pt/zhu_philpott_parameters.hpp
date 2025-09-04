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
/** @file
Parameters for potential ZhuPhilpott.
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
@see zhu_philpott.hpp
*/

namespace forcefields {
/// Parameters for ZhuPhilpott.
/// Parameters for Zhu and Philpott potential for interaction between platinum
/// and water.
namespace zhu_philpott_parameters {
/** Standard parameters for ZhuPhilpott.
This is the parameterisation published in the original paper by Zhu and
Phillpott. The binding energy of one molecule on a (100) surface is 0.50 eV
*/
struct Standard {
  static double const sigmaO_;
  static double const epsilonO_;
  static double const sigmaH_;
  static double const epsilonH_;
  static double const C10_O_;
  static double const C10_H_;
  static double const alpha_;
  static double const sigmaHPt_;
  static double const sigmaOPt_;
  static double const epsilonHPt_;
  static double const epsilonOPt_;
};
/** Iceland parameterisation.
In this parameterisation, the Lennard-Jones interactions between platinum and
water have been weakened to match the binding energy found by DFT calculation
which is 0.3 eV for one molecule on a (100) surface.
*/
struct Iceland {
  static double const sigmaO_;
  static double const epsilonO_;
  static double const sigmaH_;
  static double const epsilonH_;
  static double const C10_O_;
  static double const C10_H_;
  static double const alpha_;
  static double const sigmaHPt_;
  static double const sigmaOPt_;
  static double const epsilonHPt_;
  static double const epsilonOPt_;
};
} // namespace zhu_philpott_parameters
} // namespace forcefields
