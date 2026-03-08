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
#include "Parameters.h"

namespace eonc {


class GlobalOptimization {
public:
  GlobalOptimization(const Parameters &params);
  ~GlobalOptimization(void);
  void run(void);

private:
  const Parameters &parameters;
};

} // namespace eonc

using eonc::GlobalOptimization;
