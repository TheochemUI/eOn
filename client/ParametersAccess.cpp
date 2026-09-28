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
#include "ParametersImpl.h"

namespace eonc {

const Parameters::potential_options_t &Parameters::potential_options() const {
  return impl_->potential_options_;
}

} // namespace eonc
