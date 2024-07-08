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
#include <string>

namespace helper_functions {

struct PDef {
  std::string param;
  std::string defval;
  PDef(const std::string &param_value, const std::string &default_value)
      : param{param_value},
        defval{default_value} {}
};

std::string get_value_from_env_or_param(const char *env_variable,
                                        const PDef val_def,
                                        const std::string &warning_message = "",
                                        const bool is_mandatory = false);
} // namespace helper_functions
