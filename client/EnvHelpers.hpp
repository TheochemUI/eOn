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
std::string get_value_from_env_or_param(const char *env_variable,
                                        const std::string &param_value,
                                        const std::string &default_value = "",
                                        const std::string &warning_message = "",
                                        const bool is_mandatory = false);
}
