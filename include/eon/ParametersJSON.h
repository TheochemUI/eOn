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

#include <nlohmann/json.hpp>
#include <string>

namespace eonc {

class Parameters;

namespace config {

/// Serialize all Parameters to JSON.
nlohmann::json to_json(const Parameters &params);

/// Deserialize JSON into Parameters, then resolve computed fields.
void from_json(const nlohmann::json &j, Parameters &params);

/// Load Parameters from a JSON string. Returns 0 on success.
int load_json(const std::string &json_str, Parameters &params);

} // namespace config
} // namespace eonc
