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
#include "eon/JobResult.h"

#include "eon/Parameters.h"
#include "version.h"

#include <algorithm>
#include <cctype>
#include <string>

#ifdef WITH_XTSCI
#include <xts.h>
#endif

namespace eonc {
namespace {

std::string lower_name(std::string name) {
  std::ranges::transform(name, name.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return name;
}

} // namespace

JobResultProvenance provenanceForJob(const Parameters &params) {
  JobResultProvenance provenance;
  provenance.backend = lower_name(
      std::string(magic_enum::enum_name(params.optimizer_options().method)));
  provenance.engine_version = VERSION;
  provenance.engine_build_identity = GIT_HASH;
  provenance.rgpot_name = lower_name(
      std::string(magic_enum::enum_name(params.potential_options().potential)));
#ifdef WITH_RGPOT
  provenance.rgpot_version = std::string(JobResultProvenance::rgpot_pin);
#endif
#ifdef WITH_XTSCI
  if (params.optimizer_options().method == OptType::XTSCI) {
    provenance.xtsci = true;
    const auto stamp = xts_abi_stamp();
    provenance.xts_abi_major = stamp.abi_major;
    provenance.xts_abi_minor = stamp.abi_minor;
    provenance.xts_abi_layout = stamp.layout_revision;
  }
#endif
  return provenance;
}

} // namespace eonc
