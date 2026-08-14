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

#include <cstdlib>
#include <stdexcept>
#include <string>

#ifndef _WIN32
#include <sys/wait.h>
#endif

namespace eonc::pot {

/// \brief Describe how a command run through system() ended.
///
/// \param status The value system() returned.
/// \return An empty string when the command ran to completion with exit status
/// 0, otherwise a human readable reason for the failure.
[[nodiscard]] inline std::string commandFailure(int status) {
  if (status == -1) {
    return "could not be started";
  }
#ifndef _WIN32
  if (WIFSIGNALED(status)) {
    return "died on signal " + std::to_string(WTERMSIG(status));
  }
  if (!WIFEXITED(status)) {
    return "terminated abnormally";
  }
  if (WEXITSTATUS(status) != 0) {
    return "exited with status " + std::to_string(WEXITSTATUS(status));
  }
#else
  if (status != 0) {
    return "exited with status " + std::to_string(status);
  }
#endif
  return {};
}

/// \brief Run a shell command, throwing when it does not succeed.
///
/// External potentials feed whatever the command leaves behind straight into
/// the optimizer, so a failed command has to stop the run rather than let the
/// caller parse a stale or partial result file.
inline void runOrThrow(const std::string &command) {
  const std::string reason = commandFailure(std::system(command.c_str()));
  if (!reason.empty()) {
    throw std::runtime_error("Command \"" + command + "\" " + reason);
  }
}

} // namespace eonc::pot
