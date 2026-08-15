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

#include "eon/potentials/ExtPot/ExtPot.h"
#include "eon/potentials/ExternalCommand.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <filesystem>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

namespace {

// The names the external program reads and writes. They are part of the
// documented interface and stay fixed; the directory holding them is what
// varies from one client to the next.
constexpr const char *kToExtPot = "from_eon_to_extpot";
constexpr const char *kFromExtPot = "from_extpot_to_eon";

// Names the directory eOn was started in for a wrapper that needs to reach
// back to it.
constexpr const char *kRunDirVariable = "EON_EXTPOT_RUN_DIR";

long processId() {
#ifdef _WIN32
  return static_cast<long>(_getpid());
#else
  return static_cast<long>(getpid());
#endif
}

/// \brief Wrap a path for the shell system() hands the command to, so spaces
/// in it survive.
std::string shellQuote(const std::string &text) {
#ifdef _WIN32
  // cmd.exe takes double quotes and offers no escape for an embedded one,
  // which a Windows path cannot hold anyway.
  return "\"" + text + "\"";
#else
  std::string quoted{"'"};
  for (const char c : text) {
    if (c == '\'') {
      quoted += "'\\''";
    } else {
      quoted += c;
    }
  }
  quoted += "'";
  return quoted;
#endif
}

#ifdef _WIN32
bool isWindowsNativeProgram(const std::filesystem::path &path) {
  std::string ext = path.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return ext == ".exe" || ext == ".bat" || ext == ".cmd" || ext == ".com";
}

/// \brief cmd.exe will not honor a shebang. The default ext_pot wrapper is a
/// Python script with no suffix, so launch it with python.
std::string windowsLaunch(const std::string &command,
                          const std::string &quotedCommand) {
  std::error_code ec;
  const std::filesystem::path named{command};
  if (!std::filesystem::is_regular_file(named, ec) || ec ||
      isWindowsNativeProgram(named)) {
    return quotedCommand;
  }
  return std::format("python {}", quotedCommand);
}
#endif

/// \brief Build the command line that runs \p command in \p workDir with
/// \p runDir named in the environment.
///
/// Changing directory in the shell rather than in this process keeps the
/// exchange files private to one client without moving the working directory
/// out from under the threads evaluating other images.
std::string composeCommand(const std::string &workDir,
                           const std::string &runDir,
                           const std::string &command) {
  // A resolved single-file path (the default ./ext_pot after
  // resolveCommand) has to be quoted so spaces and metacharacters
  // survive system(). A command line with arguments is left alone.
  std::error_code quoteEc;
  const std::string quotedCommand =
      std::filesystem::is_regular_file(command, quoteEc) && !quoteEc
          ? shellQuote(command)
          : command;
#ifdef _WIN32
  // cmd.exe: quotes around name=value so the value does not include them.
  return std::format("cd /d {} && set \"{}={}\" && {}", shellQuote(workDir),
                     kRunDirVariable, runDir,
                     windowsLaunch(command, quotedCommand));
#else
  return std::format("cd {} && export {}={} && {}", shellQuote(workDir),
                     kRunDirVariable, shellQuote(runDir), quotedCommand);
#endif
}

/// \brief Resolve a command that names a file in the current directory.
///
/// The external command runs from the exchange directory, so a path written
/// relative to the directory eOn runs in, such as the default ./ext_pot, has
/// to be resolved while that directory is still current. A command that does
/// not name an existing file is left alone: it is a name for the shell to
/// look up on PATH, or a command line with arguments.
std::string resolveCommand(const std::string &command) {
  std::error_code ec;
  const std::filesystem::path named{command};
  if (!std::filesystem::is_regular_file(named, ec) || ec) {
    return command;
  }
  const std::filesystem::path resolved = std::filesystem::absolute(named, ec);
  if (ec) {
    return command;
  }
  return resolved.string();
}

} // namespace

ExtPot::ExtPot(const Parameters &p)
    : Potential(p),
      eon_extpot_path{resolveCommand(p.potential_options.extPotPath)} {}

void ExtPot::cleanMemory(void) {
  if (exchangeDir.empty() || retainExchangeDir) {
    return;
  }
  std::error_code ec;
  std::filesystem::remove(exchangeDir / kToExtPot, ec);
  std::filesystem::remove(exchangeDir / kFromExtPot, ec);
  // Removes the directory only when it is empty, so whatever else the
  // external program wrote there stays.
  std::filesystem::remove(exchangeDir, ec);
}

ExtPot::~ExtPot() { cleanMemory(); }

void ExtPot::prepareExchangeDir() {
  std::error_code ec;
  if (exchangeDir.empty()) {
    // One directory per potential object. A client started alongside this
    // one carries a different process id, and a second potential in this
    // process a different index, so no two of them name the same file.
    static std::atomic<long> instanceCount{0};
    const std::filesystem::path named{
        std::format("extpot_{}_{}", processId(), instanceCount.fetch_add(1))};
    exchangeDir = std::filesystem::absolute(named, ec);
    if (ec) {
      exchangeDir.clear();
      throw std::runtime_error(std::format("Could not resolve {}: {}",
                                           named.string(), ec.message()));
    }
  }
  std::filesystem::create_directory(exchangeDir, ec);
  if (ec) {
    throw std::runtime_error(std::format("Could not create {} to run {} in: {}",
                                         exchangeDir.string(), eon_extpot_path,
                                         ec.message()));
  }
  // A result left by an earlier call, or by a client that died holding this
  // process id, is read as this call's forces when the external program
  // writes nothing.
  std::filesystem::remove(exchangeDir / kFromExtPot, ec);
}

void ExtPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  if (eon_extpot_path.empty()) {
    throw std::runtime_error(
        "ExtPot needs potential_options.extPotPath to name a program to run");
  }
  prepareExchangeDir();
  retainExchangeDir = true;
  passToSystem(N, R, atomicNrs, box);
  std::error_code ec;
  const std::filesystem::path runDir = std::filesystem::current_path(ec);
  eonc::pot::runOrThrow(composeCommand(exchangeDir.string(),
                                       ec ? std::string{} : runDir.string(),
                                       eon_extpot_path));
  recieveFromSystem(N, F, U);
  retainExchangeDir = false;
  return;
}

void ExtPot::passToSystem(long N, const double *R, const int *atomicNrs,
                          const double *box)
// 'positions' of all particles and box
{
  const std::filesystem::path toExtPot = exchangeDir / kToExtPot;
  std::ofstream out(toExtPot, std::ios::trunc);
  if (!out) {
    throw std::runtime_error(
        std::format("Could not open {} for writing", toExtPot.string()));
  }

  for (int i = 0; i < 3; i++) {
    out << std::format("{:.19f}\t{:.19f}\t{:.19f}\n", box[i * 3 + 0],
                       box[i * 3 + 1], box[i * 3 + 2]);
  }

  for (long i = 0; i < N; i++) {
    out << std::format("{}\t{:.19f}\t{:.19f}\t{:.19f}\n", atomicNrs[i],
                       R[i * 3 + 0], R[i * 3 + 1], R[i * 3 + 2]);
  }

  out.close();
  if (!out) {
    throw std::runtime_error(
        std::format("Could not write the structure to {}", toExtPot.string()));
  }
  return;
}

void ExtPot::recieveFromSystem(long N, double *F, double *U)
// first line must be the total 'energy', the following lines should be the
// 'forces'
{
  const std::filesystem::path fromExtPot = exchangeDir / kFromExtPot;
  std::ifstream in(fromExtPot);
  if (!in) {
    throw std::runtime_error(
        std::format("Could not open {}; the external program left no result",
                    fromExtPot.string()));
  }

  if (!(in >> *U)) {
    throw std::runtime_error(
        std::format("Could not read the energy from {}", fromExtPot.string()));
  }

  for (long i = 0; i < N; i++) {
    if (!(in >> F[i * 3 + 0] >> F[i * 3 + 1] >> F[i * 3 + 2])) {
      throw std::runtime_error(
          std::format("{} holds forces for {} atoms, expected {}",
                      fromExtPot.string(), i, N));
    }
  }
  return;
}
