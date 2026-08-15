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

#include <cctype>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <system_error>

namespace eonc::pot {

/// \brief Wrap a path for the shell system() hands the command to, so spaces
/// in it survive.
inline std::string shellQuote(const std::string &text) {
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

/// \brief Return the interpreter name from a `#!` line, or empty.
///
/// `#!/usr/bin/env python3` and `#!/usr/bin/python3` both yield `python3`.
/// `env` flags (`-S`) are skipped. A file with no shebang yields empty.
inline std::string shebangInterpreter(const std::filesystem::path &path) {
  std::ifstream in(path);
  std::string line;
  if (!in || !std::getline(in, line)) {
    return {};
  }
  if (!line.empty() && line.back() == '\r') {
    line.pop_back();
  }
  if (line.size() < 2 || line[0] != '#' || line[1] != '!') {
    return {};
  }
  std::string rest = line.substr(2);
  auto take = [](std::string &s) -> std::string {
    const auto begin = s.find_first_not_of(" \t");
    if (begin == std::string::npos) {
      s.clear();
      return {};
    }
    const auto end = s.find_first_of(" \t", begin);
    std::string tok = end == std::string::npos ? s.substr(begin)
                                               : s.substr(begin, end - begin);
    if (end == std::string::npos) {
      s.clear();
    } else {
      s.erase(0, end);
    }
    return tok;
  };
  const std::string first = take(rest);
  if (first.empty()) {
    return {};
  }
  std::string name = std::filesystem::path(first).filename().string();
  if (name == "env") {
    std::string tok = take(rest);
    while (!tok.empty() && tok.front() == '-') {
      tok = take(rest);
    }
    if (tok.empty()) {
      return {};
    }
    name = std::filesystem::path(tok).filename().string();
  }
  return name;
}

/// \brief True when \p name is python, python3, python3.12, python.exe, ...
inline bool isPythonInterpreterName(std::string name) {
  for (char &c : name) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  constexpr std::string_view exe{".exe"};
  if (name.size() > exe.size() &&
      name.compare(name.size() - exe.size(), exe.size(), exe.data()) == 0) {
    name.resize(name.size() - exe.size());
  }
  return name.rfind("python", 0) == 0;
}

/// \brief Quote a resolved program path, and on Windows prefix a script so
/// cmd.exe can start it.
///
/// cmd.exe will not execute an extensionless file, quoted or not. A Python
/// shebang or a `.py` suffix becomes `python "path"`. A native `.exe` /
/// `.bat` / `.cmd` is quoted only. A command line that is not a single
/// existing file is returned unchanged.
inline std::string wrapResolvedProgram(const std::string &command,
                                       const std::string &python) {
  std::error_code ec;
  if (!std::filesystem::is_regular_file(command, ec) || ec) {
    return command;
  }
#ifdef _WIN32
  const std::filesystem::path named{command};
  std::string ext = named.extension().string();
  for (char &c : ext) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  if (ext == ".exe" || ext == ".com" || ext == ".bat" || ext == ".cmd") {
    return shellQuote(command);
  }
  const std::string interp = shebangInterpreter(named);
  if (isPythonInterpreterName(interp) || ext == ".py" || ext == ".pyw") {
    const std::string runner = python.empty() ? std::string{"python"} : python;
    return runner + " " + shellQuote(command);
  }
  if (!interp.empty()) {
    return interp + " " + shellQuote(command);
  }
  return shellQuote(command);
#else
  (void)python;
  return shellQuote(command);
#endif
}

} // namespace eonc::pot
