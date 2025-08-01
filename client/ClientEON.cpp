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

#include "CommandLine.h"
#include "HelperFunctions.h"
#include "Job.hpp"
#include "Log.hpp"
#include "Parameters.h"
#include "Parser.hpp"
#include "client/StructureComparisonJob.h"
#include "client/matter/MatterHelpers.hpp"
#include "potentials/PotentialCache.hpp"
#include "version.h"

#include <cstdlib>
#include <filesystem>
#include <string.h>
#include <time.h>

#if defined WITH_ASE_ORCA || EMBED_PYTHON || WITH_ASE_NWCHEM
#include <pybind11/embed.h>
#endif

// Includes for FPE trapping
#include "fpe_handler.h"

#ifdef WIN32
#include <float.h>
#endif

#ifndef WIN32
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <unistd.h>
#endif

#ifdef __APPLE__
#ifndef __aarch64__
#include <mach/mach.h>
#include <mach/task_info.h>

void print_memory_usage(std::shared_ptr<spdlog::logger> log) {
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO,
                                (task_info_t)&t_info, &t_info_count)) {
    SPDLOG_LOGGER_ERROR(log, "Failed to get task info");
    return;
  }

  unsigned int rss = t_info.resident_size;
  unsigned int vs = t_info.virtual_size;
  // SPDLOG_LOGGER_INFO("\nmemory usage:\nresident size (MB): {}\nvirtual "
  //                    "size (MB):  {}",
  //                    static_cast<double>(rss) / 1024 / 1024,
  //                    static_cast<double>(vs) / 1024 / 1024);
}
#endif
#endif

void printSystemInfo(std::shared_ptr<spdlog::logger> log) {
  SPDLOG_LOGGER_INFO(log, "EON Client");
  SPDLOG_LOGGER_INFO(log, "VERSION: {}", VERSION);
  SPDLOG_LOGGER_INFO(log, "BUILD DATE: {}\n", BUILD_DATE);
#ifndef __aarch64__
  SPDLOG_LOGGER_INFO(log, "OS: {}", OS_INFO);
  SPDLOG_LOGGER_INFO(log, "Arch: {}", ARCH);
#endif

#ifdef _WIN32
  TCHAR hostname[MAX_COMPUTERNAME_LENGTH + 1];
  DWORD size = sizeof(hostname) / sizeof(hostname[0]);
  if (GetComputerName(hostname, &size)) {
    SPDLOG_LOGGER_INFO(log, "Hostname: {}", hostname);
  } else {
    SPDLOG_LOGGER_ERROR(log, "Failed to get hostname");
  }
  SPDLOG_LOGGER_INFO(log, "PID: {}", GetCurrentProcessId());
#else
  struct utsname systemInfo;
  int status = uname(&systemInfo);
  if (status == 0) {
    SPDLOG_LOGGER_INFO(log, "Hostname: {}", systemInfo.nodename);
    SPDLOG_LOGGER_INFO(log, "PID: {}", getpid());
  } else {
    SPDLOG_LOGGER_ERROR(log, "Failed to get system information");
  }
#endif

  std::filesystem::path cwd = std::filesystem::current_path();
  SPDLOG_LOGGER_INFO(log, "DIR: {}", cwd.string());
}

int main(int argc, char **argv) {
  // --- Start Logging setup
  eonc::LogManager::getInstance();
  auto logger = spdlog::get("combi");
  // End logging setup
  eonc::Parameters parameters;

#if defined WITH_ASE_ORCA || EMBED_PYTHON || WITH_ASE_NWCHEM
  pybind11::scoped_interpreter guard{};
#endif

  auto params = eonc::commandLine(logger, argc, argv);

  eonc::enableFPE(); // from ExceptionsEON.h

  double beginTime = 0.0;

  if (not params["Main"]["quiet"].as<bool>()) {
    printSystemInfo(logger);
    eonc::helper_functions::getTime(&beginTime, NULL, NULL);
  }
  auto pcache = eonc::cache::PotentialCache::getInstance();
  auto pot = eonc::makePotential(params);
  pot->set_cache(pcache);

  // std::cout << params << std::endl;
  auto mats = eonc::mat::make_matter(params, pot);
  auto job = eonc::mkJob(params);

  bool result = eonc::runJob(job, mats);

  if (!result &&
      not std::holds_alternative<eonc::StructureComparisonJob>(job)) {
    SPDLOG_LOGGER_ERROR(logger, "Job execution failed");
    // throw std::runtime_error("Something went wrong running the job");
  }

  if (not params["Main"]["quiet"].as<bool>()) {
    // Timing Information
    double utime = 0, stime = 0, rtime = 0;
    eonc::helper_functions::getTime(&rtime, &utime, &stime);
    rtime = rtime - beginTime;

    // if (Potential::totalUserTime > 0) {
    //     printf("\ntime not in potential: %.4f%%\n",
    //     100*(1-Potential::totalUserTime/utime));
    // }

    printf(
        "timing information:\nreal %10.3f seconds\nuser %10.3f seconds\nsys  "
        "%10.3f seconds\n",
        rtime, utime, stime);
  }

  std::exit(EXIT_SUCCESS);
}
