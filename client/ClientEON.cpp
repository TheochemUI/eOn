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
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#endif
#include "eon/EonLogger.h"
#ifdef _WIN32
#include <windows.h>
#endif

#include "eon/BaseStructures.h"
#include "eon/Bundling.h"
#include "eon/CommandLine.h"
#include "eon/EpiCenters.h"
#include "eon/HelperFunctions.h"
#include "eon/Job.h"
#include "eon/JobResult.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/Runtime.h"
#include "version.h"
#include <cstdlib>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <string_view>

#include <cerrno>
#include <chrono>
#include <cstring>
#include <ctime>
#include <filesystem>

#ifdef EONMPI
#include "eon/ParametersMpi.h"
#include <Python.h>
#include <cstdlib>
#include <fcntl.h>
#include <mpi.h>
#include <sstream>
#endif

// Includes for FPE trapping
#include "eon/fpe_handler.h"

#ifdef _WIN32
#include <float.h>
#endif

#ifndef _WIN32
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <unistd.h>
#endif

#ifdef __APPLE__
#ifndef __aarch64__
#include <mach/mach.h>
#include <mach/task_info.h>

void print_memory_usage() {
  auto *log = eonc::log::get();
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO,
                                (task_info_t)&t_info, &t_info_count)) {
    QUILL_LOG_ERROR(log, "Failed to get task info");
    return;
  }

  unsigned int rss = t_info.resident_size;
  unsigned int vs = t_info.virtual_size;
  QUILL_LOG_INFO(
      log,
      "\nmemory usage:\nresident size (MB): {:8.2f}\nvirtual size (MB):  "
      "{:8.2f}",
      static_cast<double>(rss) / 1024 / 1024,
      static_cast<double>(vs) / 1024 / 1024);
}
#endif
#endif

void printSystemInfo() {
  auto *log = eonc::log::get();
  QUILL_LOG_INFO(log, "eOn Client");
  QUILL_LOG_INFO(log, "{}", VERSION_STRING);
#ifndef __aarch64__
  QUILL_LOG_INFO(log, "OS: {}", OS_INFO);
  QUILL_LOG_INFO(log, "Arch: {}", ARCH);
#endif

#ifdef _WIN32
  TCHAR hostname[MAX_COMPUTERNAME_LENGTH + 1];
  DWORD size = sizeof(hostname) / sizeof(hostname[0]);
  if (GetComputerName(hostname, &size)) {
    QUILL_LOG_INFO(log, "Hostname: {}", hostname);
  } else {
    QUILL_LOG_ERROR(log, "Failed to get hostname");
  }
  QUILL_LOG_INFO(log, "PID: {}", GetCurrentProcessId());
#else
  struct utsname systemInfo;
  int status = uname(&systemInfo);
  if (status == 0) {
    QUILL_LOG_INFO(log, "Hostname: {}", systemInfo.nodename);
    QUILL_LOG_INFO(log, "PID: {}", getpid());
  } else {
    QUILL_LOG_ERROR(log, "Failed to get system information");
  }
#endif

  std::filesystem::path cwd = std::filesystem::current_path();
  QUILL_LOG_INFO(log, "DIR: {}", cwd.string());
}

static int eonClientMain(int argc, char **argv) {
  eonc::Parameters parameters;

  // Help, version, and one-shot argv jobs must not pay logger setup first.
  // commandLine starts the backend only after the flag parse commits to work.
  // In an MPI build only a rank an eOn server launched (EON_SERVER_PATH set,
  // EON_CLIENT_STANDALONE unset) skips the parse: the server passes it no
  // flags. A standalone or hand-started MPI client reads its flags like the
  // serial client; rgpot initialises MPI if a calculator group needs it.
#ifdef EONMPI
  const bool serverRank = getenv("EON_CLIENT_STANDALONE") == nullptr &&
                          getenv("EON_SERVER_PATH") != nullptr;
#else
  const bool serverRank = false;
#endif
  if (argc > 1 && !serverRank) {
    eonc::commandLine(argc, argv);
    quill::Backend::stop();
    return 0;
  }

  // A missing parameter file is not a job. Reject it before the log
  // thread and the file sinks. Help and version already returned above.
  if (!serverRank) {
    std::string configName = parameters.main_options().iniFilename;
    // A bundle is config_0.ini plus its siblings. That set is a parameter
    // file even when config.ini itself is absent.
    if (eonc::helpers::existsFile("config_0.ini")) {
      configName = "config_0.ini";
    } else {
      configName = eonc::helpers::getRelevantFile(configName);
    }
    if (!eonc::helpers::existsFile(configName)) {
      std::cerr << "Can't load INI file: " << configName << '\n'
                << "problem loading parameter file, stopping\n";
      return 1;
    }
  }

  // Quill backend, file sinks, and the combi logger. Deferred until a job
  // path actually runs so process start is not dominated by the log thread.
  auto *logger = eonc::log::init_client();
  if (!logger) {
    logger = eonc::log::get();
  }
  // File sinks open relative to this directory. MPI jobs chdir later.
  const auto logHome = std::filesystem::current_path();

#ifdef EONMPI
  // The same rule as the flag parse above: only a rank an eOn server
  // launched (EON_SERVER_PATH set, EON_CLIENT_STANDALONE unset) runs as a
  // server-driven client. Any other MPI client runs config.ini standalone.
  const bool client_standalone = !serverRank;
  int number_of_clients;
  if (!client_standalone) {
    if (getenv("EON_NUMBER_OF_CLIENTS") == nullptr) {
      QUILL_LOG_ERROR(logger,
                      "error: must set the env var EON_NUMBER_OF_CLIENTS");
      logger->flush_log();
      return 1;
    }
    number_of_clients = atoi(getenv("EON_NUMBER_OF_CLIENTS"));
  } else {
    number_of_clients = 1;
  }

  int eon_mpi_inited = 0;
  MPI_Initialized(&eon_mpi_inited);
  if (!eon_mpi_inited) {
    int eon_mpi_provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &eon_mpi_provided);
  }

  int error;
  std::string config_file = "config.ini";
  if (client_standalone) {
    if (eonc::helpers::existsFile("config_0.ini")) {
      config_file = "config_0.ini";
    }
    QUILL_LOG_INFO(logger, "Loading parameter file {}", config_file);
    error = parameters.load(config_file);
  } else {
    QUILL_LOG_INFO(logger, "Loading parameter file {}",
                   parameters.main_options().iniFilename);
    error = parameters.load(parameters.main_options().iniFilename);
  }
  if (error) {
    QUILL_LOG_ERROR(logger, "problem loading parameter file");
    logger->flush_log();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // All ranks must have loaded Parameters before the process-type Allgather.
  MPI_Barrier(MPI_COMM_WORLD);

  int irank;
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);
  int isize;
  MPI_Comm_size(MPI_COMM_WORLD, &isize);

  std::vector<int> process_types(isize);
  int process_type;

  process_type = 1;

  MPI_Allgather(&process_type, 1, MPI_INT, &process_types[0], 1, MPI_INT,
                MPI_COMM_WORLD);

  int i, servers = 0, clients = 0, potentials = 0;
  int server_rank = -1;
  int my_client_number = -1;
  std::vector<int> client_ranks;
  for (i = 0; i < isize; i++) {
    switch (process_types[i]) {
    case 0:
      servers++;
      break;
    case 1:
      if (i == irank) {
        my_client_number = clients;
      }
      clients++;
      client_ranks.push_back(i);
      break;
    case 2:
      potentials++;
      break;
    }
  }

  if (clients < number_of_clients) {
    QUILL_LOG_ERROR(logger,
                    "didn't launch as many mpi client ranks as specified in "
                    "EON_NUMBER_OF_CLIENTS");
    logger->flush_log();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  clients = number_of_clients;

  if (parameters.potential_options().potential == eonc::PotType::MPI) {
    if (potentials == 0 || potentials % clients != 0) {
      QUILL_LOG_ERROR(logger,
                      "the MPI potential needs a nonzero number of potential "
                      "ranks divisible by EON_NUMBER_OF_CLIENTS ({} potential "
                      "ranks, {} clients)",
                      potentials, clients);
      logger->flush_log();
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    std::vector<int> potential_ranks(potentials);
    int j;
    for (i = 0, j = 0; i < isize; i++) {
      if (process_types[i] == 2) {
        potential_ranks[j] = i;
        j++;
      }
    }
    int potential_group_size = potentials / clients;

    for (i = 0; i < clients; i++) {
      MPI_Group orig_group, new_group;
      MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
      int offset = i * potential_group_size;
      MPI_Group_incl(orig_group, potential_group_size, &potential_ranks[offset],
                     &new_group);
      MPI_Comm pot_comm;
      MPI_Comm_create(MPI_COMM_WORLD, new_group, &pot_comm);
    }

    if (my_client_number < number_of_clients) {
      parameters.set_mpi_potential_rank(
          potential_ranks[my_client_number * potential_group_size]);
    }
  }

  // LAMMPS MPI communicator setup (runtime check, not compile-time)
  if (parameters.potential_options().potential == eonc::PotType::LAMMPS) {
    for (i = 0; i < static_cast<int>(client_ranks.size()); i++) {
      MPI_Group world_group, new_group;
      MPI_Comm_group(MPI_COMM_WORLD, &world_group);
      int r = client_ranks[i];
      MPI_Group_incl(world_group, 1, &r, &new_group);
      MPI_Comm new_comm;
      MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
      if (new_comm != MPI_COMM_NULL) {
        eonc::setMpiClientComm(parameters, new_comm);
      }
      QUILL_LOG_INFO(logger, "creating group with ranks: {}", r);
    }
  }

  if (!client_standalone) {
    server_rank = client_ranks[number_of_clients];
    if (my_client_number == number_of_clients) {
      std::ostringstream oss;
      oss << client_ranks.at(0);
      for (i = 1; i < number_of_clients; i++) {
        oss << ":" << client_ranks.at(i);
      }
      setenv("EON_CLIENT_RANKS", oss.str().c_str(), 1);

      wchar_t **py_argv =
          static_cast<wchar_t **>(malloc(sizeof(wchar_t *) * 2));
      py_argv[0] = Py_DecodeLocale(argv[0], nullptr);
      char *program = getenv("EON_SERVER_PATH");
      py_argv[1] = Py_DecodeLocale(program, nullptr);
      QUILL_LOG_INFO(logger, "rank: {} becoming {}", irank, program);
      Py_Initialize();
      Py_Main(2, py_argv);
      Py_FinalizeEx();
      // GH
      MPI_Finalize();
      return 0;
    } else if (my_client_number > number_of_clients) {
      MPI_Finalize();
      return 0;
    }
  }
#endif

  eonc::enableFPE(); // from ExceptionsEON.h

#ifdef EONMPI
  // Server sends a path starting with STOPCAR to end this loop.
  char logfilename[1024];
  snprintf(logfilename, 1024, "eonclient_%i.log", my_client_number);

  auto orig_path = std::filesystem::current_path();
  // In client/server mode a failed job is logged and its directory handed
  // back without results.dat; the server skips it and this rank stays in
  // the pool. A standalone rank still exits with the error.
  const bool keepServing = !client_standalone;
  while (true) {
    std::filesystem::current_path(orig_path);
    std::string path(1024, '\0');
    int ready = 1;
    if (!client_standalone) {
      QUILL_LOG_INFO(
          logger, "client: rank {} is ready, posting send to server rank: {}!",
          irank, server_rank);
      // Tag "1" is to interrupt the main loop and tell the communicator that a
      // client is ready
      {
        MPI_Request eon_rq;
        MPI_Isend(&ready, 1, MPI_INT, server_rank, 1, MPI_COMM_WORLD, &eon_rq);
        MPI_Request_free(&eon_rq);
      }

      // Get the path we should run in from the server
      MPI_Recv(&path[0], 1024, MPI_CHAR, server_rank, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      if (path.starts_with("STOPCAR")) {
        QUILL_LOG_INFO(logger, "rank {} got STOPCAR", irank);
        MPI_Finalize();
        return 0;
      }
      QUILL_LOG_INFO(logger, "client: rank: {} chdir to {}", irank, path);

      if (const auto chdirError =
              eonc::helpers::enterJobDirectory(path.c_str())) {
        QUILL_LOG_ERROR(logger, "error: chdir: {}", *chdirError);
        logger->flush_log();
        MPI_Send(&path[0], 1024, MPI_CHAR, server_rank, 0, MPI_COMM_WORLD);
        continue;
      }
    }
#else
  constexpr bool keepServing = false;
#endif
    // Flushes the logs into the job directory so a failed job still shows
    // its error to whoever reads the returned directory.
    auto stageFailedJobLogs = [&] {
      logger->flush_log();
      if (auto *trace =
              quill::Frontend::get_logger(std::string{"_traceback"})) {
        trace->flush_log();
      }
      for (const std::string_view logName :
           {std::string_view{"client_quill.log"},
            std::string_view{"client_traceback.log"}}) {
        eonc::helpers::stageReturnLog(logHome.string(), logName);
      }
    };

    printSystemInfo();

    bool bundlingEnabled = false;
    int bundleSize = eonc::getBundleSize();
    if (bundleSize <= 0) {
      bundleSize = 1;
      bundlingEnabled = false;
    } else {
      bundlingEnabled = true;
    }

    std::vector<std::string> bundledFilenames;
    for (int i = 0; i < bundleSize; i++) {
      // This job only. A clock above the bundle or MPI loop also counts
      // earlier jobs and the idle wait between them.
      const auto start_time = std::chrono::steady_clock::now();
      if (bundleSize > 1)
        QUILL_LOG_INFO(logger, "Beginning Job {} of {}", i + 1, bundleSize);
      std::vector<std::string> unbundledFilenames;
      if (bundlingEnabled) {
        unbundledFilenames = eonc::unbundle(i);
      }

      // check to see if parameters file exists before loading
      int error = 0;
      std::string config_file =
          eonc::helpers::getRelevantFile(parameters.main_options().iniFilename);
      QUILL_LOG_INFO(logger, "Loading parameter file {}", config_file);
      error = parameters.load(config_file);

      if (error) {
        QUILL_LOG_ERROR(logger, "problem loading parameter file, stopping");
        logger->flush_log();
        if (keepServing) {
          stageFailedJobLogs();
          continue;
        }
        return 1;
      }

      // Determine what type of job we are running according to the parameters
      // file.
      eonc::Runtime rt;
      auto job = eonc::helpers::makeJob(
          std::make_unique<eonc::Parameters>(parameters), rt);
      if (job == nullptr) {
        QUILL_LOG_ERROR(logger, "error: Unknown job: {}",
                        std::string{magic_enum::enum_name<eonc::JobType>(
                            parameters.main_options().job)});
        logger->flush_log();
        if (keepServing) {
          stageFailedJobLogs();
          continue;
        }
        return 1;
      }

      std::vector<std::string> filenames;
      try {
        filenames = job->run();
      } catch (int e) {
        QUILL_LOG_CRITICAL(logger, "[ERROR] job exited on error {}", e);
        logger->flush_log();
        if (keepServing) {
          stageFailedJobLogs();
          continue;
        }
        return EXIT_FAILURE;
      } catch (const std::exception &e) {
        QUILL_LOG_CRITICAL(logger, "[ERROR] unhandled exception: {}", e.what());
        logger->flush_log();
        std::cerr << "[ERROR] unhandled exception: " << e.what() << "\n";
        if (keepServing) {
          stageFailedJobLogs();
          continue;
        }
        return EXIT_FAILURE;
      }

      job->releasePotential();
      rt.pots().write_summary();
      job.reset();
      filenames.push_back(std::string("_potcalls.json"));

      // Finalize Timing Information
      auto end_time = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed = end_time - start_time;

      double utime = 0, stime = 0, rtime = 0;
      eonc::helpers::getTime(&rtime, &utime, &stime);

      QUILL_LOG_INFO(logger, "Timing Information:");
      QUILL_LOG_INFO(logger, "  Real time: {:.3f} seconds", elapsed.count());
      QUILL_LOG_INFO(logger, "  User time: {:.3f} seconds", utime);
      QUILL_LOG_INFO(logger, "  System time: {:.3f} seconds", stime);

      // results.dat contract is "<value> <key>" (same as all job writers and
      // eon.fileio.parse_results / eon_schema.jobs adapters).
      std::ofstream result_file("results.dat", std::ios::app);
      if (result_file.is_open()) {
        result_file << std::format("{:.12e} time_seconds\n", elapsed.count());
#ifndef _WIN32
        result_file << std::format("{:.12e} user_time\n", utime);
        result_file << std::format("{:.12e} system_time\n", stime);
#endif
        const eonc::JobResultProvenance provenance =
            eonc::provenanceForJob(parameters);
        bool have_backend = false;
        {
          std::ifstream prior("results.dat");
          std::string line;
          while (std::getline(prior, line)) {
            if (line.find("optimizer_backend") != std::string::npos) {
              have_backend = true;
              break;
            }
          }
        }
        if (!have_backend) {
          result_file << provenance.text();
        }
      } else {
        QUILL_LOG_ERROR(logger, "Failed to write timing to results.dat");
      }
      result_file.close();

      logger->flush_log();
      if (auto *trace =
              quill::Frontend::get_logger(std::string{"_traceback"})) {
        trace->flush_log();
      }
      // Quill opens these before an MPI job chdir. Copy them into the job
      // directory so return_files.dat does not name a missing log.
      constexpr std::string_view jobLogs[] = {"client_quill.log",
                                              "client_traceback.log"};
      for (const std::string_view logName : jobLogs) {
        if (eonc::helpers::stageReturnLog(logHome.string(), logName)) {
          filenames.push_back(std::string{logName});
        }
      }

      {
        std::ofstream manifest("return_files.dat");
        if (manifest) {
          for (const auto &fn : filenames) {
            manifest << fn << "\n";
          }
        }
        filenames.push_back(std::string("return_files.dat"));
      }

      if (bundlingEnabled) {
        eonc::bundle(i, filenames, &bundledFilenames);
        eonc::deleteUnbundledFiles(unbundledFilenames);
      } else {
        bundledFilenames = filenames;
      }
    }

#ifdef EONMPI
    if (client_standalone) {
      break;
    }
    {
      MPI_Request eon_rq;
      MPI_Isend(&path[0], 1024, MPI_CHAR, server_rank, 0, MPI_COMM_WORLD,
                &eon_rq);
      // path is destroyed at the end of this iteration. Request_free does
      // not complete the send, so the buffer stays live until Wait returns.
      MPI_Wait(&eon_rq, MPI_STATUS_IGNORE);
    }

    // End of MPI while loop
  }
#endif

#ifdef OSX
#ifndef __aarch64__
  print_memory_usage();
#endif
#endif

#ifdef EONMPI
  // Clean shutdown for both the standalone single-rank case and the
  // client/server case; the standalone path previously aborted MPI_COMM_WORLD.
  MPI_Finalize();
#endif

  // Ensure all queued log messages are flushed before exiting
  quill::Backend::stop();
  return EXIT_SUCCESS;
}

namespace {

int reportFatal(std::string_view what) {
  std::cerr << "eonclient: fatal error: " << what << '\n';
  // Drain the queued CRITICAL lines that say why: an escaping exception
  // reaches std::terminate without them, and so does std::exit.
  quill::Backend::stop();
#ifdef EONMPI
  int mpiReady = 0;
  MPI_Initialized(&mpiReady);
  if (mpiReady) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
#endif
  return EXIT_FAILURE;
}

} // namespace

int main(int argc, char **argv) {
  try {
    return eonClientMain(argc, argv);
  } catch (const std::exception &e) {
    return reportFatal(e.what());
  } catch (...) {
    return reportFatal("exception of unknown type");
  }
}
