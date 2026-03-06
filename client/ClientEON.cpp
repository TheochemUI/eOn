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
#include "EonLogger.h"
#include <windows.h>
#endif

#include "BaseStructures.h"
#include "Bundling.h"
#include "CommandLine.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "Job.h"
#include "Parameters.h"
#include "Potential.h"
#include "version.h"

#include <chrono>
#include <errno.h>
#include <filesystem>
#include <string.h>
#include <time.h>

#ifdef EONMPI
#include <Python.h>
#include <fcntl.h>
#include <mpi.h>
#include <sstream>
#include <stdlib.h>
#endif

#if defined WITH_ASE_ORCA || EMBED_PYTHON || WITH_ASE_NWCHEM
#include "PyGuard.h"
#endif

#ifdef EONMPIBGP
#include <libgen.h>
#endif

// Includes for FPE trapping
#include "fpe_handler.h"

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
  auto *log = quill::Frontend::get_logger("combi");
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO,
                                (task_info_t)&t_info, &t_info_count)) {
    LOG_ERROR(log, "Failed to get task info");
    return;
  }

  unsigned int rss = t_info.resident_size;
  unsigned int vs = t_info.virtual_size;
  LOG_INFO(log,
           "\nmemory usage:\nresident size (MB): {:8.2f}\nvirtual size (MB):  "
           "{:8.2f}",
           (double)rss / 1024 / 1024, (double)vs / 1024 / 1024);
}
#endif
#endif

void printSystemInfo() {
  auto *log = quill::Frontend::get_logger("combi");
  LOG_INFO(log, "eOn Client");
  LOG_INFO(log, "{}", VERSION_STRING);
#ifndef __aarch64__
  LOG_INFO(log, "OS: {}", OS_INFO);
  LOG_INFO(log, "Arch: {}", ARCH);
#endif

#ifdef _WIN32
  TCHAR hostname[MAX_COMPUTERNAME_LENGTH + 1];
  DWORD size = sizeof(hostname) / sizeof(hostname[0]);
  if (GetComputerName(hostname, &size)) {
    LOG_INFO(log, "Hostname: {}", hostname);
  } else {
    LOG_ERROR(log, "Failed to get hostname");
  }
  LOG_INFO(log, "PID: {}", GetCurrentProcessId());
#else
  struct utsname systemInfo;
  int status = uname(&systemInfo);
  if (status == 0) {
    LOG_INFO(log, "Hostname: {}", systemInfo.nodename);
    LOG_INFO(log, "PID: {}", getpid());
  } else {
    LOG_ERROR(log, "Failed to get system information");
  }
#endif

  std::filesystem::path cwd = std::filesystem::current_path();
  LOG_INFO(log, "DIR: {}", cwd.string());
}

int main(int argc, char **argv) {
  // --- Start Logging setup
  // Configure backend for optimal performance (see BackendOptions.h)
  quill::BackendOptions backend_options;
  // Use 10μs sleep for balanced performance (10× faster than 100μs default)
  // Set to 0 for maximum throughput at cost of 100% CPU on backend thread
  backend_options.sleep_duration = std::chrono::microseconds{10};
  // Reduce timestamp ordering grace period for lower latency
  backend_options.log_timestamp_ordering_grace_period =
      std::chrono::microseconds{1};
  // Flush more frequently for better responsiveness
  backend_options.sink_min_flush_interval = std::chrono::milliseconds{100};
  quill::Backend::start(backend_options);
  auto console_sink =
      quill::Frontend::create_or_get_sink<quill::ConsoleSink>("console");
  auto file_sink = quill::Frontend::create_or_get_sink<quill::FileSink>(
      "client_quill.log",
      []() {
        quill::FileSinkConfig cfg;
        cfg.set_open_mode('w');
        return cfg;
      }(),
      quill::FileEventNotifier{});
  auto *logger = quill::Frontend::create_or_get_logger(
      "combi", {std::move(console_sink), std::move(file_sink)},
      quill::PatternFormatterOptions{quill::PatternFormatterOptions{
          quill::PatternFormatterOptions{"%(message)"}}},
      quill::ClockSourceType::System);
  logger->set_log_level(quill::LogLevel::TraceL3);
  // Traceback logger
  auto trace_csink =
      quill::Frontend::create_or_get_sink<quill::ConsoleSink>("trace_console");
  auto trace_fsink = quill::Frontend::create_or_get_sink<quill::FileSink>(
      "client_traceback.log",
      []() {
        quill::FileSinkConfig cfg;
        cfg.set_open_mode('w');
        return cfg;
      }(),
      quill::FileEventNotifier{});
  quill::Frontend::create_or_get_logger(
      "_traceback", {std::move(trace_csink), std::move(trace_fsink)},
      quill::PatternFormatterOptions{quill::PatternFormatterOptions{
          " [%(log_level)] [%(source_location)] [%(caller_function)] \n "
          "%(message)\n[end %(log_level)]"}},
      quill::ClockSourceType::System);
  //--- End logging setup
  Parameters parameters;

#if defined WITH_ASE_ORCA || EMBED_PYTHON || WITH_ASE_NWCHEM
  eonc::ensure_interpreter();
#endif

#ifdef EONMPI
  bool client_standalone = false;
  if (getenv("EON_CLIENT_STANDALONE") != NULL) {
    client_standalone = true;
  }
  int number_of_clients;
  if (!client_standalone) {
    if (getenv("EON_SERVER_PATH") == NULL) {
      LOG_ERROR(logger, "error: must set the env var EON_SERVER_PATH");
      logger->flush_log();
      return 1;
    }
    if (getenv("EON_NUMBER_OF_CLIENTS") == NULL) {
      LOG_ERROR(logger, "error: must set the env var EON_NUMBER_OF_CLIENTS");
      logger->flush_log();
      return 1;
    }
    number_of_clients = atoi(getenv("EON_NUMBER_OF_CLIENTS"));
  } else {
    number_of_clients = 1;
  }

  if (MPI::Is_initialized() == false) {
    MPI::Init();
  }

  int error;
  string config_file = "config.ini";
  if (client_standalone) {
    if (helper_functions::existsFile("config_0.ini")) {
      config_file = "config_0.ini";
    }
    LOG_INFO(logger, "Loading parameter file {}", config_file);
    error = parameters.load(config_file);
  } else {
    LOG_INFO(logger, "Loading parameter file {}",
             parameters.main_options.iniFilename);
    error = parameters.load(parameters.main_options.iniFilename);
  }
  if (error) {
    LOG_ERROR(logger, "problem loading parameter file");
    logger->flush_log();
    MPI::COMM_WORLD.Abort(1);
  }

  // XXX: Barrier for gpaw-python
  MPI::COMM_WORLD.Barrier();

  int irank = MPI::COMM_WORLD.Get_rank();
  int isize = MPI::COMM_WORLD.Get_size();

  int *process_types = new int[isize];
  int process_type;

  process_type = 1;

  MPI::COMM_WORLD.Allgather(&process_type, 1, MPI::INT, &process_types[0], 1,
                            MPI::INT);

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
    LOG_ERROR(logger, "didn't launch as many mpi client ranks as specified in "
                      "EON_NUMBER_OF_CLIENTS");
    logger->flush_log();
    MPI::COMM_WORLD.Abort(1);
  }
  clients = number_of_clients;

  if (parameters.potential_options.potential == "mpi") {
    int *potential_ranks = new int[potentials];
    int j;
    for (i = 0, j = 0; i < isize; i++) {
      if (process_types[i] == 2) {
        potential_ranks[j] = i;
        j++;
      }
    }
    int potential_group_size = potentials / clients;

    for (i = 0; i < clients; i++) {
      MPI::Group orig_group, new_group;
      orig_group = MPI::COMM_WORLD.Get_group();
      int offset = i * potential_group_size;
      new_group =
          orig_group.Incl(potential_group_size, &potential_ranks[offset]);
      MPI::COMM_WORLD.Create(new_group);
    }

    if (my_client_number < number_of_clients) {
      parameters.potential_options.MPIPotentialRank =
          potential_ranks[my_client_number * potential_group_size];
    }
  }

#ifdef LAMMPS_POT
  for (i = 0; i < int(client_ranks.size()); i++) {
    MPI_Group world_group, new_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    int r = client_ranks[i];
    MPI_Group_incl(world_group, 1, &r, &new_group);
    MPI_Comm new_comm;
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
    if (new_comm != MPI_COMM_NULL) {
      parameters.potential_options.MPIClientComm = new_comm;
    }
    LOG_INFO(logger, "creating group with ranks: {}", r);
  }
#endif

  if (!client_standalone) {
    server_rank = client_ranks[number_of_clients];
    if (my_client_number == number_of_clients) {
      std::ostringstream oss;
      oss << client_ranks.at(0);
      for (i = 1; i < number_of_clients; i++) {
        oss << ":" << client_ranks.at(i);
      }
      setenv("EON_CLIENT_RANKS", oss.str().c_str(), 1);

      wchar_t **py_argv = (wchar_t **)malloc(sizeof(wchar_t *) * 2);
      py_argv[0] = Py_DecodeLocale(argv[0], NULL);
      char *program = getenv("EON_SERVER_PATH");
      py_argv[1] = Py_DecodeLocale(program, NULL);
      LOG_INFO(logger, "rank: {} becoming {}", irank, program);
      Py_Initialize();
      Py_Main(2, py_argv);
      Py_FinalizeEx();
      // GH
      MPI::Finalize();
      return 0;
    } else if (my_client_number > number_of_clients) {
      MPI::Finalize();
      return 0;
    }
  }
#endif

#ifndef EONMPI
  if (argc > 1) {
    commandLine(argc, argv);
    return 0;
  }
#endif

  eonc::enableFPE(); // from ExceptionsEON.h

  auto start_time = std::chrono::steady_clock::now();

#ifdef EONMPI
  // XXX: When do we stop? The server should probably tell everyone when to
  // stop.
  char logfilename[1024];
  snprintf(logfilename, 1024, "eonclient_%i.log", my_client_number);

  char *orig_path = new char[1024];
  getcwd(orig_path, 1024);
  while (true) {
    chdir(orig_path);
    char *path = new char[1024];
    int ready = 1;
    if (!client_standalone) {
      LOG_INFO(logger,
               "client: rank {} is ready, posting send to server rank: {}!",
               irank, server_rank);
      // Tag "1" is to interrupt the main loop and tell the communicator that a
      // client is ready
      MPI::COMM_WORLD.Isend(&ready, 1, MPI::INT, server_rank, 1);

      // Get the path we should run in from the server
      MPI::COMM_WORLD.Recv(&path[0], 1024, MPI::CHAR, server_rank, 0);
      if (strncmp("STOPCAR", path, 1024) == 0) {
        LOG_INFO(logger, "rank {} got STOPCAR", irank);
        MPI::Finalize();
        return 0;
      }
      LOG_INFO(logger, "client: rank: {} chdir to {}", irank, path);

      if (chdir(path) == -1) {
        LOG_ERROR(logger, "error: chdir: {}", strerror(errno));
      }
    }
#endif

    printSystemInfo();

    // XXX(rg): Be more gentle here
    bool bundlingEnabled = false;
    int bundleSize = -1; // getBundleSize();
    if (bundleSize == 0) {
      bundleSize = 1;
    } else if (bundleSize == -1) {
      // Not using bundling
      bundleSize = 1;
      bundlingEnabled = false;
    }

    std::vector<std::string> bundledFilenames;
    for (int i = 0; i < bundleSize; i++) {
      if (bundleSize > 1)
        LOG_INFO(logger, "Beginning Job {} of {}", i + 1, bundleSize);
      std::vector<std::string> unbundledFilenames;
      if (bundlingEnabled) {
        unbundledFilenames = unbundle(i);
      }

      // check to see if parameters file exists before loading
      int error = 0;
      string config_file = helper_functions::getRelevantFile(
          parameters.main_options.iniFilename);
      LOG_INFO(logger, "Loading parameter file {}", config_file);
      error = parameters.load(config_file);

      if (error) {
        LOG_ERROR(logger, "problem loading parameter file, stopping");
        logger->flush_log();
        exit(1);
        abort();
      }

      // Determine what type of job we are running according to the parameters
      // file.
      auto job =
          helper_functions::makeJob(std::make_unique<Parameters>(parameters));
      if (job == nullptr) {
        LOG_ERROR(logger, "error: Unknown job: {}",
                  std::string{magic_enum::enum_name<JobType>(
                      parameters.main_options.job)});
        logger->flush_log();
        return 1;
      }

      std::vector<std::string> filenames;
      try {
        filenames = job->run();
      } catch (int e) {
        LOG_CRITICAL(logger, "[ERROR] job exited on error {}", e);
        logger->flush_log();
      } catch (const std::exception &e) {
        LOG_CRITICAL(logger, "[ERROR] unhandled exception: {}", e.what());
        logger->flush_log();
        return EXIT_FAILURE;
      }

      filenames.push_back(std::string("client.log"));

      // Finalize Timing Information
      auto end_time = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed = end_time - start_time;

      double utime = 0, stime = 0, rtime = 0;
      helper_functions::getTime(&rtime, &utime, &stime);

      LOG_INFO(logger, "Timing Information:");
      LOG_INFO(logger, "  Real time: {:.3f} seconds", elapsed.count());
      LOG_INFO(logger, "  User time: {:.3f} seconds", utime);
      LOG_INFO(logger, "  System time: {:.3f} seconds", stime);

      std::ofstream result_file("results.dat", std::ios::app);
      if (result_file.is_open()) {
        result_file << "time_seconds " << elapsed.count() << "\n";
#ifndef _WIN32
        result_file << "user_time " << utime << "\n";
        result_file << "system_time " << stime << "\n";
#endif
      } else {
        LOG_ERROR(logger, "Failed to write timing to results.dat");
      }

      if (bundlingEnabled) {
        bundle(i, filenames, &bundledFilenames);
        deleteUnbundledFiles(unbundledFilenames);
      } else {
        bundledFilenames = filenames;
      }
    }

#ifdef EONMPI
    if (client_standalone) {
      break;
    }
    MPI::COMM_WORLD.Isend(&path[0], 1024, MPI::CHAR, server_rank, 0);

    // End of MPI while loop
  }
#endif

#ifdef OSX
#ifndef __aarch64__
  print_memory_usage();
#endif
#endif

#ifdef EONMPI
  if (client_standalone) {
    MPI::COMM_WORLD.Abort(0);
  } else {
    MPI::Finalize();
  }
#endif

  // Ensure all queued log messages are flushed before exiting
  quill::Backend::stop();
  exit(0);
}
