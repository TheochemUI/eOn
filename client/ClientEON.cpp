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
#include "client/matter/MatterHelpers.hpp"
#include "version.h"

#include <cstdlib>
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

#ifdef WITH_ASE_ORCA
#include <pybind11/embed.h>
#endif

#ifdef EONMPIBGP
#include <libgen.h>
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
  SPDLOG_LOGGER_INFO("\nmemory usage:\nresident size (MB): {:8.2f}\nvirtual "
                     "size (MB):  {:8.2f}",
                     static_cast<double>(rss / 1024 / 1024),
                     static_cast<double>(vs / 1024 / 1024));
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

#if defined WITH_ASE_ORCA || WITH_PYTHON
  pybind11::scoped_interpreter guard{};
#endif

#ifdef EONMPI
  bool client_standalone = false;
  if (getenv("EON_CLIENT_STANDALONE") != NULL) {
    client_standalone = true;
  }
  int number_of_clients;
  if (!client_standalone) {
    if (getenv("EON_SERVER_PATH") == NULL) {
      fprintf(stderr, "error: must set the env var EON_SERVER_PATH\n");
      return 1;
    }
    if (getenv("EON_NUMBER_OF_CLIENTS") == NULL) {
      fprintf(stderr, "error: must set the env var EON_NUMBER_OF_CLIENTS\n");
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
  std::string config_file = "config.toml";
  if (client_standalone) {
    if (helper_functions::existsFile("config_0.toml")) {
      config_file = "config_0.toml";
    }
    printf("Loading parameter file %s\n", config_file.c_str());
    error = parameters.load(config_file);
  } else {
    printf("Loading parameter file %s\n", parameters.inpFilename.c_str());
    error = parameters.load(parameters.inpFilename);
  }
  if (error) {
    fprintf(stderr, "\nproblem loading parameter file\n");
    MPI::COMM_WORLD.Abort(1);
  }
  printf("\n");

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
    fprintf(stderr, "didn't launch as many mpi client ranks as"
                    "specified in EON_NUMBER_OF_CLIENTS\n");
    MPI::COMM_WORLD.Abort(1);
  }
  clients = number_of_clients;

  if (parameters.potential == "mpi") {
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
      parameters.MPIPotentialRank =
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
      parameters.MPIClientComm = new_comm;
    }
    printf("creating group with ranks: %i\n", r);
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
      /* GH
      char **py_argv = (char **)malloc(sizeof(char **)*2);
      py_argv[0] = argv[0];
      py_argv[1] = getenv("EON_SERVER_PATH");
      fprintf(stderr, "rank: %i becoming %s\n", irank, py_argv[1]);
      Py_Initialize();
      Py_Main(2, py_argv);
      Py_Finalize();
      */
      wchar_t **py_argv = (wchar_t **)malloc(sizeof(wchar_t *) * 2);
      py_argv[0] = Py_DecodeLocale(argv[0], NULL);
      char *program = getenv("EON_SERVER_PATH");
      py_argv[1] = Py_DecodeLocale(program, NULL);
      fprintf(stderr, "rank: %i becoming %s\n", irank, program);
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
    eonc::commandLine(argc, argv);
    return 0;
  }
#endif

  eonc::enableFPE(); // from ExceptionsEON.h

  double beginTime = 0.0;
  eonc::helper_functions::getTime(&beginTime, NULL, NULL);

#ifdef EONMPI
  // XXX: When do we stop? The server should probably tell everyone when to
  // stop.
  char logfilename[1024];
  snprintf(logfilename, 1024, "eonclient_%i.log", my_client_number);
  // int outFd = open("/dev/null", O_WRONLY);

  if (!client_standalone) {
    //    int outFd = open(logfilename, O_WRONLY|O_CREAT|O_TRUNC, 0644);
    //    dup2(outFd, 1);
    //    dup2(outFd, 2);
  }
  char *orig_path = new char[1024];
  getcwd(orig_path, 1024);
  while (true) {
    chdir(orig_path);
    char *path = new char[1024];
    int ready = 1;
    if (!client_standalone) {
      fprintf(stderr,
              "client: rank %i is ready, posting send to server rank: %i!\n",
              irank, server_rank);
      // Tag "1" is to interrupt the main loop and tell the communicator that a
      // client is ready
      MPI::COMM_WORLD.Isend(&ready, 1, MPI::INT, server_rank, 1);

      // Get the path we should run in from the server
      MPI::COMM_WORLD.Recv(&path[0], 1024, MPI::CHAR, server_rank, 0);
      if (strncmp("STOPCAR", path, 1024) == 0) {
        fprintf(stderr, "rank %i got STOPCAR\n", irank);
        MPI::Finalize();
        return 0;
      }
      fprintf(stderr, "client: rank: %i chdir to %s\n", irank, path);

      if (chdir(path) == -1) {
        fprintf(stderr, "error: chdir: %s\n", strerror(errno));
      }
    }
#endif

    printSystemInfo(logger);

    // check to see if parameters file exists before loading
    std::string config_file =
        eonc::helper_functions::getRelevantFile(parameters.main.inpFilename);
    SPDLOG_LOGGER_INFO(logger, "Loading parameter file {}", config_file);
    auto params = eonc::loadTOML(config_file);

    auto CACHELOT_CMD = cachelot::cache::Cache::Create(
        eonc::cache::cache_memory, eonc::cache::page_size,
        eonc::cache::hash_initial, true);
    auto pcache = eonc::cache::PotentialCache();
    pcache.set_cache(&CACHELOT_CMD);
    auto pot = eonc::makePotential(params);
    pot->set_cache(&pcache);

    auto mats = eonc::mat::make_matter(params, pot);
    auto job = eonc::mkJob(params);

    bool result = eonc::JobRunner(job, mats[0]);
    if (!result) {
      throw std::runtime_error("Something went wrong running the job");
    }

#ifdef EONMPI
    if (client_standalone) {
      break;
    }
    MPI::COMM_WORLD.Isend(&path[0], 1024, MPI::CHAR, server_rank, 0);

    // End of MPI while loop
  }
#endif

  // Timing Information
  double utime = 0, stime = 0, rtime = 0;
  eonc::helper_functions::getTime(&rtime, &utime, &stime);
  rtime = rtime - beginTime;

  // if (Potential::totalUserTime > 0) {
  //     printf("\ntime not in potential: %.4f%%\n",
  //     100*(1-Potential::totalUserTime/utime));
  // }

  printf("timing information:\nreal %10.3f seconds\nuser %10.3f seconds\nsys  "
         "%10.3f seconds\n",
         rtime, utime, stime);

#ifdef OSX
#ifndef __aarch64__
  // TODO(rg) :: Either make it more generic or just remove it..
  print_memory_usage(logger);
#endif
#endif

#ifdef EONMPI
  if (client_standalone) {
    MPI::COMM_WORLD.Abort(0);
  } else {
    MPI::Finalize();
  }
#endif

  std::exit(0);
}
