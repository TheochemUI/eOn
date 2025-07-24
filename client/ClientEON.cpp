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
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <string.h>
#include <time.h>
#include <filesystem>

#ifdef EONMPI
#include <Python.h>
#include <fcntl.h>
#include <mpi.h>
#include <sstream>
#include <stdlib.h>
#endif

#if defined WITH_ASE_ORCA || EMBED_PYTHON || WITH_ASE_NWCHEM
#include <pybind11/embed.h>
#endif

#ifdef EONMPIBGP
#include <libgen.h>
#endif

//Includes for FPE trapping
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

  void print_memory_usage() {
      struct task_basic_info t_info;
      mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

      if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count)) {
          printf("Failed to get task info\n");
          return;
      }

      unsigned int rss = t_info.resident_size;
      unsigned int vs = t_info.virtual_size;
      printf(
          "\nmemory usage:\nresident size (MB): %8.2f\nvirtual size (MB):  %8.2f\n",
          (double)rss / 1024 / 1024, (double)vs / 1024 / 1024);
  }
  #endif
#endif

void printSystemInfo() {
  spdlog::info("EON Client");
  spdlog::info("VERSION: {}", VERSION);
  spdlog::info("BUILD DATE: {}\n", BUILD_DATE);
#ifndef __aarch64__
  spdlog::info("OS: {}", OS_INFO);
  spdlog::info("Arch: {}", ARCH);
#endif

#ifdef _WIN32
  TCHAR hostname[MAX_COMPUTERNAME_LENGTH + 1];
  DWORD size = sizeof(hostname) / sizeof(hostname[0]);
  if (GetComputerName(hostname, &size)) {
    spdlog::info("Hostname: {}", hostname);
  } else {
    spdlog::error("Failed to get hostname");
  }
  spdlog::info("PID: {}", GetCurrentProcessId());
#else
  struct utsname systemInfo;
  int status = uname(&systemInfo);
  if (status == 0) {
    spdlog::info("Hostname: {}", systemInfo.nodename);
    spdlog::info("PID: {}", getpid());
  } else {
    spdlog::error("Failed to get system information");
  }
#endif

  std::filesystem::path cwd = std::filesystem::current_path();
  spdlog::info("DIR: {}", cwd.string());
}

int main(int argc, char **argv) {
  // --- Start Logging setup
  // Sinks
  spdlog::flush_every(std::chrono::seconds(3));
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
      "client_spdlog.log", true); // Overwrite existing
  auto logger = std::make_shared<spdlog::logger>(
      "combi", spdlog::sinks_init_list({console_sink, file_sink}));
  spdlog::register_logger(logger);
  logger->set_pattern("%v");
  spdlog::set_default_logger(logger);
  // Traceback logger
  spdlog::set_level(spdlog::level::trace);
  auto trace_csink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  auto trace_fsink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
      "client_traceback.log", true); // Overwrite existing
  auto _traceback = std::make_shared<spdlog::logger>(
      "_traceback", spdlog::sinks_init_list({trace_csink, trace_fsink}));
  _traceback->set_pattern("%^ [%l] [%s:%#] [%!] \n %v\n[end %l]");
  spdlog::register_logger(_traceback);
  //--- End logging setup
  Parameters parameters;

#if defined WITH_ASE_ORCA || EMBED_PYTHON || WITH_ASE_NWCHEM
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
  string config_file = "config.ini";
  if (client_standalone) {
    if (helper_functions::existsFile("config_0.ini")) {
      config_file = "config_0.ini";
    }
    printf("Loading parameter file %s\n", config_file.c_str());
    error = parameters.load(config_file);
  } else {
    printf("Loading parameter file %s\n", parameters.iniFilename.c_str());
    error = parameters.load(parameters.iniFilename);
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
    commandLine(argc, argv);
    return 0;
  }
#endif

    eonc::enableFPE();  // from ExceptionsEON.h

  double beginTime = 0.0;
  helper_functions::getTime(&beginTime, NULL, NULL);

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
      // Potential::fcalls = 0;
      // Potential::fcallsTotal = 0;
      if (bundleSize > 1)
        printf("Beginning Job %d of %d\n", i + 1, bundleSize);
      std::vector<std::string> unbundledFilenames;
      if (bundlingEnabled) {
        unbundledFilenames = unbundle(i);
      }

      // check to see if parameters file exists before loading
      int error = 0;
      string config_file =
          helper_functions::getRelevantFile(parameters.iniFilename);
      printf("Loading parameter file %s\n", config_file.c_str());
      error = parameters.load(config_file);

      if (error) {
        fprintf(stderr, "\nproblem loading parameter file, stopping\n");
        exit(1);
        abort();
      }

      // Determine what type of job we are running according to the parameters
      // file.
      auto job =
          helper_functions::makeJob(std::make_unique<Parameters>(parameters));
      if (job == nullptr) {
        printf("error: Unknown job: %s\n",
               std::string{magic_enum::enum_name<JobType>(parameters.job)}
                   .c_str());
        return 1;
      }

      std::vector<std::string> filenames;
      try {
        filenames = job->run();
      } catch (int e) {
        SPDLOG_CRITICAL("[ERROR] job exited on error %d\n", e);
      }

      filenames.push_back(std::string("client.log"));

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

  // Timing Information
  double utime = 0, stime = 0, rtime = 0;
  helper_functions::getTime(&rtime, &utime, &stime);
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

  exit(0);
}
