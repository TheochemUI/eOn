#include "stacktrace.h"

#ifdef STACKTRACE
#include <Python.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#if STACKTRACE == alpha
#define USE_LADEBUG
#else
#define USE_GDB
#endif

/* static char mpipython[] = "/usr/local/bin/mpipython"; */

static char mpipython[] = "/usr/local/lam-7.0.6-gcc/bin/mpipython";
static char *execname;
static char crashfilename[100];
static abortfunc *abortfunction;

static void sig_handler(int sig, siginfo_t *sip, void *extra);
static void sethandler(int signal);

void setSignalHandlers(int node, abortfunc *ab) {
  abortfunction = ab;
  if (node >= 0) {
    sprintf(crashfilename, "CRASH.%d", node);
    execname = mpipython; /* Py_GetProgramFullPath() is not realiable */
  } else {
    strcpy(crashfilename, "CRASH");
    execname = Py_GetProgramFullPath();
  }
  sethandler(SIGBUS);
  sethandler(SIGILL);
  sethandler(SIGFPE);
  sethandler(SIGSEGV);
}

static void sethandler(int signal) {
  struct sigaction sa;
  int ret;
  sa.sa_sigaction = sig_handler;
  /* sa.sa_handler = sig_handler; */
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_SIGINFO;
  ret = sigaction(signal, &sa, 0);
  if (ret) {
    perror("sigaction");
    exit(1);
  }
}

static void sig_handler(int sig, siginfo_t *sip, void *extra) {
  char pid[100];
  int fpid, i;
  FILE *debug;
  FILE *crash;
  char debugname[100];
  char *name;
  char *r;

  switch (sig) {
  case SIGBUS:
    name = "SIGBUS: Bus error";
    break;
  case SIGCHLD:
    name = "SIGCHLD: Child process event";
    break;
  case SIGILL:
    name = "SIGILL: Illegal operation";
    break;
  case SIGFPE:
    name = "SIGFPE: Arithmetic error";
    break;
  case SIGSEGV:
    name = "SIGSEGV: Segmentation fault";
    break;
  default:
    name = "Unknown signal";
  }

  fprintf(stderr, "\n\nASAP signal handler: An unexpected error occurred.\n");
  fprintf(stderr,
          "Process %s got signal %d (%s).\n  si_code = %d.  si_addr = 0x%lx\n",
          execname, sig, name, sip->si_code, (long)sip->si_addr);

  crash = fopen(crashfilename, "w");
  fprintf(
      crash,
      "\nProcess %s got signal %d (%s).\n  si_code = %d.  si_addr = 0x%lx\n",
      execname, sig, name, sip->si_code, (long)sip->si_addr);

  r = "Unknown reason";

  switch (sig) {
  case SIGBUS:
    switch (sip->si_code) {
    case BUS_ADRALN:
      r = "invalid address alignment";
      break;
    case BUS_ADRERR:
      r = "non-existent physical address";
      break;
    case BUS_OBJERR:
      r = "object specific hardware error";
      break;
    }
    break;
  case SIGCHLD:
    switch (sip->si_code) {
    case CLD_EXITED:
      r = "child has exited";
      break;
    case CLD_KILLED:
      r = "child was killed";
      break;
    case CLD_DUMPED:
      r = "child terminated abnormally";
      break;
    case CLD_TRAPPED:
      r = "traced child has trapped";
      break;
    case CLD_STOPPED:
      r = "child has stopped";
      break;
    case CLD_CONTINUED:
      r = "stopped child has continued";
      break;
#if STACKTRACE == alpha
    case CLD_SIGEXITING:
      r = "child is about to exit because it received a fatal signal";
      break;
#endif
    }
    break;
  case SIGILL:
    switch (sip->si_code) {
    case ILL_ILLOPC:
      r = "illegal opcode";
      break;
    case ILL_ILLOPN:
      r = "illegal operand";
      break;
    case ILL_ILLADR:
      r = "illegal addressing mode";
      break;
    case ILL_ILLTRP:
      r = "illegal trap";
      break;
    case ILL_PRVOPC:
      r = "privileged opcode";
      break;
    case ILL_PRVREG:
      r = "privileged register";
      break;
    case ILL_COPROC:
      r = "coprocessor error";
      break;
    case ILL_BADSTK:
      r = "internal stack error";
      break;
    }
    break;
  case SIGFPE:
    switch (sip->si_code) {
    case FPE_INTDIV:
      r = "integer divide by zero";
      break;
    case FPE_INTOVF:
      r = "integer overflow";
      break;
    case FPE_FLTDIV:
      r = "floating point divide by zero";
      break;
    case FPE_FLTOVF:
      r = "floating point overflow";
      break;
    case FPE_FLTUND:
      r = "floating point underflow";
      break;
    case FPE_FLTRES:
      r = "floating point inexact result";
      break;
    case FPE_FLTINV:
      r = "invalid floating point operation";
      break;
    case FPE_FLTSUB:
      r = "subscript out of range";
      break;
    }
    break;
  case SIGSEGV:
    switch (sip->si_code) {
    case SEGV_MAPERR:
      r = "address not mapped to object";
      break;
    case SEGV_ACCERR:
      r = "invalid permissions for mapped object";
      break;
    }
    break;
  }
  fprintf(stderr, "Cause of signal: %s\n\n", r);
  fprintf(crash, "Cause of signal: %s\n\n", r);
  fprintf(stderr, "\nA traceback will be written to the file %s\n",
          crashfilename);

  fprintf(stderr,
          "Please email this file to the programmers for debugging usage.\n\n");

  sprintf(pid, "%d", getpid());
#ifdef USE_LADEBUG
  sprintf(debugname, "/tmp/dump.ladebug.%d", getpid());
#else
  sprintf(debugname, "/tmp/dump.gdb.%d", getpid());
#endif
  debug = fopen(debugname, "w");
  if (debug == 0) {
    perror("fopen failed");
    exit(3);
  }
  fprintf(debug, "where\ndetach\nquit\n");
  fclose(debug);
  fflush(crash);
  fpid = fork();
  if (fpid == -1) {
    perror("fork failed");
    exit(1);
  }
  if (fpid == 0) {
    /* The child */
    dup2(fileno(crash), fileno(stderr));
    dup2(fileno(crash), fileno(stdout));
#if defined(USE_LADEBUG)
    execlp("ladebug", "ladebug", "-pid", pid, "-c", debugname, execname, 0);
#elif defined(USE_IDB_GDB)
    execlp("idb", "idb", "-gdb", "-pid", pid, "-c", debugname, execname, 0);
#else
    execlp("gdb", "gdb", "-x", debugname, execname, pid, 0);
#endif
  } else {
    sleep(5);
#ifdef USE_LADEBUG
    fprintf(stderr, "Transferring control to debugger using SIGINT\n");
    kill(getpid(), SIGINT);
#endif
    sleep(3);
    unlink(debugname);
    if (abortfunction) {
      fprintf(stderr, "Exiting using abort function (MPI)....\n");
      abortfunction();
    } else {
      fprintf(stderr, "Exiting using exit....\n");
      exit(3);
    }
  }
  fprintf(stderr,
          "\n\nReached end of signal handler: This should not happen!!!!\n");
  exit(42);
}

#endif /* STACKTRACE */
