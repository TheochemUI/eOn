#include "LAMMPSPot.h"
#include <map>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// LAMMPS library inclusion
#ifdef EONMPI
#define LAMMPS_LIB_MPI
#endif
#include "library.h"

LAMMPSPot::LAMMPSPot(std::shared_ptr<Parameters> p)
    : Potential(p) {
  numberOfAtoms = 0;
  LAMMPSObj = nullptr;
  std::fill_n(oldBox, 9, 0.0);
}

LAMMPSPot::~LAMMPSPot() { cleanMemory(); }

void LAMMPSPot::cleanMemory(void) {
  if (LAMMPSObj != nullptr) {
    lammps_close(LAMMPSObj);
    LAMMPSObj = nullptr;
  }
  return;
}

void LAMMPSPot::force(long N, const double *R, const int *atomicNrs, double *F,
                      double *U, double *variance, const double *box) {
  variance = nullptr;
  int i;
  bool newLammps = false;
  for (int i = 0; i < 9; i++) {
    if (oldBox[i] != box[i])
      newLammps = true;
  }
  if (numberOfAtoms != N)
    newLammps = true;
  if (newLammps) {
    makeNewLAMMPS(N, R, atomicNrs, box);
  }
  if (!LAMMPSObj) {
    throw std::runtime_error("Should have a LAMMPS instance by now");
  }

  lammps_scatter_atoms(LAMMPSObj, "x", 1, 3, const_cast<double *>(R));
  lammps_command(LAMMPSObj, "run 1 pre no post no");

  double *pe = (double *)lammps_extract_variable(LAMMPSObj, "pe", nullptr);
  *U = *pe;
  free(pe);

  double *fx = (double *)lammps_extract_variable(LAMMPSObj, "fx", "all");
  double *fy = (double *)lammps_extract_variable(LAMMPSObj, "fy", "all");
  double *fz = (double *)lammps_extract_variable(LAMMPSObj, "fz", "all");

  for (i = 0; i < N; i++) {
    F[3 * i + 0] = fx[i];
    F[3 * i + 1] = fy[i];
    F[3 * i + 2] = fz[i];
  }

  // convert from kCal/mol -> eV if LAMMPS is using real units

  if (realunits) {
    *U = *U / 23.0609;
    for (i = 0; i < 3 * N; i++) {
      F[i] = F[i] / 23.0609;
    }
  }

  free(fx);
  free(fy);
  free(fz);
}

void LAMMPSPot::makeNewLAMMPS(long N, const double *R, const int *atomicNrs,
                              const double *box) {
  numberOfAtoms = N;
  for (int i = 0; i < 9; i++)
    oldBox[i] = box[i];

  if (LAMMPSObj != nullptr) {
    cleanMemory();
  }

  std::map<int, int> type_map;
  int ntypes = 0;
  for (int i = 0; i < N; i++) {
    if (type_map.count(atomicNrs[i]) == 0) {
      ntypes += 1;
      type_map.insert(std::pair<int, int>(atomicNrs[i], ntypes));
    }
  }


#ifdef EONMPI
  // TODO(rg):  Suffix should be a parameter
  const char *lmpargv[] = { "liblammps", "-log", "none", "-echo", "log", "-screen", "none", "-suffix", "omp"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPSObj = lammps_open(lmpargc, (char **)lmpargv, m_params->MPIClientComm, nullptr);
#else
  const char *lmpargv[] = { "liblammps", "-log", "none", "-echo", "log", "-screen", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);
  LAMMPSObj = lammps_open_no_mpi(lmpargc, (char **)lmpargv, NULL);
#endif

  char cmd[200];

  if (m_params->LAMMPSThreads > 0) {
    snprintf(cmd, 200, "package omp %i force/neigh", m_params->LAMMPSThreads);
    lammps_command(LAMMPSObj, cmd);
  }

  // Gives units in Angstoms and eV
  //    lammps_command(LAMMPSObj, "units metal");

  // We need to allow for 'real' units for reaxff
  realunits = false;
  FILE *file;
  file = fopen("in.lammps", "r");
  if (!file) {
    fprintf(stderr, "couldn't open in.lammps: %s\n", strerror(errno));
    return;
  }
  char line[256];
  while (fgets(line, sizeof(line), file)) {
    if (strcmp(line, "#!units real\n") == 0) {
      realunits = true;
    }
  }
  fclose(file);

  if (realunits) {
    lammps_command(LAMMPSObj, "units real");
  } else {
    lammps_command(LAMMPSObj, "units metal");
  }

  lammps_command(LAMMPSObj, "atom_style	charge");

  // Preserves atomic index ordering
  lammps_command(LAMMPSObj, "atom_modify map array sort 0 0");

  // Always check to see if the neighbor list must be updated
  lammps_command(LAMMPSObj, "neigh_modify delay 1");

  // Define periodic cell
  //     prism args = xlo xhi ylo yhi zlo zhi xy xz yz
  //         xlo,xhi,ylo,yhi,zlo,zhi = bounds of untilted prism
  //         xy = distance to tilt y in x direction
  //         xz = distance to tilt z in x direction
  //         yz = distance to tilt z in y direction
  snprintf(cmd, 200, "region cell prism 0 %f 0 %f 0 %f %f %f %f units box",
           box[0], box[4], box[8], box[3], box[6], box[7]);

  lammps_command(LAMMPSObj, cmd);
  snprintf(cmd, 200, "create_box %i cell", ntypes);
  lammps_command(LAMMPSObj, cmd);

  // Initialize the atoms and their types
  for (int i = 0; i < N; i++) {
    snprintf(cmd, 200, "create_atoms %i single %f %f %f units box",
             type_map[atomicNrs[i]], 0.0, 0.0, 0.0);
    lammps_command(LAMMPSObj, cmd);
  }

  // We don't care about mass but have to set it
  lammps_command(LAMMPSObj, "mass * 1.0");

  // Read in user commands from in.lammps files
  struct stat buffer;
  if (stat("in.lammps", &buffer) == -1) {
    fprintf(stderr, "couldn't open in.lammps: %s\n", strerror(errno));
    exit(1);
  } else {
    lammps_file(LAMMPSObj, "in.lammps");
  }

  // Define variables for force and energy so they can be extracted
  lammps_command(LAMMPSObj, "variable fx atom fx");
  lammps_command(LAMMPSObj, "variable fy atom fy");
  lammps_command(LAMMPSObj, "variable fz atom fz");
  lammps_command(LAMMPSObj, "variable pe equal pe");
}
