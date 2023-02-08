#include "LAMMPS.h"
#include <map>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "library.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

// General Functions
lammps::lammps(Parameters *p){
    parameters = p;
    numberOfAtoms = 0;
    LAMMPSObj=NULL;
    for (int i=0;i<9;i++) oldBox[i] = 0.0;
}

lammps::~lammps(){
    cleanMemory();
}

void lammps::cleanMemory(void){
    if(LAMMPSObj != NULL){
        lammps_close(LAMMPSObj);
        LAMMPSObj = NULL;
    }
    return;
}

void lammps::force(long N, const double *R, const int *atomicNrs,
                       double *F, double *U, const double *box){

    int i;
    bool newLammps=false;
    for (int i=0;i<9;i++) {
        if (oldBox[i] != box[i]) newLammps = true;
    }
    if (numberOfAtoms != N) newLammps = true;
    if (newLammps) {
        makeNewLAMMPS(N, R, atomicNrs, box);
    }

    lammps_scatter_atoms(LAMMPSObj, "x", 1, 3, R);
    lammps_command(LAMMPSObj,"run 1 pre no post no");

    double *pe = (double *)lammps_extract_variable(LAMMPSObj, "pe", NULL);
    *U = *pe;
    free(pe);

    double *fx = (double *)lammps_extract_variable(LAMMPSObj, "fx", "all");
    double *fy = (double *)lammps_extract_variable(LAMMPSObj, "fy", "all");
    double *fz = (double *)lammps_extract_variable(LAMMPSObj, "fz", "all");

    for (i=0;i<N;i++) {
        F[3*i+0] = fx[i];
        F[3*i+1] = fy[i];
        F[3*i+2] = fz[i];
    }

    // convert from kCal/mol -> eV if LAMMPS is using real units

    if(realunits){
        *U = *U / 23.0609;
        for (i=0;i<3*N;i++) {
            F[i] = F[i] / 23.0609;
        }
    }

    free(fx);
    free(fy);
    free(fz);
}

void lammps::makeNewLAMMPS(long N, const double *R, const int *atomicNrs, const double *box){
    numberOfAtoms = N;
    for (int i=0;i<9;i++) oldBox[i] = box[i];

    if(LAMMPSObj != NULL){
        cleanMemory();
    }

    std::map<int,int> type_map;
    int ntypes=0;
    for (int i=0;i<N;i++) {
        if (type_map.count(atomicNrs[i]) == 0) {
            ntypes += 1;
            type_map.insert(std::pair<int,int>(atomicNrs[i],ntypes));
        }
    }

    char *lammps_argv[9];
    int nargs=7;
    lammps_argv[0] = "";
    lammps_argv[1] = "-log";
    if (parameters->LAMMPSLogging) {
        lammps_argv[2] = "log.lammps";
    }else{
        lammps_argv[2] = "none";
    }
    lammps_argv[3] = "-echo";
    lammps_argv[4] = "log";
    lammps_argv[5] = "-screen";
    lammps_argv[6] = "none";


    if (parameters->LAMMPSThreads > 0) {
        lammps_argv[7] = "-suffix";
        lammps_argv[8] = "omp";
        nargs += 2;
    }
    #ifdef EONMPI
        lammps_open(nargs, lammps_argv, parameters->MPIClientComm, &LAMMPSObj);
    #else
        lammps_open_no_mpi(nargs, lammps_argv, &LAMMPSObj);
    #endif
    void *ptr = LAMMPSObj;

    char cmd[200];

    if (parameters->LAMMPSThreads > 0) {
        snprintf(cmd, 200, "package omp %i force/neigh", parameters->LAMMPSThreads);
        lammps_command(ptr, cmd);
    }

    //Gives units in Angstoms and eV
//    lammps_command(ptr, "units metal");

    // We need to allow for 'real' units for reaxff
    realunits = false;
    FILE *file;
    file = fopen("in.lammps", "r");
    if (!file) {
        fprintf(stderr, "couldn't open in.lammps: %s\n",strerror(errno));
        return;
    }
    char line[256];
    while (fgets(line, sizeof(line), file)) {
      if (strcmp(line, "#!units real\n") == 0){
         realunits = true;
      }
    }
    fclose(file);

    if (realunits){
        lammps_command(ptr, "units real");
    } else {
        lammps_command(ptr, "units metal");
    }

    lammps_command(ptr, "atom_style	charge");

    //Preserves atomic index ordering
    lammps_command(ptr, "atom_modify map array sort 0 0");

    //Always check to see if the neighbor list must be updated
    lammps_command(ptr, "neigh_modify delay 1");


    //Define periodic cell
    //    prism args = xlo xhi ylo yhi zlo zhi xy xz yz
    //        xlo,xhi,ylo,yhi,zlo,zhi = bounds of untilted prism
    //        xy = distance to tilt y in x direction
    //        xz = distance to tilt z in x direction
    //        yz = distance to tilt z in y direction
    snprintf(cmd, 200, "region cell prism 0 %f 0 %f 0 %f %f %f %f units box",
             box[0], box[4], box[8], box[3], box[6], box[7]);

    lammps_command(ptr, cmd);
    snprintf(cmd, 200, "create_box %i cell", ntypes);
    lammps_command(ptr, cmd);

    //Initialize the atoms and their types
    for (int i=0;i<N;i++) {
        snprintf(cmd, 200, "create_atoms %i single %f %f %f units box",
                 type_map[atomicNrs[i]], 0.0, 0.0, 0.0);
        lammps_command(ptr, cmd);
    }

    //We don't care about mass but have to set it
    lammps_command(ptr, "mass * 1.0");

    //Read in user commands from in.lammps files
    struct stat buffer;
    if (stat("in.lammps", &buffer) == -1) {
        fprintf(stderr, "couldn't open in.lammps: %s\n",strerror(errno));
        exit(1);
    }else{
        lammps_file(ptr, "in.lammps");
    }

    //Define variables for force and energy so they can be extracted
    lammps_command(ptr, "variable fx atom fx");
    lammps_command(ptr, "variable fy atom fy");
    lammps_command(ptr, "variable fz atom fz");
    lammps_command(ptr, "variable pe equal pe");
}
