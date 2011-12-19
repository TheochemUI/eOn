#include "LAMMPS_EON.h"
#include <map>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "string.h"
#include "library.h"

// General Functions
lammps_eon::lammps_eon(Parameters *p){
    parameters = p;
    numberOfAtoms = 0;
    box0 = box4 = box8 = 0;
    LAMMPSObj=NULL;
}

void lammps_eon::cleanMemory(void){
    if(LAMMPSObj != NULL){
        lammps_close(LAMMPSObj);
        LAMMPSObj = NULL;
    }
    return;
}

void lammps_eon::force(long N, const double *R, const int *atomicNrs,
                       double *F, double *U, const double *box){

    int i;

    if (numberOfAtoms != N || box0 != box[0] || box4 != box[4] || 
        box8 != box[8]){
        makeNewLAMMPS(N, R, atomicNrs, box);
    }    

    lammps_put_coords(LAMMPSObj, (double *)R);
    lammps_command(LAMMPSObj,"run 0");  
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

    free(fx);
    free(fy);
    free(fz);
}

void lammps_eon::makeNewLAMMPS(long N, const double *R, const int *atomicNrs, const double *box){

    numberOfAtoms = N;
    box0 = box[0];
    box4 = box[4];
    box8 = box[8];
    
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

    char *lammps_argv[7];
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
    lammps_open_no_mpi(7, lammps_argv, &LAMMPSObj);
    void *ptr = LAMMPSObj;

    char cmd[80];

    //Gives units in Angstoms and eV
    lammps_command(ptr, "units metal");
    lammps_command(ptr, "atom_style	atomic");

    //Preserves atomic index ordering
    lammps_command(ptr, "atom_modify map array sort 0 0");


    //Define periodic cell
    snprintf(cmd, 80, "region box block 0 %.8f 0 %.8f 0 %.8f units box", 
             box0, box4, box8);    
    lammps_command(ptr, cmd);
    snprintf(cmd, 80, "create_box %i box", ntypes);
    lammps_command(ptr, cmd);

    //Initialize the atoms and their types
    for (int i=0;i<N;i++) {
        snprintf(cmd, 80, "create_atoms %i single %f %f %f units box", 
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
