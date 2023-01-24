#include <time.h>

#include "Potential.h"
#include "Parameters.h"
#include "Log.h"
#include "HelperFunctions.h"

#ifdef IMD_POT
    #include "potentials/IMD/IMD.h"
#endif

#ifdef WITH_GPRD
    #include "potentials/GPRPotential/GPRPotential.h"
#endif

#include "potentials/EMT/EffectiveMediumTheory.h"
#include "potentials/ExtPot/ExtPot.h"
#include "potentials/Morse/Morse.h"
#include "potentials/LJ/LJ.h"
#include "potentials/LJCluster/LJCluster.h"
#include "potentials/EAM/EAM.h"

#ifdef WITH_FORTRAN
  #include "potentials/Aluminum/Aluminum.h"
  #include "potentials/EDIP/EDIP.h"
  #include "potentials/FeHe/FeHe.h"
  #include "potentials/Lenosky/Lenosky.h"
  #include "potentials/SW/SW.h"
  #include "potentials/Tersoff/Tersoff.h"
#endif

#ifdef WITH_PYTHON
 #ifdef PYAMFF_POT
  #include "potentials/PyAMFF/PyAMFF.h"
 #endif
  #include "potentials/QSC/QSC.h"
#endif

#ifdef EONMPI
    #include "potentials/MPIPot/MPIPot.h"
#endif

#ifdef LAMMPS_POT
    #include "potentials/LAMMPS/LAMMPS.h"
#endif

#ifdef NEW_POT
    #include "potentials/NewPot/NewPot.h"
#endif

// TODO: This should be guarded by WITH_FORTRAN as well
#ifdef CUH2_POT
    #include "potentials/CuH2/CuH2.h"
#endif

#ifndef WIN32
#ifdef WITH_VASP
    #include "potentials/VASP/VASP.h"
#endif
#endif

#ifdef WITH_AMS
    #include "potentials/AMS/AMS.h"
    #include "potentials/AMS_IO/AMS_IO.h"
#endif

#ifdef WITH_WATER
  #include "potentials/Water/Water.hpp"
 #ifdef WITH_FORTRAN
  #include "potentials/Water_H/Tip4p_H.h"
 #endif
  #include "potentials/Water_Pt/Tip4p_Pt.hpp"
#endif

#include <cstdlib>

const std::string Potential::POT_EMT =         "emt"s;
const std::string Potential::POT_EXT =         "ext_pot"s;
const std::string Potential::POT_LJ =          "lj"s;
const std::string Potential::POT_LJCLUSTER =   "lj_cluster"s;
const std::string Potential::POT_MORSE_PT =    "morse_pt"s;
const std::string Potential::POT_NEW =         "new_pot"s;

#ifdef CUH2_POT
const char Potential::POT_CUH2[] =        "cuh2_pot";
#endif

#ifdef IMD_POT
const std::string Potential::POT_IMD =         "imd"s;
#endif

#ifdef WITH_GPRD
const std::string Potential::POT_GPR =         "gpr_pot"s;
#endif

#ifdef WITH_WATER
const std::string Potential::POT_TIP4P =       "tip4p"s;
const std::string Potential::POT_TIP4P_PT =    "tip4p_pt"s;
#ifdef WITH_FORTRAN
const std::string Potential::POT_TIP4P_H =     "tip4p_h"s;
#endif
const std::string Potential::POT_SPCE =        "spce"s;
#endif

#ifdef WITH_FORTRAN
const std::string Potential::POT_EAM_AL =      "eam_al"s;
const std::string Potential::POT_EDIP =        "edip"s;
const std::string Potential::POT_FEHE =        "fehe"s;
const std::string Potential::POT_LENOSKY_SI =  "lenosky_si"s;
const std::string Potential::POT_SW_SI =       "sw_si"s;
const std::string Potential::POT_TERSOFF_SI =  "tersoff_si"s;
#endif

#ifdef LAMMPS_POT
const std::string Potential::POT_LAMMPS =      "lammps"s;
#endif

#ifdef EONMPI
const std::string Potential::POT_MPI =         "mpi"s;
#endif

#ifdef WITH_PYTHON
 #ifdef PYAMFF_POT
const std::string Potential::POT_PYAMFF =      "pyamff"s;
 #endif
const std::string Potential::POT_QSC =         "qsc"s;
#endif

#ifdef WITH_AMS
const std::string Potential::POT_AMS =         "ams"s;
const std::string Potential::POT_AMS_IO =      "ams_io"s;
#endif

#ifdef WITH_VASP
const std::string Potential::POT_VASP =        "vasp"s;
#endif

Potential* Potential::pot = NULL;

Potential *Potential::getPotential(Parameters *parameters)
{
    if(pot) {
        return pot;
    }
    if(parameters->potential == PotentialStrings::POT_LJ)
        pot = new LJ();
    else if(parameters->potential == PotentialStrings::POT_LJCLUSTER)
        pot = new LJCluster();
    else if (parameters->potential == PotentialStrings::POT_EXT)
        pot = new ExtPot(parameters);
    else if (parameters->potential == PotentialStrings::POT_MORSE_PT)
        pot = new Morse();
    else if (parameters->potential == PotentialStrings::POT_EMT)
        pot = new EffectiveMediumTheory(parameters);

#ifdef IMD_POT
  else if (parameters->potential == PotentialStrings::POT_IMD)
    pot = new IMD();
#endif

#ifdef WITH_PYTHON
#ifdef PYAMFF_POT
    else if(parameters->potential == PotentialStrings::POT_PYAMFF)
        pot = new PyAMFF();
#endif
  else if (parameters->potential == PotentialStrings::POT_QSC)
    pot = new QSC();
#endif

#ifdef WITH_AMS
    else if(parameters->potential == PotentialStrings::POT_AMS)
        pot = new AMS(parameters);
    else if(parameters->potential == PotentialStrings::POT_AMS_IO)
        pot = new AMS_IO(parameters);
#endif

#ifdef WITH_WATER
  else if (parameters->potential == PotentialStrings::POT_TIP4P)
    pot = new Tip4p();
  else if (parameters->potential == PotentialStrings::POT_SPCE)
    pot = new SpceCcl();
  else if (parameters->potential == PotentialStrings::POT_TIP4P_PT)
    pot = new Tip4p_Pt();
  else if (parameters->potential == PotentialStrings::POT_TIP4P_H)
    pot = new Tip4p_H();
#endif

#ifdef WITH_FORTRAN
  else if (parameters->potential == PotentialStrings::POT_EAM_AL)
    pot = new Aluminum();
  else if (parameters->potential == PotentialStrings::POT_LENOSKY_SI)
    pot = new Lenosky();
  else if (parameters->potential == PotentialStrings::POT_SW_SI)
    pot = new SW();
  else if (parameters->potential == PotentialStrings::POT_TERSOFF_SI)
    pot = new Tersoff();
  else if (parameters->potential == PotentialStrings::POT_EDIP)
    pot = new EDIP();
  else if (parameters->potential == PotentialStrings::POT_FEHE)
    pot = new FeHe();
#endif

#ifdef EONMPI
    else if(parameters->potential == PotentialStrings::POT_MPI)
        pot = new MPIPot(parameters);
#endif

#ifdef LAMMPS_POT
    else if(parameters->potential == PotentialStrings::POT_LAMMPS)
        pot = new lammps(parameters);
#endif

#ifdef NEW_POT
    else if(parameters->potential == PotentialStrings::POT_NEW)
        pot = new NewPot(parameters);
#endif

#ifdef CUH2_POT
    else if(parameters->potential == POT_CUH2)
        pot = new CuH2(parameters);
#endif

#ifndef WIN32
#ifdef WITH_VASP
    else if(parameters->potential == PotentialStrings::POT_VASP)
        pot = new VASP();
#endif
#endif

    else {
        printf("Unknown Potential: %s\n", parameters->potential.c_str());
        std::exit(1);
    }
    pot->initialize();
    pot->params = parameters;
    return pot;
}

int Potential::fcalls = 0;
int Potential::fcallsTotal = 0;
double Potential::totalUserTime=0;

AtomMatrix Potential::force(long nAtoms, AtomMatrix positions,
                            VectorXi atomicNrs, double *energy, Matrix3d box, int nImages)
{
    AtomMatrix forces(nAtoms,3);

    double start, userStart, sysStart;
    if (params->LogPotential) {
        helper_functions::getTime(&start, &userStart, &sysStart);
    }
      int *atnrs = atomicNrs.data();
      double *pos = positions.data();
      double *frcs = forces.data();
      double *bx = box.data();
    // TODO: Be better with the number of images
    force(nAtoms, pos, atnrs, frcs, energy,
          bx,1);

    double finish, userFinish, sysFinish;
    if (params->LogPotential) {
        helper_functions::getTime(&finish, &userFinish, &sysFinish);

        log_file("[Potential] fcall#: %4d  real: %.6e  user: %.6e  sys: %.6e seconds\n",
                 fcalls, finish - start, userFinish - userStart, sysFinish - sysStart);
        totalUserTime += userFinish - userStart;
    }

    fcalls += 1;
    fcallsTotal += 1;
//    printf("nAtoms %f\n",nAtoms);
//    printf("AtomMatrix positions %f\n",positions);
//    printf("VectorXi atomicNmrs %f\n",atomicNrs);
//    printf("energy %f\n",energy);
//    printf("Matrix3d box %f\n",box);

    if (params->maxForceCalls != 0) {
        if (fcallsTotal > params->maxForceCalls) {
            throw 1017;
        }
    }

    return forces;
}

void Potential::setParams(Parameters *params){
  this->params = params;
}
