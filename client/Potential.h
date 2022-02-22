#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Parameters.h"
#include "Eigen.h"

class Potential
{

    public:

        Potential(){}
        virtual ~Potential(){}

        static const std::string POT_EMT;
        static const std::string POT_EXT;
        static const std::string POT_LJ;
        static const std::string POT_LJCLUSTER;
        static const std::string POT_MORSE_PT;
        static const std::string POT_NEW;

#ifdef CUH2_POT
        static const char POT_CUH2[];
#endif

#ifdef IMD_POT
        static const std::string POT_IMD;
#endif

#ifdef WITH_GPRD
        static const std::string POT_GPR;
#endif

#ifdef WITH_WATER
        static const std::string POT_TIP4P;
        static const std::string POT_TIP4P_PT;
 #ifdef WITH_FORTRAN
        static const std::string POT_TIP4P_H;
 #endif
        static const std::string POT_SPCE;
#endif

#ifdef WITH_FORTRAN
        static const std::string POT_EAM_AL;
        static const std::string POT_EDIP;
        static const std::string POT_FEHE;
        static const std::string POT_LENOSKY_SI;
        static const std::string POT_SW_SI;
        static const std::string POT_TERSOFF_SI;
#endif

#ifndef WIN32
#ifdef WITH_VASP
        static const std::string POT_VASP;
#endif
#endif

        // ??? Unused
        static const std::string POT_BOPFOX;
        static const std::string POT_BOP;

#ifdef LAMMPS_POT
        static const std::string POT_LAMMPS;
#endif

#ifdef EONMPI
        static const std::string POT_MPI;
#endif

#ifdef WITH_PYTHON
 #ifdef PYAMFF_POT
        static const std::string POT_PYAMFF;
 #endif
        static const std::string POT_QSC;
#endif

#ifdef WITH_AMS
        static const std::string POT_AMS;
        static const std::string POT_AMS_IO;
#endif

        static Potential* getPotential(Parameters *parameters);

        static int fcalls;
        static int fcallsTotal;
        static double totalUserTime;
        
        Parameters *params;

        AtomMatrix force(long nAtoms, AtomMatrix positions,
                         VectorXi atomicNrs, double *energy, Matrix3d box, int nImages);

        void virtual initialize() = 0;
        void virtual force(long nAtoms, const double *positions,
                           const int *atomicNrs, double *forces, double *energy,
                           const double *box, int nImages) = 0;

        static Potential* pot;

};

#endif
