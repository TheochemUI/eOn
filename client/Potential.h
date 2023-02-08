#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Parameters.h"
#include "Eigen.h"

class Potential
{

    public:

        Potential(){}
        virtual ~Potential(){}

        static const char POT_EMT[];
        static const char POT_EXT[];
        static const char POT_LJ[];
        static const char POT_LJCLUSTER[];
        static const char POT_MORSE_PT[];
        static const char POT_NEW[];

#ifdef CUH2_POT
        static const char POT_CUH2[];
#endif

#ifdef IMD_POT
        static const char POT_IMD[];
#endif

#ifdef WITH_GPRD
        static const char POT_GPR[];
#endif

#ifdef WITH_WATER
        static const char POT_TIP4P[];
        static const char POT_TIP4P_PT[];
 #ifdef WITH_FORTRAN
        static const char POT_TIP4P_H[];
 #endif
        static const char POT_SPCE[];
#endif

#ifdef WITH_FORTRAN
        static const char POT_EAM_AL[];
        static const char POT_EDIP[];
        static const char POT_FEHE[];
        static const char POT_LENOSKY_SI[];
        static const char POT_SW_SI[];
        static const char POT_TERSOFF_SI[];
#endif

#ifndef WIN32
#ifdef WITH_VASP
        static const char POT_VASP[];
#endif
#endif

        // ??? Unused
        static const char POT_BOPFOX[];
        static const char POT_BOP[];

#ifdef LAMMPS_POT
        static const char POT_LAMMPS[];
#endif

#ifdef EONMPI
        static const char POT_MPI[];
#endif

#ifdef WITH_PYTHON
 #ifdef PYAMFF_POT
        static const char POT_PYAMFF[];
 #endif
        static const char POT_QSC[];
#endif

#ifdef WITH_AMS
        static const char POT_AMS[];
        static const char POT_AMS_IO[];
#endif

        static Potential* getPotential(Parameters *parameters);

        static int fcalls;
        static int fcallsTotal;
        static int wu_fcallsTotal;
        static double totalUserTime;
        
        Parameters *params;

        AtomMatrix force(long nAtoms, AtomMatrix positions,
                         VectorXi atomicNrs, double *energy, Matrix3d box);

        void virtual initialize() = 0;
        void virtual force(long nAtoms, const double *positions,
                           const int *atomicNrs, double *forces, double *energy,
                           const double *box) = 0;

        static Potential* pot;

        std::string getName() const;

};

#endif
