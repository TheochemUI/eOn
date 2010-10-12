/*
 *===============================================
 *  EON Potential.h
 *===============================================
 */

#include "Parameters.h"
#include "PotentialsInterface.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

#define POT_USER 0
#define POT_LJ 1
#define POT_MORSE 2
#define POT_EMT 3
#define POT_EDIP 4
#define POT_VASP 5
#define POT_TERSOFF 6
#define POT_SW 7
#define POT_LENOSKY 8
#define POT_LJBINARY 9
#define POT_ALUMINUM 10
#define POT_EAM 11
#define POT_QSC 12
#define POT_ZPICE 13
#define POT_TIP4P 14

/* Class serving as a wrapper between the force calculator and the Matter object */
class Potentials{
public:
    Potentials(Parameters *parameters); ///< Constructor
 
    ~Potentials(); ///< Destructor
 
    /* this function must be provided by each force calculator
    @param[in]    nAtoms      the number of atoms
    @param[in]    *positions  pointer to the array of 3N atoms positions
    @param[in]    *atomicNrs  pointer to the array of N atomic numbers
    @param[out]   *forces     pointer to the array of 3N forces
    @param[out]   *energy     pointer to the total energy
    @param[in]    *box        pointer to the array containing the 3 lengths of the supercell */
    void force(long nAtoms, Matrix<double, Eigen::Dynamic, 3> positions, Matrix<int, Eigen::Dynamic, 1> atomicNrs, Matrix<double, Eigen::Dynamic, 3> forces, double *energy, Matrix<double, 3, 3> box);
 
private:
    PotentialsInterface *interface_;
    Parameters *parameters_;
};
