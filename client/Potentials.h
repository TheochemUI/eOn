/*
 *===============================================
 *  EON Potential.h
 *===============================================
 */

#include "Parameters.h"
#include "PotentialsInterface.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

enum{
    POT_USER=0,
    POT_LJ,
    POT_MORSE,
    POT_EMT,
    POT_EDIP,
    POT_VASP,
    POT_TERSOFF,
    POT_SW,
    POT_LENOSKY,
    POT_LJBINARY,
    POT_ALUMINUM,
    POT_EAM,
    POT_QSC,
    POT_ZPICE,
    POT_TIP4P,
    POT_BOPFOX,
    N_POTS
};

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
    Matrix<double, Eigen::Dynamic, 3> force(long nAtoms, Matrix<double, Eigen::Dynamic, 3> positions, Matrix<int, Eigen::Dynamic, 1> atomicNrs, double *energy, Matrix<double, 3, 3> box);
 
private:
    PotentialsInterface *interface_;
    Parameters *parameters_;
};
