/*
 *===============================================
 *  EON Prefactors.h
 *===============================================
 */
#ifndef PREFCTORS_H
#define PREFCTORS_H

#include <math.h>
#include <cassert>
#include <cstdlib>

#include "Constants.h"
#include "Matter.h"
//#include "IO.h"
#include "HelperFunctions.h"
#include "Parameters.h"

#include "Eigen/Eigen"
USING_PART_OF_NAMESPACE_EIGEN

/** Functionality to determine prefactors for the forward and reverse proces. It is assume that the potential energy surface can be expressed with harmonic functions at the vicinity of the minima and the saddle point.*/ 
class Prefactors {
public:
    Prefactors();///< The object shall be initialized later with Prefactors::initialize
    
    /** Constructor where object is initialized.
    @param[in]   *saddle       Pointer to the Matter object describing the saddle point.
    @param[in]   *min1         Pointer to the Matter object describing one of the minima being connected to the saddle point.
    @param[in]   *min2         Pointer to the Matter object describing the other minima being connected to the saddle point.
    @param[in]   *parameters   Pointer to the Parameter object containing the runtime parameters.*/    
    Prefactors(const Matter *saddle, const Matter *min1, const Matter *min2, Parameters *parameters);
    
    ~Prefactors();///< Destructor.
        
    /** Initialized the object.
    @param[in]   *saddle       Pointer to the Matter object describing the saddle point.
    @param[in]   *min1         Pointer to the Matter object describing one of the minima being connected to the saddle point.
    @param[in]   *min2         Pointer to the Matter object describing the other minima being connected to the saddle point.
    @param[in]   *parameters   Pointer to the Parameter object containing the runtime parameters.*/            
    void initialize(const Matter *saddle, const Matter *min1, const Matter *min2, Parameters *parameters);
    
    bool compute(double *prefactors);///< Determine the prefactors and store the result in \a prefactors. Value at \a prefactors[0] is for process min1->saddle, value at \a prefactors[1] is for process min2->saddle.
    long totalForceCalls;
    
private:
    long sizeHessian;///< Size of hessian. Value set by the Prefactors::atomsToAccountForInHessian
    long nAtoms;///< Number of atoms.

    bool *coordinatesToAccountFor;///< Array where index is true if there is going to be accounted for the coordinate.Value set by the Prefactors::atomsToAccountForInHessian
    VectorXd eigenValMin1;///< Eigenvalues for configuration min1_ for the degrees of freedom being accounted for.
    VectorXd eigenValMin2;///< Eigenvalues for configuration min2_ for the degrees of freedom being accounted for.
    VectorXd eigenValSaddle;///< Eigenvalues for configuration saddle_ for the degrees of freedom being accounted for.
    
    const Matter *min1;///< Pointer to object defined outside this objects scope. Should correspond to one of the minima being connected to the saddle point (saddle_).
    const Matter *min2;///< Pointer to object defined outside this objects scope. Should correspond to the other minima being connected to the saddle point (saddle_).
    const Matter *saddle;///< Pointer to object defined outside this objects scope. Should correspond to the saddle point configuration.
    
    Parameters *parameters;///< Pointer to the Parameter object containing the runtime parameters.

    void clean();///< Clean up dynamical allocated memory
        
    long atomsToAccountForInHessian();///< Determines the number of atoms to be accounted.
	long atomsMovedMoreThan(double minDisplacement);///< Determines which atoms to account for. In the analysis atoms being displaced more than minDisplacement is considered. If an atom is considered displaced the neighbors within the radius Parameters::getWithinRadiusDisplaced_Hessian are also considered.
    bool getEigenValues();///< Determine the eigenvalues. The calculation terminates if either a negative mode exists in one of the minima or if there is not only one negative mode in the saddle point. Note that there is accounted for the masses in Prefactors::massScaleHessian.
    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> determineHessian(const Matter *matter);///< Filling in the values in \a hessian for the coordinates being accounted for the configuration \a matter.
    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> massScaleHessian(Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian);///< Scale the \a hessian according to the masses of the atoms for the coordinates.
};
#endif
