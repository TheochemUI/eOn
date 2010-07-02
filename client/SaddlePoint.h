/*
 *===============================================
 *  SaddlePoint.h
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 10/24/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *-----------------------------------------------
 *  Heavily inspired of codes by:
 *      Graeme Henkelman
 *      Roar Olsen
 *===============================================
*/
#ifndef SADDLE_POINT_H
#define SADDLE_POINT_H

class Matter;
class Parameters;
class LowestEigenmodeInterface;

/** Rely on the dimer method for saddle point determination. The object rely on an object being able to determine the lowest eigenmode and a function to determine where a displacement prior the saddle point search should be centered.*/
class SaddlePoint {
public:
    
    SaddlePoint();///< The object shall be initialized later with SaddlePoint::initialize
    
    /** Constructor where object is initialized.
    @param[in]   initial        Pointer to where the initial state (minimum) shall be stored.
    @param[in]   saddle      Pointer to where the conformation where to start the saddle point search. It will also be used to return the saddle point.
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
    SaddlePoint(Matter * initial, Matter *saddle, Parameters *parameters);

    ~SaddlePoint();///< Destructor.
    
    /** Initialized the object.
    @param[in]    *initial    Initial state (minimum)
    @param[in]   *saddle   Where to start the saddle point search (may be equal to initial).
    @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.
    */
    void initialize(Matter * initial, Matter *saddle, Parameters *parameters);

    /** Try to determine a nearby saddle point
    @param[out]  *min1     Pointer to an Matter object containing one of the minima connected to the saddle point.
    @param[out]  *min2     Pointer to an Matter object containing the other minima connected to the saddle point.
    The values returned is true if the calculation converged.*/
    long locate(Matter *min1, Matter *min2);
    LowestEigenmodeInterface const * getLowestEigenmode() const;
    long getnFreeCoord() const;
    double const *const getEigenMode() const;
private:
    Matter * initial_;
    Matter *saddle_;///< Pointer to atom object outside the scope of the class.
    Parameters *parameters_;///< Pointer to a structure outside the scope of the class containing runtime parameters. 
    LowestEigenmodeInterface *lowestEigenmode_;///< Pointer to the method used to determine the lowest eigenmode.
        
    double eigenValue_;///< Containing an estimate for the lowest eigenvalue
    double *eigenMode_;///< Double array for the lowest eigenmode, its size equals the number of free atoms times 3.  
    double *initialDisplacement_;///RT: used to keep track of the initial displacement.  
    long nFreeCoord_;///< Number of free coordinates.
    long state_;///< To keep track of where eventual problems occured.

    void clean();///< Clean up dynamical allocated memory
        
    void displaceState(Matter *matter);///< Displacing the atomic positions in /a matter according to values in Parameters, being centered on the atom determined by one of the EpiCenter functions. 
    void correctingForces(double *force);///< Modifying force. If the eigenvalue is positive (Convex) the eigenmode is followed by zeroing all other forces. If negative (Concave) the force parallel to the eigenmode is reversed. 

    /** Determine the two minima connected to the saddle point, by displacing the positions in the saddle point by either adding or subtracting a part of the lowest eigenmode
    @param[out]  *min1     Pointer to an Matter object containing one of the minima connected to the saddle point.
    @param[out]  *min2     Pointer to an Matter object containing other of the minima connected to the saddle point.*/
    void relaxFromSaddle(Matter *min1, Matter *min2);
    
    void jumpToConvexRegion();
    void displaceInConcaveRegion();
    
    void searchForSaddlePoint(double initialEnergy);
};
#endif
