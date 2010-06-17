/*
 *===============================================
 *  QMBox.h
 *  eon
 *-----------------------------------------------
 *  Created by Rye Terrell on 6/7/10.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef QMBOX_H
#define QMBOX_H

#include "MinimizersInterface.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Constants.h"
#include "Parameters.h"
#include "../common/Quickmin.h"

class QMBox : public MinimizersInterface
{

    public:
        /** Constructor to be used when a structure is minimized.
        @param[in]   *matter        Pointer to the Matter object to be relaxed.
        @param[in]   *parameters    Pointer to the Parameter object containing the runtime parameters.*/
        QMBox(Matter *matter, Parameters *parameters);

        ~QMBox();///< Destructor.

        void oneStep();///< Do one iteration.
//        void oneStepPart1(double *freeForces);
//        void oneStepPart2(double *freeForces);
        
        void fullRelax();///< Relax the Matter object corresponding to the pointer that was passed with the constructor.
        bool isItConverged(double convergeCriterion);///< Determine if the norm of the force vector is bellow the \a convergeCriterion.
        
    private:
        long nFreeCoord_;///< Number of free coordinates.
        void increment_velocity();

        Matter *matter_;///< Pointer to atom object \b outside the scope of the class.    
        Parameters *parameters_;///< Pointer to a structure outside the scope of the class containing runtime parameters. 

        double *tempListDouble_;///< Double array, its size equals the number of atoms times 3.
        double *boxforce_;
        double *boxv_;
        double dR;
        double dT;
        Quickmin *qmBox_;
        double *forces_;///< Double array, its size equals the number of \b free atoms times 3.
};

#endif
