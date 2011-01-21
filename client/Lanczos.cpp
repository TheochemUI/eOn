//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

// Lanczos method taken from Andri Arnaldsson implementation in the vtstcode
// using the structure of the ImprovedDimer code

#include "Lanczos.h"

using namespace helper_functions;

Lanczos::Lanczos(Matter const *matter, Parameters *params)
{
    parameters    = params;
    x0            = new Matter(parameters);
    x1            = new Matter(parameters);
    *x0           = *matter;
    *x1           = *matter;

//    tau.setZero(matter->numberOfAtoms(), 3);
    w.setZero(matter->numberOfAtoms(), 3);
    q.setZero(matter->numberOfAtoms(), 3);
    qold.setZero(matter->numberOfAtoms(), 3);
    z.setZero(matter->numberOfAtoms(), 3);

    d.setZero(parameters->lanczosMaxIterations);
    e.setZero(parameters->lanczosMaxIterations);
    a.setZero(parameters->lanczosMaxIterations);
    b.setZero(parameters->lanczosMaxIterations);

    PP = new Matrix<double, Eigen::Dynamic, 3>[parameters->lanczosMaxIterations];
    for(long i=0; i<parameters->lanczosMaxIterations; i++)
    {
        PP[i].setZero(matter->numberOfAtoms(), 3);
    }

    totalForceCalls = 0;
}

Lanczos::~Lanczos()
{
    delete x0;
    delete x1;
    delete [] PP;
}

void Lanczos::initialize(Matter const *matter, Matrix<double, Eigen::Dynamic, 3> displacement)
{
    *x0 = *matter;
    *x1 = *matter;
    tau = displacement.cwise() * matter->getFree();
    tau.normalize();
    
    Matrix<double, Eigen::Dynamic, 3> x0_r = x0->getPositions();
    x1->setPositions(x0_r + tau * parameters->lanczosSeparation);
}

void Lanczos::compute(Matter const *matter)
{
    
    *x0 = *matter;
    *x1 = *matter;
    Matrix<double, Eigen::Dynamic, 3> x0_r = x0->getPositions();
    Matrix<double, Eigen::Dynamic, 3> f0 = x0->getForces();
    Matrix<double, Eigen::Dynamic, 3> f1;

    double dR = parameters->lanczosSeparation;
    long i = 0;
    q.setZero(matter->numberOfAtoms(), 3);
    beta = 1;

    do // while we have not reached tolerance or maximum iterations
    {
        qold = q;
        q = w / beta;
        P[i] = q;
        
        x1->setPositions(x0_r + q * parameters->lanczosSeparation);
        f1 = x1->getForces();

        z = f1 - f0;
        w = z - beta * qold;
        alpha = (w.cwise() * q).sum();
        d(i) = alpha;
        w = w - alpha * q;
        beta = w.norm();
        e(i) = beta;
        // Check the eigenvalues
        if(i > 1){
            a = -d / dR;
            b = -e / dR;
            // Replace with a proper tridiagonalizer such as: 
            // CALL dsterf(it,aa(1:it),bb(1:it),info)
            fullMatrix.setZero(iteration,iteration);
            for(long j = 0, j < i, j++){
                fullMatrix(j,j) = a(j);
                if(j>0){
                    fullMatrix(j,j-1) = b(j-1);
                    fullMatrix(j-1,j) = b(j-1);
                }
            }
            Eigen::SelfAdjointEigenSolver<MatrixXd> es(fullMatrix);
            VectorXd eigval = es.eigenvalues();
            eigvalMin = eigval(1);
            eigvalChange = abs((eigvalMin-eigvalMinOld)/eigvalMinOld);
        }

        if(i > 0){
            eigvalMinOld = eigvalMin;
        }else{
            eigValOld = - alpha / dR;
        }
        #ifndef NDEBUG
//            printf("LANCZOS   -----   ---------  % 9.3e   ---------  % 9.3e  % 9.3e  %9ld   ---------\n",
//            F_R.norm(), C_tau, phi_min*(180.0/M_PI), statsRotations);
        #endif
        
    } while(eigvalChange > parameters->lanczosTolerance and iteration < parameters->lanczosMaxIterations);

    Eigen::SelfAdjointEigenSolver<MatrixXd> es(fullMatrix);
    Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigvec = es.eigenvectors();
    tau = eigvec(1);
}

double Lanczos::getEigenvalue()
{
    return C_tau;
}

void Lanczos::setEigenvector(Matrix<double, Eigen::Dynamic, 3> const eigenvector)
{
    tau   = eigenvector;
    C_tau = 0.0;
}

Matrix<double, Eigen::Dynamic, 3> Lanczos::getEigenvector()
{
    return tau;
}

