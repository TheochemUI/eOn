//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Lanczos.h"

using namespace helper_functions;

Lanczos::Lanczos(Matter const *matter, Parameters *params)
{
    parameters = params;
    totalForceCalls = 0;
    lowestEv.resize(matter->numberOfAtoms(),3);
    lowestEv.setZero();
    lowestEw = 0.0;
}

Lanczos::~Lanczos()
{
}

void Lanczos::initialize(Matter const *matter, AtomMatrix displacement)
{
    r.resize(3*matter->numberOfFreeAtoms());
    int j=0;
    for (int i=0;i<matter->numberOfAtoms();i++) {
        if (!matter->getFixed(i)) {
            r.segment<3>(j) = displacement.row(i);
            j++;
        }
    }
    r.normalize();
    lowestEv = displacement.normalized();
}

void Lanczos::compute(Matter const *matter)
{
    int size = 3*matter->numberOfFreeAtoms();
    MatrixXd T(size,size), Q(size,size);
    T.setZero();
    VectorXd u(size);
    Matter *tmpMatter = new Matter(parameters);
    *tmpMatter = *matter;
    double dr = parameters->lanczosSeparation;
    double alpha, beta=r.norm();
    VectorXd force1, force2;
    double ew=0.0, ewOld=0.0, ewChange;
    VectorXd evT;
    force1 = tmpMatter->getForcesFreeV();

    int i=0;
    for (i=0;i<size;i++) {
        Q.col(i) = r/beta;

        // Full Gram-Schmidt orthogonalization
        // FOR TESTING ONLY
        if (i>0) {
            VectorXd h(size);
            h.setZero();
            for (int j=0;j<i;j++) {
                h += Q.col(j).dot(Q.col(i))*Q.col(j);
            }
            Q.col(i) -= h;
        }

        tmpMatter->setPositionsFreeV(matter->getPositionsFreeV()+dr*Q.col(i));
        force2 = tmpMatter->getForcesFreeV();

        u = -(force2-force1)/dr;
        if (i==0) {
            r = u;
        }else{
            r = u-beta*Q.col(i-1);
        }
        alpha = Q.col(i).dot(r);
        r = r-alpha*Q.col(i);
        beta = r.norm();

        //Add to Tridiagonal Matrix
        T(i,i) = alpha;
        if (i>0) {
            T(i-1,i) = beta;
            T(i,i-1) = beta;
        }

        //Check Eigenvalues
        if (i>1) {
            Eigen::SelfAdjointEigenSolver<MatrixXd> es(T.block(0,0,i+1,i+1));
            ew = es.eigenvalues()(0); 
            printf("lz ew: %f\n",ew);
            evT = es.eigenvectors().col(0);
            ewChange = fabs((ew-ewOld)/ewOld);
            ewOld = ew;
            if (ewChange < parameters->lanczosTolerance) {
                break;
            }
        }else{
            ewOld = -alpha/dr;
            ewChange = ewOld;
        }

        if (i >= parameters->lanczosMaxIterations) {
            printf("Hit max lanczos iterations\n");
            break;
        }
    }
    VectorXd ev = Q.block(0,0,size,i+1)*evT;
    ev.normalize();
    r=ev;
    lowestEw = ew;

    lowestEv.resize(matter->numberOfAtoms(),3);
    int j=0;
    for (i=0;i<matter->numberOfAtoms();i++) {
        if (!matter->getFixed(i)) {
            lowestEv.row(i) = ev.segment<3>(j);
            j++;
        }
    }

    printf("Building exact hessian\n");
    VectorXd posDisplace(size);
    MatrixXd hessian(size,size);
    hessian.setZero();
    for (i=0;i<size;i++) {
        posDisplace = matter->getPositionsFreeV();
        posDisplace(i) += dr;
        tmpMatter->setPositionsFreeV(posDisplace);
        force2 = tmpMatter->getForcesFreeV();
        hessian.col(i) = -(force2-force1)/dr;
    }
    hessian = (hessian+hessian.transpose())/2;
    printf("Solving eigensystem\n");
    Eigen::SelfAdjointEigenSolver<MatrixXd> esH(hessian);
    double exactEw = esH.eigenvalues()(0);
    printf("lanczosEw: %f\n", ew);
    printf("exactEw: %f\n", exactEw);
    VectorXd exactEv = esH.eigenvectors().col(0);
    printf("ew diff: %e ev dot: %e\n\n", exactEw-ew, exactEv.dot(ev));
    delete tmpMatter;
}

double Lanczos::getEigenvalue()
{
    return lowestEw;
}

void Lanczos::setEigenvector(AtomMatrix const eigenvector)
{
    printf("Not implemented\n");
}

AtomMatrix Lanczos::getEigenvector()
{
    return lowestEv;
}
