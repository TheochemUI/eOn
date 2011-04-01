//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "ExactMinMode.h"
#include <cmath>

ExactMinMode::ExactMinMode(Matter const *matter, Parameters *params)
{
    parameters = params;
    lowestEv.resize(matter->numberOfAtoms(),3);
    lowestEv.setZero();
    lowestEw = 0.0;
}

ExactMinMode::~ExactMinMode()
{
}

void ExactMinMode::compute(Matter const *matter, AtomMatrix direction)
{
    int size = 3*matter->numberOfFreeAtoms();
    Matter *tmpMatter = new Matter(parameters);
    *tmpMatter = *matter;
    double dr = parameters->lanczosFiniteDist;
    VectorXd force1, force2;


    VectorXd posDisplace(size);
    MatrixXd hessian(size,size);
    hessian.setZero();
    force1 = tmpMatter->getForcesFreeV();
    for (int i=0;i<size;i++) {
        posDisplace = matter->getPositionsFreeV();
        posDisplace(i) += dr;
        tmpMatter->setPositionsFreeV(posDisplace);
        force2 = tmpMatter->getForcesFreeV();
        hessian.col(i) = -(force2-force1)/dr;
    }
    hessian = (hessian+hessian.transpose())/2;
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(hessian);
    lowestEw = es.eigenvalues()(0);
    VectorXd exactEv = es.eigenvectors().col(0);

    lowestEv.resize(matter->numberOfAtoms(),3);
    int j=0;
    for (int i=0;i<matter->numberOfAtoms();i++) {
        if (!matter->getFixed(i)) {
            lowestEv.row(i) = exactEv.segment<3>(j);
            j++;
        }
    }

    delete tmpMatter;
}

double ExactMinMode::getEigenvalue()
{
    return lowestEw;
}

AtomMatrix ExactMinMode::getEigenvector()
{
    return lowestEv;
}
