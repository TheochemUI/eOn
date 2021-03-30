// This Lanczos algorithm is implemented as described in this paper:
// R. A. Olsen, G. J. Kroes, G. Henkelman, A. Arnaldsson, and H. Jónsson, 
// Comparison of methods for finding saddle points without knowledge of the final states,
// J. Chem. Phys. 121, 9776-9792 (2004).

#include "Lanczos.h"
#include "Log.h"
#include <cmath>

Lanczos::Lanczos(Matter *matter, Parameters *params)
{
    parameters = params;
    lowestEv.resize(matter->numberOfAtoms(),3);
    lowestEv.setZero();
    lowestEw = 0.0;
}

Lanczos::~Lanczos()
{
}

// The 1 character variables in this method match the variables in the
// equations in the paper given at the top of this file.
void Lanczos::compute(Matter *matter, AtomMatrix direction)
{
    int size = 3*matter->numberOfFreeAtoms();
    MatrixXd T(size,parameters->lanczosMaxIterations), Q(size,parameters->lanczosMaxIterations);
    T.setZero();
    VectorXd u(size), r(size);

    // Convert the AtomMatrix of all the atoms into
    // a single column vector with just the free coordinates.
    int i,j;
    for (i=0,j=0;i<matter->numberOfAtoms();i++) {
        if (!matter->getFixed(i)) {
            r.segment<3>(j) = direction.row(i);
            j+=3;
        }
    }

    double alpha, beta=r.norm();
    double ew=0, ewOld=0, ewAbsRelErr;
    double dr = parameters->finiteDifference;
    VectorXd evEst, evT, evOldEst;

    VectorXd force1, force2;
    Matter *tmpMatter = new Matter(parameters);
    *tmpMatter = *matter;
    force1 = tmpMatter->getForcesFreeV();

    for (i=0;i<size;i++) {
        statsRotations = i;
        Q.col(i) = r/beta;

        // Finite difference force in the direction of the ith Lanczos vector
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

        T(i,i) = alpha;
        if (i>0) {
            T(i-1,i) = beta;
            T(i,i-1) = beta;
        }

        beta = r.norm();

        if (beta <= 1e-10*fabs(alpha)) {
            /* If Q(0) is an eigenvector (or a linear combination of a subset of eignevectors)
            then the lanczos cannot complete the basis of vector Q.*/
            if (i == 0) {
                ew = alpha;
                evEst = Q.col(0);
            }
            log_file("[ILanczos] ERROR: linear dependence\n");
            break;
        }
        //Check Eigenvalues
        if (i >= 1) {
            Eigen::SelfAdjointEigenSolver<MatrixXd> es(T.block(0,0,i+1,i+1));
            ew = es.eigenvalues()(0); 
            evT = es.eigenvectors().col(0);
            ewAbsRelErr = fabs((ew-ewOld)/ewOld);
            ewOld = ew;

            //Convert eigenvector of T matrix to eigenvector of full Hessian
            evEst = Q.block(0,0,size,i+1)*evT;
            evEst.normalize();
            statsAngle = acos(fabs(evEst.dot(evOldEst)))*(180/M_PI);
            statsTorque = ewAbsRelErr;
            evOldEst = evEst;
            log_file("[ILanczos] %9s %9s %10s %14s %9.4f %10.6f %7.3f %5i\n", 
                "----", "----", "----", "----", ew, ewAbsRelErr, statsAngle, i);
            if (ewAbsRelErr < parameters->lanczosTolerance) {
                log_file("[ILanczos] Tolerence reached: %f\n", parameters->lanczosTolerance);
                break;
            }
        }else{
            ew = alpha;
            ewOld = ew;
            evEst = Q.col(0);
            evOldEst = Q.col(0);
            if (lowestEw != 0.0 && parameters->lanczosQuitEarly) {
                double Cprev = lowestEw;
                double Cnew = u.dot(Q.col(i));
                ewAbsRelErr = fabs((Cnew-Cprev)/Cprev);
                if (ewAbsRelErr <= parameters->lanczosTolerance) {
                    statsAngle = 0.0;
                    statsTorque = ewAbsRelErr;
                    log_file("[ILanczos] Tolerence reached: %f\n", parameters->lanczosTolerance);
                    break;
                }
            }
        }

        if (i >= parameters->lanczosMaxIterations-1) {
            log_file("[ILanczos] Max iterations\n");
            break;
        }
    }

    lowestEw = ew;

    // Convert back from free atom coordinate column vector
    // to AtomMatrix style.
    lowestEv.resize(matter->numberOfAtoms(),3);
    for (i=0,j=0;i<matter->numberOfAtoms();i++) {
        if (!matter->getFixed(i)) {
            lowestEv.row(i) = evEst.segment<3>(j);
            j+=3;
        }
    }

    delete tmpMatter;
}

double Lanczos::getEigenvalue()
{
    return lowestEw;
}

AtomMatrix Lanczos::getEigenvector()
{
    return lowestEv;
}
