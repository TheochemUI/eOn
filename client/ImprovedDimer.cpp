//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

// An implementation of Johannes Kästner and Paul Sherwood's improved dimer.
// An attempt to keep to the variable names in their 2008 paper has been made.

#include "HelperFunctions.h"
#include "ImprovedDimer.h"
#include <cmath>
#include <cassert>

using namespace helper_functions;

const char ImprovedDimer::OPT_SD[] = "sd";
const char ImprovedDimer::OPT_CG[] = "cg";
const char ImprovedDimer::OPT_LBFGS[] = "lbfgs";

ImprovedDimer::ImprovedDimer(Matter const *matter, Parameters *params)
{
    parameters    = params;
    x0            = new Matter(parameters);
    x1            = new Matter(parameters);
    *x0           = *matter;
    *x1           = *matter;
    tau.resize(matter->numberOfAtoms(), 3);
    tau.setZero();
    totalForceCalls = 0;

    if(parameters->dimerOptMethod == OPT_CG){
        init_cg = true;
    }
}

ImprovedDimer::~ImprovedDimer()
{
    delete x0;
    delete x1;
}

void ImprovedDimer::compute(Matter const *matter, AtomMatrix initialDirection)
{
    tau = initialDirection.cwise() * matter->getFree();
    tau.normalize();
    *x0 = *matter;
    *x1 = *matter;
    AtomMatrix x0_r = x0->getPositions();

//GH    x1->setPositions(x0_r + tau * parameters->dimerSeparation);
    x1->setPositions(x0_r + tau * parameters->dimerSeparation * 0.5);

    // other vectors
    AtomMatrix x1_rp;
    AtomMatrix x1_r;
    AtomMatrix tau_prime;
    AtomMatrix g1_prime;

    double delta = parameters->dimerSeparation * 0.5;
    double phi_tol = M_PI * (parameters->dimerConvergedAngle/180.0);
    double phi_prime = 0.0;
    double phi_min = 0.0;

    statsRotations = 0;

    Matter *x1p = new Matter(parameters);

    // Calculate the gradients on x0 and x1, g0 and g1, respectively.
    AtomMatrix g0 = -x0->getForces();
    AtomMatrix g1 = -x1->getForces();

    do // while we have not reached phi_tol or maximum rotations.
    {

        // Calculate the rotational force, F_R.
        F_R = -2.0 * (g1 - g0) + 2.0 * ((g1 - g0).cwise() * tau).sum() * tau;

//GH        statsTorque = F_R.norm()/delta;
        statsTorque = F_R.norm()/parameters->dimerSeparation;

        // Determine the step direction, theta
        if(parameters->dimerOptMethod == OPT_SD) // steepest descent
        {
            theta = F_R / F_R.norm();
        }
        else if(parameters->dimerOptMethod == OPT_CG) // conjugate gradients
        {
            if(init_cg){
                init_cg = false;
                gamma = 0.0;
            }else{
                a = fabs((F_R.cwise() * F_R_Old).sum());
                b = F_R_Old.squaredNorm();
                if(a < 0.5*b) {
                    gamma = (F_R.cwise() * (F_R - F_R_Old)).sum()/b;
                }else{ 
                    gamma = 0.0;
                }
            }

            if(gamma == 0.0) {
                theta = F_R;
            }else{
                theta = F_R + thetaOld * F_R_Old.norm() * gamma;
            }

            theta = theta - (theta.cwise() * tau).sum() * tau;
            theta = theta / theta.norm();

            F_R_Old = F_R;
            thetaOld = theta;
        }
        else if(parameters->dimerOptMethod == OPT_LBFGS) // quasi-newton
        {
            // LBFGS not implemented yet
        }

        // Calculate the curvature along tau, C_tau.
        C_tau = ((g1 - g0).cwise() * tau).sum() / delta;

        // Calculate a rough estimate (phi_prime) of the optimum rotation angle.
        double d_C_tau_d_phi = 2.0 * ((g1 - g0).cwise() * theta).sum() / delta;
        phi_prime = -0.5 * atan(d_C_tau_d_phi / (2.0 * abs(C_tau)));
        statsAngle = phi_prime * (180.0 / M_PI);
/*
        cout <<"initial angle: "<<statsAngle;
        if(phi_prime > phi_tol){ cout <<" full rotation\n"; }
        else{ cout <<" no full rotation\n"; }
*/
        if(phi_prime > phi_tol)
        {
            double b1 = 0.5 * d_C_tau_d_phi;

            // Calculate g1_prime. 
            x0_r = x0->getPositions();
            x1_rp = x0_r + (tau * cos(phi_prime) + theta * sin(phi_prime)) * delta;
            *x1p = *x1;
            x1p->setPositions(x1_rp);
            g1_prime = -x1p->getForces();

            // Calculate C_tau_prime.
            tau_prime = (x1_rp - x0_r) / (x1_rp - x0_r).norm();
            double C_tau_prime = ((g1_prime - g0).cwise() * tau_prime).sum() / delta;

            // Calculate phi_min.
            double a1 = C_tau - C_tau_prime + b1 * sin(2.0 * phi_prime) / (1.0 - cos(2.0 * phi_prime));
            double a0 = 2.0 * (C_tau - a1);
            phi_min = 0.5 * atan(b1 / a1);

//GH: test showing that a smaller angle avoids oscillations
//            phi_min = phi_min/2.0;

            // Determine the curvature for phi_min.
            double C_tau_min = 0.5 * a0 + a1 * cos(2.0 * phi_min) + b1 * sin(2.0 * phi_min);

            // If the curvature is being maximized, push it over pi/2.
            if(C_tau_min > C_tau)
            {
//                cout<<"flip\n";
                phi_min += M_PI * 0.5;
                C_tau_min = 0.5 * a0 + a1 * cos(2.0*phi_min) + b1 * sin(2.0*phi_min);
            }

            statsAngle = phi_min * (180.0 / M_PI);

//            cout    <<"b1: "<<b1 <<" a1: "<<a1 <<" a0: "<<a0 <<" phi: "<<phi_min*180.0/M_PI<<endl;

            // Update x1, tau, and C_tau.
            x1_r = x0_r + (tau * cos(phi_min) + theta * sin(phi_min)) * delta;
            x1->setPositions(x1_r);
            tau = (x1_r - x0_r) / (x1_r - x0_r).norm();
//            cout <<tau<<endl;

            if(parameters->saddleMaxLocalizedAtoms > 0)
            {
//                cout <<"localize\n";
                tau = localize(tau, parameters->saddleMaxLocalizedAtoms);
                tau.normalize();
            }
            C_tau = C_tau_min;

            // Calculate the new g1.
            g1 = g1 * (sin(phi_prime - phi_min)/sin(phi_prime)) + g1_prime*(sin(phi_min)/sin(phi_prime)) + 
                 g0 * (1.0 - cos(phi_min) - sin(phi_min) * tan(phi_prime * 0.5));

//GH            statsTorque = F_R.norm()/delta;
            statsTorque = F_R.norm()/parameters->dimerSeparation;
            statsRotations += 1;

            FILE *fp = fopen("saddlesearch.dat", "a");
            fprintf(fp, "IDIMERROT  -----   ---------   ----------------   ---------  % 9.3e  % 9.3e  % 9.3e   % 9ld\n",
                        C_tau, statsTorque, statsAngle, statsRotations);
            printf("IDIMERROT  -----   ---------   ----------------   ---------  % 9.3e  % 9.3e  % 9.3e   % 9ld\n",
                   C_tau, statsTorque, statsAngle, statsRotations);
            fclose(fp);
        }
        else
        {
            FILE *fp = fopen("saddlesearch.dat", "a");
            fprintf(fp, "IDIMERROT  -----   ---------   ----------------   ---------  % 9.3e  % 9.3e   ---------   ---------\n",
                        C_tau, F_R.norm()/delta);
            printf("IDIMERROT  -----   ---------   ----------------   ---------  % 9.3e  % 9.3e   ---------   ---------\n",
                   C_tau, F_R.norm()/delta);
            fclose(fp);
        }

    } while(phi_prime > phi_tol and phi_min > phi_tol and statsRotations < parameters->dimerRotationsMax);

    delete x1p;

}

double ImprovedDimer::getEigenvalue()
{
    return C_tau;
}

AtomMatrix ImprovedDimer::getEigenvector()
{
    return tau;
}


