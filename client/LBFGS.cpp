//Based on the LBFGS minimizer written in ASE.

#include "LBFGS.h"
#include "Log.h"
#include <cassert>
#include <cmath>

LBFGS::LBFGS(ObjectiveFunction *objfPassed, Parameters *parametersPassed)
{
    objf = objfPassed;
    parameters = parametersPassed;

    iteration = 0;

    //Shouldn't have a memory longer than the number of degrees of freedom.
    memory = min(objf->degreesOfFreedom(), (int)parameters->optLBFGSMemory);
}

LBFGS::~LBFGS()
{
    return;
}

VectorXd LBFGS::getStep(double maxMove, VectorXd f)
{
    double H0 = parameters->optLBFGSInverseCurvature;
    VectorXd r = objf->getPositions();

    if (iteration > 0) {
        VectorXd dr = objf->difference(r,rPrev);
        //double C = dr.dot(fPrev-f)/dr.dot(dr);
        double C = (fPrev-f).dot(fPrev-f)/dr.dot(fPrev-f);
        if (C<0) {
//            log_file("[LBFGS] Negative curvature: %.4f eV/A^2 take max move step\n",C);
            reset();
            return helper_functions::maxAtomMotionAppliedV(1000*f, maxMove);
        }

        if (parameters->optLBFGSAutoScale) {
            H0 = 1./C;
//            log_file("[LBFGS] Curvature: %.4e eV/A^2\n", C); 
        }
    }

    if (iteration == 0 && parameters->optLBFGSAutoScale) {
        objf->setPositions(r+parameters->finiteDifference*f.normalized());
        VectorXd dg = objf->getGradient(true)+f;
        double C = dg.dot(f.normalized())/parameters->finiteDifference;
        H0 = 1.0/C;
        objf->setPositions(r);
        if (H0 < 0) {
//            log_file("[LBFGS] Negative curvature calculated via FD: %.4e eV/A^2, take max move step\n", C); 
            reset();
            return helper_functions::maxAtomMotionAppliedV(1000*f, maxMove);
        }else{
//            log_file("[LBFGS] Curvature calculated via FD: %.4e eV/A^2\n", C); 
        }
    }

    int loopmax = s.size();
    double a[loopmax];

    VectorXd q = -f;

    for (int i=loopmax-1;i>=0;i--) {
        a[i] = rho[i] * s[i].dot(q);
        q -= a[i] * y[i];
    }

    VectorXd z = H0 * q;

    for (int i=0;i<loopmax;i++) {
        double b = rho[i] * y[i].dot(z);
        z += s[i] * (a[i] - b);
    }

    VectorXd d = -z;

    double distance = helper_functions::maxAtomMotionV(d);
    if (distance >= maxMove && parameters->optLBFGSDistanceReset) {
//        log_file("[LBFGS] reset memory, proposed step too large: %.4f\n", distance);
        reset();
        return helper_functions::maxAtomMotionAppliedV(H0*f, maxMove);
    }

    double vd = d.normalized().dot(f.normalized());
    if (vd>1.0) vd=1.0;
    if (vd<-1.0) vd=-1.0;
    double angle = acos(vd) * (180.0 / M_PI);
    if (angle > 90.0 && parameters->optLBFGSAngleReset) {
//        log_file("[LBFGS] reset memory, angle between LBFGS angle and force too large: %.4f\n", angle);
        reset();
        return helper_functions::maxAtomMotionAppliedV(H0*f, maxMove);
    }

    return helper_functions::maxAtomMotionAppliedV(d, maxMove);
}

void LBFGS::reset(void)
{
    s.clear();
    y.clear();
    rho.clear();
}

int LBFGS::update(VectorXd r1, VectorXd r0, VectorXd f1, VectorXd f0)
{
    VectorXd s0 = objf->difference(r1, r0);

    //y0 is the change in the gradient, not the force
    VectorXd y0 = f0 - f1;

    // GH: added to prevent crashing
    if (abs(s0.dot(y0)) < LBFGS_EPS) {
        cout <<"Error in LBFGS\n";
        log_file("[LBFGS] error, s0.y0 is too small: %.4f\n", s0.dot(y0));
        return -1;
    }

    s.push_back(s0);
    y.push_back(y0);
    rho.push_back(1.0/(s0.dot(y0)));

    if ((int)s.size() > memory) {
        s.erase(s.begin());
        y.erase(y.begin());
        rho.erase(rho.begin());
    }
    return 0;
}

int LBFGS::step(double maxMove)
{
    int status = 0;
    VectorXd r = objf->getPositions();
    VectorXd f = -objf->getGradient();

    if (iteration > 0) {
        status = update(r, rPrev, f, fPrev);
    }
    if(status < 0) return -1;

    VectorXd dr = getStep(maxMove,f);

    objf->setPositions(r+dr);

    rPrev = r;
    fPrev = f;

    iteration++;

//    return objf->isConverged();
    if(objf->isConverged()) 
        return 1;
    return 0;
}


int LBFGS::run(int maxSteps, double maxMove)
{
    int status;
    while(!objf->isConverged() && iteration < maxSteps) {
        status = step(maxMove);
        if(status < 0) return -1;
    }
    if(objf->isConverged()) return 1;
    return 0;

//    return objf->isConverged();
}
