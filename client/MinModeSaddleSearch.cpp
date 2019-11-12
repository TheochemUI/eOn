#include "MinModeSaddleSearch.h"
#include "ConjugateGradients.h"
#include "HelperFunctions.h"
#include "Lanczos.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "EpiCenters.h"
#include "ObjectiveFunction.h"
#include "Log.h"

#include <cstdlib>

using namespace helper_functions;

class MinModeObjectiveFunction : public ObjectiveFunction
{
    private:
        AtomMatrix eigenvector;
        LowestEigenmode *minModeMethod;
        Matter *matter;
        Parameters *parameters;
        int iteration;
    
    public:
        MinModeObjectiveFunction(Matter *matterPassed, LowestEigenmode *minModeMethodPassed,
                                 AtomMatrix modePassed, Parameters *parametersPassed)
        {
            matter = matterPassed;
            minModeMethod = minModeMethodPassed;
            eigenvector = modePassed;
            parameters = parametersPassed;
        }

        ~MinModeObjectiveFunction(void)
        { }

        VectorXd getGradient(bool fdstep=false)
        {
            AtomMatrix proj;
            AtomMatrix force = matter->getForces();

            if (!fdstep || iteration == 0) {
                minModeMethod->compute(matter, eigenvector);
                iteration++;
            }

            eigenvector = minModeMethod->getEigenvector();
            double eigenvalue = minModeMethod->getEigenvalue();

            proj = (force.cwise() * eigenvector).sum() * eigenvector.normalized();

            if (0 < eigenvalue) {
                if (parameters->saddlePerpForceRatio > 0.0) {
                    // reverse force parallel to eigenvector, and reduce perpendicular force
                    double const d = parameters->saddlePerpForceRatio;
                    force = d*force - (1.+d)*proj;

                // zero out the smallest forces to keep displacement confined
                }else if(parameters->saddleConfinePositive) {
                    if (parameters->saddleBowlBreakout){
                        AtomMatrix forceTemp = matter->getForces();
                        double *indicies_max;
                        indicies_max = new double[parameters->saddleBowlActive];

                        // determine the force for the x largest component
                        double f_max;
                        int i_max;
                        for (int j=0; j < parameters->saddleBowlActive; j++){
                            f_max = forceTemp.row(0).norm();
                            i_max = 0;
                            for (int i=0; i<matter->numberOfAtoms(); i++) {
                                if (f_max  < forceTemp.row(i).norm()){
                                    f_max = forceTemp.row(i).norm();
                                    i_max = i;
                                }
                            }
                            forceTemp[3*i_max+0] = 0;
                            forceTemp[3*i_max+1] = 0;
                            forceTemp[3*i_max+2] = 0;
                            indicies_max[j] = i_max;
                        }
                        for (int i=0; i<matter->numberOfAtoms(); i++) {
                            forceTemp[3*i+0] = 0.0;
                            forceTemp[3*i+1] = 0.0;
                            forceTemp[3*i+2] = 0.0;
                        }
                        // only set the projected forces corresponding to the atoms subject to the largest forces
                        for (int j=0; j< parameters->saddleBowlActive; j++) {
                            forceTemp[3*indicies_max[j]+0] = -proj[3*indicies_max[j]+0];
                            forceTemp[3*indicies_max[j]+1] = -proj[3*indicies_max[j]+1];
                            forceTemp[3*indicies_max[j]+2] = -proj[3*indicies_max[j]+2];
                        }
                        force = forceTemp;
                        delete[] indicies_max;
                    }
                    else{
                        int sufficientForce = 0;
                        double minForce = parameters->saddleConfinePositiveMinForce;
                        while(sufficientForce < parameters->saddleConfinePositiveMinActive){
                            sufficientForce = 0;
                            force = matter->getForces();
                            for (int i=0; i<3*matter->numberOfAtoms(); i++) {
                                if (fabs(force[i]) < minForce)
                                    force[i] = 0;
                                else{
                                    sufficientForce = sufficientForce + 1;
                                    force[i] = -parameters->saddleConfinePositiveBoost*proj[i];
                                }
                            }
                            minForce *= parameters->saddleConfinePositiveScaleRatio;
                        }
                    }
                }else{
                    // follow eigenmode
                    force = -proj;
                }
            }else{
                // reversing force parallel to eigenmode
                force += -2.*proj;
            }

            VectorXd forceV = VectorXd::Map(force.data(), 3*matter->numberOfAtoms());
            return -forceV;
        }
        double getEnergy() { return matter->getPotentialEnergy(); }
        void setPositions(VectorXd x) { matter->setPositionsV(x); }
        VectorXd getPositions() { return matter->getPositionsV(); }
        int degreesOfFreedom() { return 3*matter->numberOfAtoms(); }
        bool isConverged() { return getConvergence() < parameters->saddleConvergedForce; }

        double getConvergence() {
            if (parameters->optConvergenceMetric == "norm") {
                return matter->getForcesFreeV().norm(); 
            } else if (parameters->optConvergenceMetric == "max_atom") {
                return matter->maxForce(); 
            } else if (parameters->optConvergenceMetric == "max_component") {
                return matter->getForces().maxCoeff(); 
            } else {
                log("[MinModeSaddleSearch] unknown opt_convergence_metric: %s\n", parameters->optConvergenceMetric.c_str());
                exit(1);
            }
        }

        VectorXd difference(VectorXd a, VectorXd b) {
            return matter->pbcV(a-b);
        }
};

MinModeSaddleSearch::MinModeSaddleSearch(Matter *matterPassed, AtomMatrix modePassed,
                          double reactantEnergyPassed, Parameters *parametersPassed)
{
    reactantEnergy = reactantEnergyPassed;
    matter = matterPassed;
    mode = modePassed;
    parameters = parametersPassed;
    status = STATUS_GOOD;
    iteration = 0;
    
    if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
        if (parameters->dimerImproved) {
            minModeMethod = new ImprovedDimer(matter, parameters);
        }else{
            minModeMethod = new Dimer(matter, parameters);
        }
    }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
        minModeMethod = new Lanczos(matter, parameters);
    }
}

MinModeSaddleSearch::~MinModeSaddleSearch()
{
    delete minModeMethod;
}

int MinModeSaddleSearch::run()
{
    log("Saddle point search started from reactant with energy %f eV.\n", reactantEnergy);

    int optStatus;
    int firstIteration = 1;

    const char *forceLabel = parameters->optConvergenceMetricLabel.c_str(); 
    if(parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
        log("[Dimer]  %9s   %9s   %10s   %18s   %9s   %7s   %6s   %4s\n", 
            "Step", "Step Size", "Delta E", forceLabel, "Curvature", 
            "Torque", "Angle", "Rots");
    }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
        log("[Lanczos]  %9s %9s %10s %18s %9s %10s %7s %5s\n", 
            "Step", "Step Size", "Delta E", forceLabel, "Curvature", "Rel Change", "Angle", "Iters");
    }

    ostringstream climb;
    climb << "climb";
    if(parameters->writeMovies)
    {
        matter->matter2con(climb.str(), false);
    }
    
    AtomMatrix initialPosition = matter->getPositions();

    MinModeObjectiveFunction objf(matter, minModeMethod, mode, parameters);
    //objf.getGradient();
    if (parameters->saddleNonnegativeDisplacementAbort) {
        objf.getGradient();
        if (minModeMethod->getEigenvalue() > 0) {
            printf("%f\n", minModeMethod->getEigenvalue());
            return STATUS_NONNEGATIVE_ABORT;
        }
    }

    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);

    while (!objf.isConverged() || iteration == 0) {

        if(!firstIteration){

            // Abort if negative mode becomes delocalized
            if(parameters->saddleNonlocalCountAbort != 0) {
                long nm = numAtomsMoved(initialPosition - matter->getPositions(), 
                                    parameters->saddleNonlocalDistanceAbort);
                if(nm >= parameters->saddleNonlocalCountAbort) {
                    status = STATUS_NONLOCAL_ABORT;
                    break;
                }
            }

            // Abort if negative mode becomes zero
            //cout << "curvature: "<<fabs(minModeMethod->getEigenvalue())<<"\n";
            if(fabs(minModeMethod->getEigenvalue()) < parameters->saddleZeroModeAbortCurvature){
                printf("%f\n", minModeMethod->getEigenvalue());
                status = STATUS_ZEROMODE_ABORT;
                break;
            }
        }
        firstIteration = 0;
        
        if (iteration >= parameters->saddleMaxIterations) {
            status = STATUS_BAD_MAX_ITERATIONS;
            break;
        }

        AtomMatrix pos = matter->getPositions();
        
        if (parameters->saddleBowlBreakout){
            // use negative step to communicate that the system is the negative region and a max step should be performed
          if ((minModeMethod->getEigenvalue() > 0) and (parameters->optMethod == "cg")){
              optStatus = optimizer->step(-parameters->optMaxMove);
          }
          else{
              optStatus = optimizer->step(parameters->optMaxMove);
          }
        }
        else{
            optStatus = optimizer->step(parameters->optMaxMove);
        }

        if(optStatus < 0){
            status = STATUS_OPTIMIZER_ERROR;
            break;
        }

        double de = objf.getEnergy()-reactantEnergy;
        // should be the total displacement of the system not just a single atom
        //double stepSize = helper_functions::maxAtomMotion(matter->pbc(matter->getPositions() - pos));
        
        double stepSize;
        
        //Melander, Laasonen, Jonsson, JCTC, 11(3), 1055–1062, 2015 http://doi.org/10.1021/ct501155k
        if(parameters->saddleRemoveRotation) {
            rotationRemove(pos, matter);
        }
        stepSize = (matter->pbc(matter->getPositions() - pos)).norm();
  
        iteration++;

        if(parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER)
        {
            log("[Dimer]  %9ld   %9.7f   %10.4f   %18.5e   %9.4f   %7.3f   %6.3f   %4ld\n",
                        iteration, stepSize, matter->getPotentialEnergy()-reactantEnergy,
                        objf.getConvergence(),
                        minModeMethod->getEigenvalue(),
                        minModeMethod->statsTorque,
                        minModeMethod->statsAngle,
                        minModeMethod->statsRotations);
        }else if (parameters->saddleMinmodeMethod == LowestEigenmode::MINMODE_LANCZOS) {
            log("[Lanczos]  %9i %9.6f %10.4f %18.5e %9.4f %10.6f %7.3f %5i\n", 
                iteration, stepSize, matter->getPotentialEnergy()-reactantEnergy,
                objf.getConvergence(),
                minModeMethod->getEigenvalue(),
                minModeMethod->statsTorque,
                minModeMethod->statsAngle,
                minModeMethod->statsRotations);
        }else{
            log("[MinModeSaddleSearch] Unknown min_mode_method: %s\n", parameters->saddleMinmodeMethod.c_str());
            exit(1);
        }

        if (parameters->writeMovies) {
            matter->matter2con(climb.str(), true);
        }

        if (de > parameters->saddleMaxEnergy) {
            status = STATUS_BAD_HIGH_ENERGY;
            break;
        }

        if (parameters->checkpoint) {
            matter->matter2con("displacement_cp.con", false);
            FILE *fileMode = fopen("mode_cp.dat", "wb");
            helper_functions::saveMode(fileMode, matter, 
                                       minModeMethod->getEigenvector());
            fclose(fileMode);
        }
    }

    if (iteration == 0) minModeMethod->compute(matter, mode);

    if (getEigenvalue() > 0.0 && status == STATUS_GOOD) {
        log("[MinModeSaddleSearch] eigenvalue not negative\n");
        status = STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE;
    }

    delete optimizer;

    return status;
}

double MinModeSaddleSearch::getEigenvalue()
{
    return minModeMethod->getEigenvalue();
}

AtomMatrix MinModeSaddleSearch::getEigenvector()
{
    return minModeMethod->getEigenvector();
}
