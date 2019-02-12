#include "NudgedElasticBand.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Log.h"

using namespace helper_functions;

class NEBObjectiveFunction : public ObjectiveFunction
{
    public:

        NEBObjectiveFunction(NudgedElasticBand *nebPassed, Parameters *parametersPassed)
        {
            neb = nebPassed;
            parameters = parametersPassed;
        }

        ~NEBObjectiveFunction(void){};

        VectorXd getGradient(bool fdstep=false)
        {
            VectorXd forceV;
            forceV.resize(3*neb->atoms*neb->images);
            if(neb->movedAfterForceCall) neb->updateForces();
            for(long i=1; i<=neb->images; i++){
                forceV.segment(3*neb->atoms*(i-1),3*neb->atoms) = VectorXd::Map(neb->projectedForce[i]->data(), 3*neb->atoms);
            }
            return -forceV;
        }

        double getEnergy()
        {
            double Energy=0;
            for(long i=1; i<=neb->images; i++) {
                Energy += neb->image[i]->getPotentialEnergy();
            }
            return Energy;
        }

        void setPositions(VectorXd x)
        {
            neb->movedAfterForceCall = true;
            for(long i=1; i<=neb->images; i++) {
                neb->image[i]->setPositions(MatrixXd::Map(x.segment(3*neb->atoms*(i-1),3*neb->atoms).data(),neb->atoms,3));
            }
        }

        VectorXd getPositions()
        {
            VectorXd posV;
            posV.resize(3*neb->atoms*neb->images);
            for(long i=1; i<=neb->images; i++){
                posV.segment(3*neb->atoms*(i-1),3*neb->atoms) = VectorXd::Map(neb->image[i]->getPositions().data(), 3*neb->atoms);
            }
            return posV;
        }

        int degreesOfFreedom() { return 3*neb->images*neb->atoms; }

        bool isConverged() { return getConvergence() < parameters->nebConvergedForce; }

        double getConvergence() { return neb->convergenceForce(); }

        VectorXd difference(VectorXd a, VectorXd b) {
            VectorXd pbcDiff(3*neb->images*neb->atoms);
            for (int i=1;i<=neb->images;i++) {
                int n = (i-1)*3*neb->atoms;
                int m = 3*neb->atoms;
                pbcDiff.segment(n,m) = neb->image[i]->pbcV(a.segment(n,m)-b.segment(n,m));
            }
            return pbcDiff;
        }

    private:

        NudgedElasticBand *neb;
        Parameters *parameters;
};


NudgedElasticBand::NudgedElasticBand(Matter *initialPassed, Matter *finalPassed, Parameters *parametersPassed)
{
    parameters = parametersPassed;
    images = parameters->nebImages;
    atoms = initialPassed->numberOfAtoms();
    image = new Matter *[images+2];
    tangent = new AtomMatrix *[images+2];
    projectedForce = new AtomMatrix *[images+2];
    extremumPosition = new double[2*(images+1)];
    extremumEnergy = new double[2*(images+1)];
    extremumCurvature = new double[2*(images+1)];
    numExtrema = 0;

    log("\nNEB: initialize\n");
    for(long i=0; i<=images+1; i++)
    {
        image[i] = new Matter(parameters);
        *image[i] = *initialPassed;
        tangent[i] = new AtomMatrix;
        tangent[i]->resize(atoms,3);
        projectedForce[i] = new AtomMatrix;
        projectedForce[i]->resize(atoms,3);
    }
    *image[images+1] = *finalPassed;  // final image

    AtomMatrix posInitial = image[0]->getPositions();
    AtomMatrix posFinal = image[images+1]->getPositions();
    AtomMatrix imageSep = image[0]->pbc(posFinal-posInitial)/(images+1);
    for(long i=1; i<=images; i++) {
        image[i]->setPositions(posInitial+imageSep*double(i));
    }

    movedAfterForceCall = true;

    // Make sure that the endpoints know their energy
    image[0]->getPotentialEnergy();
    image[images+1]->getPotentialEnergy();
    climbingImage = 0;

    return;
}

NudgedElasticBand::~NudgedElasticBand()
{
    clean();
    return;
}

void NudgedElasticBand::clean(void)
{
    for(long i=0; i<=images+1; i++) {
        delete image[i];
        delete tangent[i];
        delete projectedForce[i];
    }
    delete [] image;
    delete [] tangent;
    delete [] projectedForce;
    delete [] extremumPosition;
    delete [] extremumEnergy;
    delete [] extremumCurvature;
    return;
}

int NudgedElasticBand::compute(void)
{
    int status = 0;
    long iteration = 0;

    log("Nudged elastic band calculation started.\n");

    updateForces();

    NEBObjectiveFunction objf(this, parameters);

    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);

    const char *forceLabel = parameters->optConvergenceMetricLabel.c_str(); 
    log("%10s %12s %14s %11s %12s\n", "iteration", "step size", forceLabel, "max image", "max energy");
    log("---------------------------------------------------------------\n");

    char fmt[] = "%10li %12.4e %14.4e %11li %12.4f\n";
    char fmtTiny[] = "%10li %12.4e %14.4e %11li %12.4e\n";

    while (!objf.isConverged())
    {
        if (parameters->writeMovies) {
            bool append = true;
            if (iteration == 0) append = false;
            image[maxEnergyImage]->matter2con("neb_maximage.con", append);
        }
        VectorXd pos = objf.getPositions();
        if(iteration) { // so that we print forces before taking an optimizer step
            if (iteration >= parameters->nebMaxIterations) {
                status = STATUS_BAD_MAX_ITERATIONS;
                break;
            }
            optimizer->step(parameters->optMaxMove);
        }
        iteration++;

        double dE = image[maxEnergyImage]->getPotentialEnergy() - 
                    image[0]->getPotentialEnergy();
        double stepSize = helper_functions::maxAtomMotionV(image[0]->pbcV(objf.getPositions()-pos));
        if (dE > 0.01) {
            log(fmt, iteration, stepSize, convergenceForce(), maxEnergyImage, dE);
        }else{
            log(fmtTiny, iteration, stepSize, convergenceForce(), maxEnergyImage, dE);
        }
    }

    if(objf.isConverged()) {
        status = STATUS_GOOD;
        log("NEB converged\n");
    }

    printImageData();
    findExtrema();

    delete optimizer;
    return status;
}


// generate the force value that is compared to the convergence criterion
double NudgedElasticBand::convergenceForce(void)
{
    if(movedAfterForceCall) updateForces();
    double fmax = 0;

    for(long i=1; i<=images; i++) {
    
            if( parameters->nebClimbingImageConvergedOnly == true && 
                parameters->nebClimbingImageMethod && 
                climbingImage!=0 ) {
                i = climbingImage;
            }
            if (parameters->optConvergenceMetric == "norm") {
                fmax = max(fmax, projectedForce[i]->norm());
            } else if (parameters->optConvergenceMetric == "max_atom") {
                for (int j=0;j<image[0]->numberOfAtoms();j++) {
                    if (image[0]->getFixed(j)) continue;
                    fmax = max(fmax, projectedForce[i]->row(j).norm());
                }
            } else if (parameters->optConvergenceMetric == "max_component") {
                fmax = max(fmax, projectedForce[i]->maxCoeff());
            } else {
                log("[Nudged Elastic Band] unknown opt_convergence_metric: %s\n", parameters->optConvergenceMetric.c_str());
                exit(1);
            }
            if( parameters->nebClimbingImageConvergedOnly == true && 
                parameters->nebClimbingImageMethod && 
                climbingImage!=0 ) {
                break;
            }
    }
    return fmax;
}

// Update the forces, do the projections, and add spring forces
void NudgedElasticBand::updateForces(void)
{
    // variables for tangent
    double maxDiffEnergy, minDiffEnergy;
    double energyDiffPrev, energyDiffNext;
    double energy, energyPrev, energyNext;
    //bool higherEnergyPrev, higherEnergyNext;

    // variables for climbing image
    double maxEnergy;

    // variables for force projections
    AtomMatrix force(atoms,3), forcePerp(atoms,3), forcePar(atoms,3);
    AtomMatrix forceSpringPar(atoms,3), forceSpring(atoms,3), forceSpringPerp(atoms,3);
    AtomMatrix forceDNEB(atoms, 3);
    AtomMatrix pos(atoms,3), posNext(atoms,3), posPrev(atoms,3);
    double distNext, distPrev;

    // update the forces on the images and find the highest energy image
    maxEnergy = image[0]->getPotentialEnergy();
    maxEnergyImage = 0;
    for(long i=1; i<=images+1; i++) {
        image[i]->getForces();
        if(image[i]->getPotentialEnergy() > maxEnergy) {
            maxEnergy = image[i]->getPotentialEnergy();
            maxEnergyImage = i;
        }
    }

    for(long i=1; i<=images; i++)
    {
        // set local variables
        force = image[i]->getForces();
        pos = image[i]->getPositions();
        posPrev = image[i-1]->getPositions();
        posNext = image[i+1]->getPositions();
        energy = image[i]->getPotentialEnergy();
        energyPrev = image[i-1]->getPotentialEnergy();
        energyNext = image[i+1]->getPotentialEnergy();

        // determine the tangent
        if(parameters->nebOldTangent) {
            // old tangent
            *tangent[i] = image[i]->pbc(posNext - posPrev);
        }else{
          // new improved tangent
	  //higherEnergyPrev = energyPrev > energyNext;
	  //higherEnergyNext = energyNext > energyPrev;

	  if(energyNext > energy && energy > energyPrev) {
	    *tangent[i] = image[i]->pbc(posNext - pos);
	  }else if(energy > energyNext && energyPrev > energy){
	    *tangent[i] = image[i]->pbc(pos - posPrev);
	  }else{
	    // we are at an extremum
	    energyDiffPrev = energyPrev - energy;
	    energyDiffNext = energyNext - energy;
	    
	    // calculate the energy difference to neighboring images
	    minDiffEnergy = min(abs(energyDiffPrev), abs(energyDiffNext));
	    maxDiffEnergy = max(abs(energyDiffPrev), abs(energyDiffNext));
	    
	    // use these energy differences to weight the tangent
	    if(energyDiffPrev > energyDiffNext) {
	      *tangent[i] = image[i]->pbc(posNext - pos) * minDiffEnergy;
	      *tangent[i] += image[i]->pbc(pos - posPrev) * maxDiffEnergy;
	    }else{
	      *tangent[i] = image[i]->pbc(posNext - pos) * maxDiffEnergy;
	      *tangent[i] += image[i]->pbc(pos - posPrev) * minDiffEnergy;
	    }
	  }
	}
        tangent[i]->normalize();

        // project the forces and add springs
        force = image[i]->getForces();

        // calculate the force perpendicular to the tangent
        forcePerp = force - (force.cwise() * *tangent[i]).sum() * *tangent[i];
        forceSpring = parameters->nebSpring * image[i]->pbc((posNext - pos) - (pos - posPrev));

        // calculate the spring force
        distPrev = image[i]->pbc(posPrev - pos).squaredNorm();
        distNext = image[i]->pbc(posNext - pos).squaredNorm();
        forceSpringPar = parameters->nebSpring * (distNext-distPrev) * *tangent[i];

        if (parameters->nebDoublyNudged) {
            forceSpringPerp = forceSpring - (forceSpring.cwise() * *tangent[i]).sum() * *tangent[i];
            forceDNEB = forceSpringPerp - (forceSpringPerp.cwise() * forcePerp.normalized()).sum() * forcePerp.normalized();
            if (parameters->nebDoublyNudgedSwitching) {
                double switching;
                switching = 2.0/M_PI * atan(pow(forcePerp.norm(),2) / pow(forceSpringPerp.norm(),2)); 
                forceDNEB *= switching;
            }
        }else{
            forceDNEB.setZero();
        }


        if(parameters->nebClimbingImageMethod && i==maxEnergyImage)
        {
            // we are at the climbing image
            climbingImage = maxEnergyImage;
            *projectedForce[i] = force - (2.0 * (force.cwise() * *tangent[i]).sum() * *tangent[i]) + forceDNEB;
        }
        else  // all non-climbing images
        {
            // sum the spring and potential forces for the neb force
            if (parameters->nebElasticBand) {
                *projectedForce[i] = forceSpring + force;
            }else{
                *projectedForce[i] = forceSpringPar + forcePerp + forceDNEB;
            }
            //*projectedForce[i] = forceSpring + forcePerp;


            //if (parameters->nebFullSpring) {


            movedAfterForceCall = false;  // so that we don't repeat a force call
        }

        //zero net translational force
        if (image[i]->numberOfFreeAtoms() == image[i]->numberOfAtoms()) {
            for (int j=0;j<=2;j++) {
                double translationMag = projectedForce[i]->col(j).sum();
                int natoms = projectedForce[i]->col(j).size();
                projectedForce[i]->col(j).array() -= translationMag/((double)natoms);
            }
        }
    }

    return;
}

// Print NEB image data
void NudgedElasticBand::printImageData(bool writeToFile)
{
    double dist, distTotal=0;
    AtomMatrix tangentStart = image[0]->pbc(image[1]->getPositions() - image[0]->getPositions());
    AtomMatrix tangentEnd = image[images]->pbc(image[images+1]->getPositions() - image[images]->getPositions());
    AtomMatrix tang;

    log("Image data (as in neb.dat)\n");

    FILE *fh=NULL;
    if (writeToFile) {
        fh = fopen("neb.dat", "w");
    }

    for(long i=0; i<=images+1; i++)
    {
        if(i==0){ tang = tangentStart; }
        else if (i==images+1) { tang = tangentEnd; }
        else { tang = *tangent[i]; }
        if(i>0) {
            dist = image[i]->distanceTo(*image[i-1]);
            distTotal += dist;
        }
        if (fh == NULL) {
            log("%3li %12.6f %12.6f %12.6f\n",i,distTotal,
                image[i]->getPotentialEnergy()-image[0]->getPotentialEnergy(), (image[i]->getForces().cwise()*tang).sum());
        }else{
            fprintf(fh, "%3li %12.6f %12.6f %12.6f\n",i,distTotal,
                image[i]->getPotentialEnergy()-image[0]->getPotentialEnergy(), (image[i]->getForces().cwise()*tang).sum());
        }
    }
    if (writeToFile) {
        fclose(fh);
    }
}

// Estimate the barrier using a cubic spline
void NudgedElasticBand::findExtrema(void)
{
    // calculate the cubic parameters for each interval (a,b,c,d)

    AtomMatrix tangentEndpoint;
    double a[images+1], b[images+1], c[images+1], d[images+1];
    double F1, F2, U1, U2, dist;

    for(long i=0; i<=images; i++)
    {
        dist = image[i]->distanceTo(*image[i+1]);
        if(i==0) {
            tangentEndpoint = image[i]->pbc(image[1]->getPositions() - image[0]->getPositions());
            tangentEndpoint.normalize();
            F1 = (image[i]->getForces().cwise()*tangentEndpoint).sum()*dist;
        } else {
            F1 = (image[i]->getForces().cwise()*(*tangent[i])).sum()*dist;
        }
        if(i==images) {
            tangentEndpoint =  image[i+1]->pbc(image[images+1]->getPositions() - image[images]->getPositions());
            tangentEndpoint.normalize();
            F2 = (image[i+1]->getForces().cwise()*tangentEndpoint).sum()*dist;
        } else {
            F2 = (image[i+1]->getForces().cwise()*(*tangent[i+1])).sum()*dist;
        }
        U1 = image[i]->getPotentialEnergy();
        U2 = image[i+1]->getPotentialEnergy();
        a[i] = U1;
        b[i] = -F1;
        c[i] = 3.*(U2 - U1) + 2.*F1 + F2;
        d[i] = -2.*(U2 - U1) - (F1 + F2);
    }

    // finding extrema along the MEP

//    long numExtrema = 0;
//    double extremaEnergy[2*(images+1)]; // the maximum number of extrema
//    double extremaPosition[2*(images+1)];
    double discriminant, f;

    for(long i=0; i<=images; i++)
    {
        discriminant = pow(c[i],2) - 3.*b[i]*d[i];
        if(discriminant >= 0) {
            f = -1;

            // quadratic case
            if( (d[i] == 0) && (c[i] != 0) ) {
                f = ( -b[i]/(2.*c[i]) );
            }
            // cubic case 1
            else if( d[i] != 0 ) {
                f = -(c[i] + sqrt(discriminant))/(3.*d[i]);
            }
            if( (f >= 0) && (f <= 1) ) {
                extremumPosition[numExtrema] = i + f;
                extremumEnergy[numExtrema] = d[i]*pow(f,3) + c[i]*pow(f,2) + b[i]*f + a[i];
                extremumCurvature[numExtrema] = 6.0*d[i]*f + 2*c[i];
                numExtrema ++;
            }
            // cubic case 2
            if( d[i] != 0 ) {
                f = ( -(c[i] - sqrt(discriminant))/(3.*d[i]) );
            }
            if( (f >= 0) && (f <= 1) ) {
                extremumPosition[numExtrema] = i + f;
                extremumEnergy[numExtrema] = d[i]*pow(f,3) + c[i]*pow(f,2) + b[i]*f + a[i];
                extremumCurvature[numExtrema] = 6*d[i]*f + 2*c[i];
                numExtrema ++;
            }
        }
    }

    log("\nFound %li extrema\n",numExtrema);
    log("Energy reference: %f\n",image[0]->getPotentialEnergy());
    for(long i=0; i<numExtrema; i++) {
        log("extrema #%li at image position %f with energy %f and curvature %f\n",i+1,extremumPosition[i],extremumEnergy[i]-image[0]->getPotentialEnergy(), extremumCurvature[i]);
    }
}
