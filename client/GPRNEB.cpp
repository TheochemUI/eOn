#include "GPRNEB.h"

#include<algorithm>
#include<cmath>

class GPRNEBObjectiveFunction : public ObjectiveFunction
{
    public:

        GPRNEBObjectiveFunction(GPRNEB *gpnebPassed, Parameters *parametersPassed):
            gpneb{gpnebPassed},
            params{parametersPassed},
            threshold{parametersPassed->nebConvergedForce}
        {
        }

        ~GPRNEBObjectiveFunction(void){};

        VectorXd getGradient(bool fdstep=false)
        {
            VectorXd forceV;
            forceV.resize(3*gpneb->nfree*gpneb->nimages);
            if(gpneb->movedAfterForceCall){
                gpneb->updateForces();
            }
            for(size_t idx{1}; idx<=gpneb->nimages; idx++){
                forceV.segment(3*gpneb->nfree*(idx-1), 3*gpneb->nfree) =
                    VectorXd::Map(gpneb->projectedForceArray.at(idx).data(), 3*gpneb->nfree);
            }
            return -forceV;
        }

        double getEnergy()
        {
            double Energy=0;
            for (size_t idx{1}; idx < gpneb->imageArray.size()-1; idx++){
                auto pe_forces = gpneb->imageArray[idx].gpr_energy_forces();
                Energy += std::get<double>(pe_forces);
            }
            return Energy;
        }

        VectorXd getPositions()
        {
            VectorXd posV;
            posV.resize(3 * gpneb->nfree * gpneb->nimages);
            for (size_t idx{1}; idx < gpneb->imageArray.size()-1; idx++){
                auto image = gpneb->imageArray[idx];
                auto pe_forces = image.gpr_energy_forces();
                // NOTE: Free positions ONLY?
                posV.segment(3*gpneb->nfree*(idx-1), 3*gpneb->nfree) = VectorXd::Map(image.truePotMatter.getPositionsFree().data(), 3*gpneb->nfree);
            }
            return posV;
        }

        void setPositions(VectorXd x){
            gpneb->movedAfterForceCall = true;
            for (size_t idx{1}; idx < gpneb->imageArray.size()-1; idx++){
                // NOTE: Free positions ONLY?
                gpneb->imageArray[idx].truePotMatter.setPositionsFree(MatrixXd::Map(x.segment(3*gpneb->nfree*(idx-1),3*gpneb->nfree).data(),gpneb->nfree,3));
            }
        }

        int degreesOfFreedom() { return 3*gpneb->nimages*gpneb->nfree; }

        bool isConverged() { return getConvergence() < this->threshold; }

        double getConvergence() { return gpneb->convergenceForce(); }

        VectorXd difference(VectorXd a, VectorXd b){
            VectorXd pbcDiff(3*gpneb->nimages*gpneb->nfree);
            for (size_t idx{1}; idx < gpneb->imageArray.size()-1; idx++){
                int n = (idx-1)*3*gpneb->nfree;
                int m = 3*gpneb->nfree;
                pbcDiff.segment(n, m) = gpneb->imageArray[idx].truePotMatter.pbcV(a.segment(n,m)-b.segment(n,m));
            }
            return pbcDiff;
        }

    private:
        GPRNEB *gpneb;
        double threshold; // for converged forces
        Parameters* params;
};

GPRNEB::GPRNEB(std::vector<GPRMatter> initPath, Parameters params):
    params{params},
    imageArray{initPath},
    nebImages{initPath.begin()+1, initPath.end()-1},
    init_path_length{initPath.front().truePotMatter.distanceTo(
        initPath.back().truePotMatter)},
    nimages{params.nebImages},
    threshold{params.nebConvergedForce},
    natoms{initPath.front().truePotMatter.numberOfAtoms()},
    nfree{initPath.front().truePotMatter.numberOfFreeAtoms()},
    totImages{nimages+2}
{
    log("\nNEB: initialize\n");
    tangentArray.resize(totImages);
    projectedForceArray.resize(totImages);
    tangentArrayTrue.resize(totImages);
    projectedForceArrayTrue.resize(totImages);
    extremumEnergies.resize(2*(nimages+1));
    extremumPositions.resize(2*(nimages+1));
    extremumCurvatures.resize(2*(nimages+1));
    numExtrema = 0;
    movedAfterForceCall = true;
    climbingImage = 0;
    return;
}

GPRNEB::~GPRNEB()
{
}


int GPRNEB::compute()
{
    int status = 0;
    long iteration = 0;

    log("Nudged elastic band calculation started.\n");

    updateForces();

    GPRNEBObjectiveFunction objf(this, &this->params);

    Optimizer *optimizer = Optimizer::getOptimizer(&objf, &this->params);

    const char *forceLabel = this->params.optConvergenceMetricLabel.c_str();
    log("%10s %12s %14s %11s %12s\n", "iteration", "step size", forceLabel, "max image", "max energy");
    log("---------------------------------------------------------------\n");

    char fmt[] = "%10li %12.4e %14.4e %11li %12.4f\n";
    char fmtTiny[] = "%10li %12.4e %14.4e %11li %12.4e\n";

    while (!objf.isConverged())
    {
        if (this->params.writeMovies) {
            bool append = true;
            if (iteration == 0) append = false;
            this->imageArray.at(maxEnergyImage).truePotMatter.matter2con("neb_maximage.con", append);
        }
        // NOTE: Free positions only
        VectorXd pos = objf.getPositions();
        if(iteration) { // so that we print forces before taking an optimizer step
            if (iteration >= this->params.nebMaxIterations) {
                status = STATUS_BAD_MAX_ITERATIONS;
                break;
            }
            optimizer->step(this->params.optMaxMove);
        }
        iteration++;

        auto pef_maxei = this->imageArray[maxEnergyImage].gpr_energy_forces();
        auto pef_init = this->imageArray.front().gpr_energy_forces();
        double dE = std::get<double>(pef_maxei) -
                    std::get<double>(pef_init);
        double stepSize = helper_functions::maxAtomMotionV(imageArray.front().truePotMatter.pbcV(objf.getPositions()-pos));
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
double GPRNEB::convergenceForce()
{
    if(movedAfterForceCall) {
        updateForces();
    }
    double fmax = 0;

    for(size_t idx{1}; idx < imageArray.size()-1; idx++) {

            if( this->params.nebClimbingImageConvergedOnly == true &&
                this->params.nebClimbingImageMethod &&
                climbingImage!=0 ) {
                idx = climbingImage;
            }
            if (this->params.optConvergenceMetric == "norm") {
                fmax = std::max(fmax, projectedForceArray[idx].norm());
            } else if (this->params.optConvergenceMetric == "max_atom") {
                for (size_t jdx{0}; jdx < imageArray.front().truePotMatter.numberOfFreeAtoms(); jdx++) {
                    fmax = std::max(fmax, this->projectedForceArray[idx].row(jdx).norm());
                }
            } else if (this->params.optConvergenceMetric == "max_component") {
                fmax = max(fmax, this->projectedForceArray[idx].maxCoeff());
            } else {
                log("[Nudged Elastic Band] unknown opt_convergence_metric: %s\n", this->params.optConvergenceMetric.c_str());
                exit(1);
            }
            if( this->params.nebClimbingImageConvergedOnly == true &&
                this->params.nebClimbingImageMethod &&
                climbingImage!=0 ) {
                break;
            }
    }
    return fmax;
}

// Update the forces, do the projections, and add spring forces
void GPRNEB::updateForces()
{
    // variables for tangent
    double maxDiffEnergy, minDiffEnergy;
    double energyDiffPrev, energyDiffNext;
    double energy, energyPrev, energyNext;
    //bool higherEnergyPrev, higherEnergyNext;

    // variables for climbing image
    double maxEnergy;

    // variables for force projections
    AtomMatrix force(nfree,3), forcePerp(nfree,3), forcePar(nfree,3);
    AtomMatrix forceSpringPar(nfree,3), forceSpring(nfree,3), forceSpringPerp(nfree,3);
    AtomMatrix forceDNEB(nfree, 3);
    AtomMatrix pos(nfree,3), posNext(nfree,3), posPrev(nfree,3);
    double distNext, distPrev;

    // update the forces on the images and find the highest energy image
    auto maxelemResult = std::max_element(imageArray.begin(), imageArray.end(),
                                        helper_functions::max_gpr_energy);
    maxEnergyImage = std::distance(imageArray.begin(), maxelemResult);
    maxEnergy = std::get<double>(imageArray[maxEnergyImage].gpr_energy_forces());
    // std::cout<<"\n Got max image "<<maxEnergyImage<<std::endl;
    // std::cout<<"\n Got max energy "<<maxEnergy<<std::endl;

    for (size_t idx{1}; idx < imageArray.size()-1; idx++){
        // set local variables
        auto image = imageArray[idx]; // TODO: Maybe remove
        auto pef_cur = image.gpr_energy_forces();
        auto pef_prev = this->imageArray[idx-1].gpr_energy_forces();
        auto pef_next = this->imageArray[idx+1].gpr_energy_forces();
        force = std::get<AtomMatrix>(pef_cur);
        pos = this->imageArray[idx].truePotMatter.getPositionsFree();
        posPrev = this->imageArray[idx-1].truePotMatter.getPositionsFree();
        posNext = this->imageArray[idx+1].truePotMatter.getPositionsFree();
        energy = std::get<double>(pef_cur);
        energyPrev = std::get<double>(pef_prev);
        energyNext = std::get<double>(pef_next);

        // determine the tangent
        if(this->params.nebOldTangent) {
            // old tangent
            tangentArray[idx] = image.truePotMatter.pbc(posNext - posPrev);
        } else {
            // new improved tangent
            //higherEnergyPrev = energyPrev > energyNext;
            //higherEnergyNext = energyNext > energyPrev;
            if(energyNext > energy && energy > energyPrev) {
                tangentArray[idx] = image.truePotMatter.pbc(posNext - pos);
            }else if(energy > energyNext && energyPrev > energy){
                tangentArray[idx] = image.truePotMatter.pbc(pos - posPrev);
            }else{
                // we are at an extremum
                energyDiffPrev = energyPrev - energy;
                energyDiffNext = energyNext - energy;

                // calculate the energy difference to neighboring images
                minDiffEnergy = std::min(std::abs(energyDiffPrev), std::abs(energyDiffNext));
                maxDiffEnergy = std::max(std::abs(energyDiffPrev), std::abs(energyDiffNext));

                // use these energy differences to weight the tangent
                if(energyDiffPrev > energyDiffNext) {
                    tangentArray.at(idx) = image.truePotMatter.pbc(posNext - pos) * minDiffEnergy;
                    tangentArray.at(idx) += image.truePotMatter.pbc(pos - posPrev) * maxDiffEnergy;
                }else{
                    tangentArray.at(idx) = image.truePotMatter.pbc(posNext - pos) * maxDiffEnergy;
                    tangentArray.at(idx) += image.truePotMatter.pbc(pos - posPrev) * minDiffEnergy;
                }
            }
        }
        tangentArray[idx].normalize();

        // project the forces and add springs

        // calculate the force perpendicular to the tangent
        forcePerp = force - (force.array() * (tangentArray.at(idx)).array()).sum() * tangentArray.at(idx);
        forceSpring = this->params.nebSpring * image.truePotMatter.pbc((posNext - pos) - (pos - posPrev));

        // calculate the spring force
        distPrev = image.truePotMatter.pbc(posPrev - pos).squaredNorm();
        distNext = image.truePotMatter.pbc(posNext - pos).squaredNorm();
        forceSpringPar = this->params.nebSpring * (distNext-distPrev) * tangentArray[idx];

        if (this->params.nebDoublyNudged) {
            forceSpringPerp = forceSpring - (forceSpring.array() * (tangentArray[idx]).array()).sum() * tangentArray[idx];
            forceDNEB = forceSpringPerp - (forceSpringPerp.array() * forcePerp.normalized().array()).sum() * forcePerp.normalized();
            if (this->params.nebDoublyNudgedSwitching) {
                double switching;
                switching = 2.0/M_PI * std::atan(std::pow(forcePerp.norm(),2) / std::pow(forceSpringPerp.norm(),2));
                forceDNEB *= switching;
            }
        }else{
            forceDNEB.setZero();
        }


        if(this->params.nebClimbingImageMethod && idx==maxEnergyImage)
        {
            // we are at the climbing image
            climbingImage = maxEnergyImage;
            this->projectedForceArray[idx] = force - (2.0 * (force.array() * (tangentArray[idx]).array()).sum() * tangentArray[idx]) + forceDNEB;
        }
        else  // all non-climbing images
        {
            // sum the spring and potential forces for the neb force
            if (this->params.nebElasticBand) {
                this->projectedForceArray[idx] = forceSpring + force;
                // std::cout<<"Spring forces "<<forceSpring+force<<"\n";
            }else{
                this->projectedForceArray[idx] = forceSpringPar + forcePerp + forceDNEB;
                // std::cout<<"Spring forces "<<forceSpring+force+forceDNEB<<"\n";
            }
                // std::cout<<"\n Projected forces "<<projectedForceArray[idx]<<"\n";
            //*projectedForce[i] = forceSpring + forcePerp;
            //if (this->params.nebFullSpring) {
            movedAfterForceCall = false;  // so that we don't repeat a force call
        }
        //zero net translational force
        if (image.truePotMatter.numberOfAtoms() == image.truePotMatter.numberOfFreeAtoms()) {
            for (size_t jdx{0};jdx <= 2; jdx++) {
                double translationMag = projectedForceArray[idx].col(jdx).sum();
                int num_atoms = projectedForceArray[idx].col(jdx).size();
                this->projectedForceArray[idx].col(jdx).array() -= translationMag/(static_cast<double>(num_atoms));
                // std::cout<<"Final Projc "<<this->projectedForceArray[idx].col(jdx).array()<<"\n";
            }
        }
    }
    // Never move the first and last images
    // NOTE: This behaviour is **different** from that of vanilla EON
    this->projectedForceArray.front().setZero();
    this->projectedForceArray.back().setZero();
    return;
}

void GPRNEB::getTrueNEBForces()
{
    // variables for tangent
    double maxDiffEnergy, minDiffEnergy;
    double energyDiffPrev, energyDiffNext;
    double energy, energyPrev, energyNext;
    //bool higherEnergyPrev, higherEnergyNext;

    // variables for climbing image
    double maxEnergy;

    // variables for force projections
    AtomMatrix force(natoms,3), forcePerp(natoms,3), forcePar(natoms,3);
    AtomMatrix forceSpringPar(natoms,3), forceSpring(natoms,3), forceSpringPerp(natoms,3);
    AtomMatrix forceDNEB(natoms, 3);
    AtomMatrix pos(natoms,3), posNext(natoms,3), posPrev(natoms,3);
    double distNext, distPrev;

    // update the forces on the images and find the highest energy image
    auto maxelemResult = std::max_element(imageArray.begin(), imageArray.end(),
                                        helper_functions::max_true_energy);
    maxEnergyImage = std::distance(imageArray.begin(), maxelemResult);
    maxEnergy = std::get<double>(imageArray[maxEnergyImage].true_energy_forces());

    for (size_t idx{1}; idx < imageArray.size()-1; idx++){
        // set local variables
        auto image = imageArray[idx]; // TODO: Maybe remove
        auto pef_cur = image.true_energy_forces();
        auto pef_prev = this->imageArray[idx-1].true_energy_forces();
        auto pef_next = this->imageArray[idx+1].true_energy_forces();
        force = std::get<AtomMatrix>(pef_cur);
        pos = this->imageArray[idx].truePotMatter.getPositions();
        posPrev = this->imageArray[idx-1].truePotMatter.getPositions();
        posNext = this->imageArray[idx+1].truePotMatter.getPositions();
        energy = std::get<double>(pef_cur);
        energyPrev = std::get<double>(pef_prev);
        energyNext = std::get<double>(pef_next);

        // determine the tangent
        if(this->params.nebOldTangent) {
            // old tangent
            tangentArrayTrue[idx] = image.truePotMatter.pbc(posNext - posPrev);
        } else {
            // new improved tangent
            //higherEnergyPrev = energyPrev > energyNext;
            //higherEnergyNext = energyNext > energyPrev;
            if(energyNext > energy && energy > energyPrev) {
                tangentArrayTrue[idx] = image.truePotMatter.pbc(posNext - pos);
            }else if(energy > energyNext && energyPrev > energy){
                tangentArrayTrue[idx] = image.truePotMatter.pbc(pos - posPrev);
            }else{
                // we are at an extremum
                energyDiffPrev = energyPrev - energy;
                energyDiffNext = energyNext - energy;

                // calculate the energy difference to neighboring images
                minDiffEnergy = std::min(std::abs(energyDiffPrev), std::abs(energyDiffNext));
                maxDiffEnergy = std::max(std::abs(energyDiffPrev), std::abs(energyDiffNext));

                // use these energy differences to weight thise tangent
                if(energyDiffPrev > energyDiffNext) {
                    tangentArrayTrue.at(idx) = image.truePotMatter.pbc(posNext - pos) * minDiffEnergy;
                    tangentArrayTrue.at(idx) += image.truePotMatter.pbc(pos - posPrev) * maxDiffEnergy;
                }else{
                    tangentArrayTrue.at(idx) = image.truePotMatter.pbc(posNext - pos) * maxDiffEnergy;
                    tangentArrayTrue.at(idx) += image.truePotMatter.pbc(pos - posPrev) * minDiffEnergy;
                }
            }
        }
        tangentArrayTrue[idx].normalize();

        // project the forces and add springs

        // calculate the force perpendicular to the tangent
        forcePerp = force - (force.array() * (tangentArrayTrue.at(idx)).array()).sum() * tangentArrayTrue.at(idx);
        forceSpring = this->params.nebSpring * image.truePotMatter.pbc((posNext - pos) - (pos - posPrev));

        // calculate the spring force
        distPrev = image.truePotMatter.pbc(posPrev - pos).squaredNorm();
        distNext = image.truePotMatter.pbc(posNext - pos).squaredNorm();
        forceSpringPar = this->params.nebSpring * (distNext-distPrev) * tangentArrayTrue[idx];

        if (this->params.nebDoublyNudged) {
            forceSpringPerp = forceSpring - (forceSpring.array() * (tangentArrayTrue[idx]).array()).sum() * tangentArrayTrue[idx];
            forceDNEB = forceSpringPerp - (forceSpringPerp.array() * forcePerp.normalized().array()).sum() * forcePerp.normalized();
            if (this->params.nebDoublyNudgedSwitching) {
                double switching;
                switching = 2.0/M_PI * std::atan(std::pow(forcePerp.norm(),2) / std::pow(forceSpringPerp.norm(),2));
                forceDNEB *= switching;
            }
        }else{
            forceDNEB.setZero();
        }


        if(this->params.nebClimbingImageMethod && idx==maxEnergyImage)
        {
            // we are at the climbing image
            climbingImage = maxEnergyImage;
            this->projectedForceArrayTrue[idx] = force - (2.0 * (force.array() * (tangentArrayTrue[idx]).array()).sum() * tangentArrayTrue[idx]) + forceDNEB;
        }
        else  // all non-climbing images
        {
            // sum the spring and potential forces for the neb force
            if (this->params.nebElasticBand) {
                this->projectedForceArrayTrue[idx] = forceSpring + force;
                // std::cout<<"Spring forces "<<forceSpring+force<<"\n";
            }else{
                this->projectedForceArrayTrue[idx] = forceSpringPar + forcePerp + forceDNEB;
                // std::cout<<"Spring forces "<<forceSpring+force+forceDNEB<<"\n";
            }
                // std::cout<<"\n Projected forces "<<projectedForceArrayTrue[idx]<<"\n";
            //*projectedForce[i] = forceSpring + forcePerp;
            //if (this->params.nebFullSpring) {
            movedAfterForceCall = false;  // so that we don't repeat a force call
        }

        //zero net translational force
        if (image.truePotMatter.numberOfAtoms() == image.truePotMatter.numberOfFreeAtoms()) {
            for (size_t jdx{0};jdx <= 2; jdx++) {
                double translationMag = projectedForceArrayTrue[idx].col(jdx).sum();
                int num_atoms = projectedForceArrayTrue[idx].col(jdx).size();
                this->projectedForceArrayTrue[idx].col(jdx).array() -= translationMag/(static_cast<double>(num_atoms));
            }
        }
    }
    return;
}
double GPRNEB::getTrueConvForce()
{
    getTrueNEBForces();
    double fmax = 0;

    for(size_t idx{1}; idx < imageArray.size()-1; idx++) {

            if( this->params.nebClimbingImageConvergedOnly == true &&
                this->params.nebClimbingImageMethod &&
                climbingImage!=0 ) {
                idx = climbingImage;
            }
            if (this->params.optConvergenceMetric == "norm") {
                fmax = std::max(fmax, projectedForceArrayTrue[idx].norm());
            } else if (this->params.optConvergenceMetric == "max_atom") {
                for (size_t jdx{0}; jdx < imageArray.front().truePotMatter.numberOfFreeAtoms(); jdx++) {
                    fmax = std::max(fmax, this->projectedForceArrayTrue[idx].row(jdx).norm());
                }
            } else if (this->params.optConvergenceMetric == "max_component") {
                fmax = max(fmax, this->projectedForceArrayTrue[idx].maxCoeff());
            } else {
                log("[Nudged Elastic Band] unknown opt_convergence_metric: %s\n", this->params.optConvergenceMetric.c_str());
                exit(1);
            }
            if( this->params.nebClimbingImageConvergedOnly == true &&
                this->params.nebClimbingImageMethod &&
                climbingImage!=0 ) {
                break;
            }
    }
    return fmax;
}

// Print NEB image data
void GPRNEB::printImageData(bool writeToFile)
{
    double dist, distTotal=0;
    AtomMatrix tangentStart = imageArray.front().truePotMatter.pbc(imageArray[1].truePotMatter.getPositionsFree() - imageArray.front().truePotMatter.getPositionsFree());
    AtomMatrix tangentEnd = imageArray[nimages].truePotMatter.pbc(imageArray.back().truePotMatter.getPositionsFree() - imageArray[nimages].truePotMatter.getPositionsFree());
    AtomMatrix tang;

    log("Image data (as in neb.dat)\n");

    FILE *fh=NULL;
    if (writeToFile) {
        fh = fopen("neb.dat", "w");
    }

    for(size_t idx{0}; idx <= nimages+1; idx++)
    {
        if(idx==0){ tang = tangentStart; }
        else if (idx==nimages+1) { tang = tangentEnd; }
        else { tang = tangentArray[idx]; }
        if(idx>0) {
            dist = imageArray[idx].truePotMatter.distanceTo(imageArray[idx-1].truePotMatter);
            distTotal += dist;
        }
        auto pef_cur = this->imageArray[idx].gpr_energy_forces();
        auto pef_init = this->imageArray.front().gpr_energy_forces();
        if (fh == NULL) {
            log("%3li %12.6f %12.6f %12.6f\n", idx, distTotal,
                std::get<double>(pef_cur)-std::get<double>(pef_init), (std::get<AtomMatrix>(pef_cur).array()*tang.array()).sum());
        }else{
            fprintf(fh, "%3li %12.6f %12.6f %12.6f\n",idx,distTotal,
                std::get<double>(pef_cur)-std::get<double>(pef_init), (std::get<AtomMatrix>(pef_cur).array()*tang.array()).sum());
        }
    }
    if (writeToFile) {
        fclose(fh);
    }
}

// Estimate the barrier using a cubic spline
void GPRNEB::findExtrema()
{
    // calculate the cubic parameters for each interval (a,b,c,d)

    AtomMatrix tangentEndpoint;
    std::vector<double> a, b, c, d;
    a.resize(nimages+1);
    b.resize(nimages+1);
    c.resize(nimages+1);
    d.resize(nimages+1);
    double F1, F2, U1, U2, dist;

    for(size_t idx{0}; idx<=nimages; idx++)
    {
        dist = imageArray[idx].truePotMatter.distanceTo(imageArray[idx+1].truePotMatter);
        auto pef_f1 = this->imageArray[idx].gpr_energy_forces();
        auto pef_f2 = this->imageArray[idx+1].gpr_energy_forces();
        if(idx==0) {
            tangentEndpoint = imageArray[idx].truePotMatter.pbc(imageArray[1].truePotMatter.getPositionsFree() -
                                                                imageArray.front().truePotMatter.getPositionsFree());
            tangentEndpoint.normalize();
            F1 = (std::get<AtomMatrix>(pef_f1).array()*tangentEndpoint.array()).sum()*dist;
        } else {
            F1 = (std::get<AtomMatrix>(pef_f1).array()*(tangentArray[idx]).array()).sum()*dist;
        }
        if(idx==nimages) {
            tangentEndpoint =  imageArray[idx+1].truePotMatter.pbc(imageArray[nimages+1].truePotMatter.getPositionsFree() - imageArray[nimages].truePotMatter.getPositionsFree());
            tangentEndpoint.normalize();
            F2 = (std::get<AtomMatrix>(pef_f2).array()*tangentEndpoint.array()).sum()*dist;
        } else {
            F2 = (std::get<AtomMatrix>(pef_f2).array()*(tangentArray[idx+1]).array()).sum()*dist;
        }
        U1 = std::get<double>(pef_f1);
        U2 = std::get<double>(pef_f2);
        a[idx] = U1;
        b[idx] = -F1;
        c[idx] = 3.*(U2 - U1) + 2.*F1 + F2;
        d[idx] = -2.*(U2 - U1) - (F1 + F2);
    }

    // finding extrema along the MEP

//    long numExtrema = 0;
//    double extremaEnergy[2*(images+1)]; // the maximum number of extrema
//    double extremaPosition[2*(images+1)];
    double discriminant, f;

    for(size_t idx{0}; idx<=nimages; idx++)
    {
        discriminant = std::pow(c[idx],2) - 3.*b[idx]*d[idx];
        if(discriminant >= 0) {
            f = -1;

            // quadratic case
            if( (d[idx] == 0) && (c[idx] != 0) ) {
                f = ( -b[idx]/(2.*c[idx]) );
            }
            // cubic case 1
            else if( d[idx] != 0 ) {
                f = -(c[idx] + std::sqrt(discriminant))/(3.*d[idx]);
            }
            if( (f >= 0) && (f <= 1) ) {
                extremumPositions[numExtrema] = idx + f;
                extremumEnergies[numExtrema] = d[idx]*std::pow(f,3) + c[idx]*std::pow(f,2) + b[idx]*f + a[idx];
                extremumCurvatures[numExtrema] = 6.0*d[idx]*f + 2*c[idx];
                numExtrema++;
            }
            // cubic case 2
            if( d[idx] != 0 ) {
                f = ( -(c[idx] - std::sqrt(discriminant))/(3.*d[idx]) );
            }
            if( (f >= 0) && (f <= 1) ) {
                extremumPositions[numExtrema] = idx + f;
                extremumEnergies[numExtrema] = d[idx]*std::pow(f,3) + c[idx]*std::pow(f,2) + b[idx]*f + a[idx];
                extremumCurvatures[numExtrema] = 6*d[idx]*f + 2*c[idx];
                numExtrema++;
            }
        }
    }

    auto pef_init = this->imageArray.front().gpr_energy_forces();
    log("\nFound %li extrema\n",numExtrema);
    log("Energy reference: %f\n", std::get<double>(pef_init));
    for(size_t idx{0}; idx<numExtrema; idx++) {
        log("extrema #%li at image position %f with energy %f and curvature %f\n",idx+1,extremumPositions[idx],extremumEnergies[idx]-std::get<double>(pef_init), extremumCurvatures[idx]);
    }
}

// TODO: Not here, but can be used to train PES
// NOTE: used to also use img.areEnergiesCloseToTrue()
bool GPRNEB::needsRetraining(double eps){
    bool retval{false};
    // NOTE: Assumes that the final and end points are relaxed
    for (size_t idx{1}; idx < this->imageArray.size()-1; idx++){
        // TODO: Use other convergence criteria
        if (imageArray[idx].isForceMaxElementLow(eps)){
            retval = false;
        } else {
            retval = true;
        }
    }
    return retval;
}

bool GPRNEB::stoppedEarly(std::vector<Matter> prevPath, double max_dist_factor){
    bool noEarlyStopping{false};
    double max_dist = max_dist_factor * this->init_path_length;
    auto curpath = this->getCurPath();
    for (size_t idx{1}; idx < prevPath.size()-1; idx++){
    auto img = prevPath[idx];
    std::vector<double> distances;
    std::transform(curpath.begin(),
                   curpath.end(),
                   std::back_inserter(distances),
                   [&](Matter mat)->double{
                       return mat.distanceTo(img);
                   });
    // Remove distance to self
    distances.erase(std::remove(distances.begin(), distances.end(), 0), distances.end());
    noEarlyStopping = std::any_of(distances.begin(),
                         distances.end(),
                         [&](const double dist)->bool{
                             return (std::abs(dist) < max_dist); });
    if (!noEarlyStopping){
        std::cout<<"\n EARLY STOPPING\n";
        return !noEarlyStopping;
    }
    }
    return !noEarlyStopping;
}

bool GPRNEB::notStoppedEarly(std::vector<Matter> prevPath, double max_dist_factor){
    return not stoppedEarly(prevPath, max_dist_factor);
}

std::vector<Matter> GPRNEB::getCurPathFull(){
    std::vector<Matter> matvec;
    for (size_t idx{0}; idx < this->imageArray.size(); idx++){
        matvec.push_back(imageArray[idx].truePotMatter);
    }
    return matvec;
}

std::vector<Matter> GPRNEB::getCurPath(){
    std::vector<Matter> matvec;
    // NOTE: Intermediates only, assumes relaxed final, end points
    for (size_t idx{1}; idx < this->imageArray.size()-1; idx++){
        matvec.push_back(imageArray[idx].truePotMatter);
    }
    return matvec;
}
