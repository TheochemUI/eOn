/*
 *===============================================
 *  EON Prefactors.cpp
 *===============================================
 */
#include "Prefactors.h"

using namespace helper_functions;


Prefactors::Prefactors(){
    eigenValMin1_ = 0;
    eigenValMin2_ = 0;
    eigenValSaddle_ = 0;
    coordinatesToAccountFor_ = 0;
    totalForceCalls = 0;
}


Prefactors::Prefactors(const Matter *saddle, const Matter *min1, 
                       const Matter *min2, Parameters *parameters){
    eigenValMin1_ = 0;
    eigenValMin2_ = 0;
    eigenValSaddle_ = 0;
    coordinatesToAccountFor_ = 0;
    initialize(saddle, min1, min2, parameters);
    return;
}


Prefactors::~Prefactors(){
    clean();
    return;
}


void Prefactors::clean(){
    // min1_, min2_ and saddle_ should not be deleted
    // copies of pointers that were passed in! 
    if(eigenValMin1_ != 0){
        delete [] eigenValMin1_;
        eigenValMin1_ = 0;
    }
    if(eigenValMin2_ != 0){
        delete [] eigenValMin2_;
        eigenValMin2_ = 0;
    }
    if(eigenValSaddle_ != 0){
        delete [] eigenValSaddle_;
        eigenValSaddle_ = 0;
    }
    if(coordinatesToAccountFor_ != 0){
        delete [] coordinatesToAccountFor_;
        coordinatesToAccountFor_ = 0;
    }
    
    return;
}


void Prefactors::initialize(const Matter *saddle, const Matter *min1, 
                            const Matter *min2, Parameters *parameters){
    clean();
    saddle_ = saddle;
    min1_ = min1;
    min2_ = min2;
    parameters_ = parameters;
    
    nAtoms_ = saddle->numberOfAtoms();
    coordinatesToAccountFor_ = new bool[3*nAtoms_];
    return;
}


bool Prefactors::compute(double *prefactors){
    bool good = false;
    
    long forceCallsSaddle;
    long forceCallsMin1;
    long forceCallsMin2;

    sizeHessian_ = atomsToAccountForInHessian();
    fprintf(stdout, "Hessian size %li\n",sizeHessian_);
    if (sizeHessian_ == 0) {
        ///XXX: This should not exit with a status of 1
        //      due to showing up as a boinc error.
        //      Users do not get credit for non-zero
        //      exit statuses.
        fprintf(stderr, "Error: size of hessian is zero. "
                "Try with smaller min_Displacement_Hessian\n");
        exit(1);
    };    
    eigenValMin1_ = new double[sizeHessian_];
    eigenValMin2_ = new double[sizeHessian_];
    eigenValSaddle_ = new double[sizeHessian_];    
    
    forceCallsSaddle = saddle_->getForceCalls();
    forceCallsMin1 = min1_->getForceCalls();
    forceCallsMin2 = min2_->getForceCalls();    
    
    good = getEigenValues();
    
    if(good){
        // Calculating the prefactor for the forward and backward process
        prefactors[0] = 1;
        prefactors[1] = 1;
        
        for(int i=0; i<sizeHessian_; i++){
            prefactors[0] = prefactors[0]*eigenValMin1_[i];
            prefactors[1] = prefactors[1]*eigenValMin2_[i];
            if(eigenValSaddle_[i]>0){
                prefactors[0] = prefactors[0]/eigenValSaddle_[i];
                prefactors[1] = prefactors[1]/eigenValSaddle_[i];
            }
        }
        // 10.18e-15 conversion factor from au to 1/s
        prefactors[0] = sqrt(prefactors[0])/(2*M_PI*10.18e-15);
        prefactors[1] = sqrt(prefactors[1])/(2*M_PI*10.18e-15);
    }
    return(good);
}

long Prefactors::atomsToAccountForInHessian(){
	long sizeHessian;
	double minDisp;
	// Account for all atoms moved more than the value specified in parameters
	if (parameters_->hessianMaxSize == 0){
		minDisp = parameters_->hessianMinDisplacement;
		sizeHessian = atomsMovedMoreThan(minDisp);
	}
	// Will ensure that there is not accounted for more
	// than a specified number of atoms
	else{
		int loop = 0;
		do {
			minDisp = parameters_->hessianMinDisplacement + loop * 0.1;		
			sizeHessian = atomsMovedMoreThan(minDisp);
			loop = loop + 1;
		} while (parameters_->hessianMaxSize < sizeHessian);
	}
	return(sizeHessian);
}

long Prefactors::atomsMovedMoreThan(double minDisplacement){
    long sizeHessian;
    double diffR1, diffR2, diffRSaddle;
    for(int i=0; i<3*nAtoms_; i++)
        coordinatesToAccountFor_[i] = false;
    //----- Initialize end -----
    //std::cout<<"determineActiveAtoms\n";
    
    // Picking out all atoms that are displaced
    for(int i=0; i<nAtoms_; i++){
        diffR1 = saddle_->distance(*min1_, i);
        diffR2 = saddle_->distance(*min2_, i);
        if(((minDisplacement<diffR1) || 
			(minDisplacement<diffR2)) &&
		   !saddle_->getFixed(i)){
            coordinatesToAccountFor_[ 3*i ] = true;
            coordinatesToAccountFor_[3*i+1] = true;
            coordinatesToAccountFor_[3*i+2] = true;
            
            // Picking out free atoms in the vicinity of a displaced atom
            for(int j=0; j<nAtoms_; j++){
                diffRSaddle = saddle_->distance(i,j);
                
                if(diffRSaddle<parameters_->hessianWithinRadiusDisplaced 
                   && (!saddle_->getFixed(j))){
                    coordinatesToAccountFor_[ 3*j ] = true;
                    coordinatesToAccountFor_[3*j+1] = true;
                    coordinatesToAccountFor_[3*j+2] = true;
                }
            }
        }
    }
    // Counting all the atoms to be accounted for in the Hessian
    sizeHessian = 0;
    for(int i=0; i<3*nAtoms_; i++){
        if(coordinatesToAccountFor_[i] == true)
            sizeHessian = sizeHessian+1;
    }
    return(sizeHessian);
}
	

void Prefactors::determineHessian(double **hessian, const Matter *matter){
    long iCoord, jCoord;
    long iHessian, jHessian;
    long nCoord = 3*nAtoms_;
    long forceCallsAtomsTemp;
    
    Matter matterTemp(parameters_);
    matterTemp = *matter;
    
    double *pos, *posTemp, *posDisplace;
    double *force1, *force2;
    
    pos = new double[nCoord];
    posTemp = new double[nCoord];
    posDisplace = new double[nCoord];
    force1 = new double[nCoord];
    force2 = new double[nCoord];
    
    double dr = 1e-6; // value used by graeme 1e-4;
    double tdr = 2*dr;
    matter->getPositions(pos);
    forceCallsAtomsTemp = matterTemp.getForceCalls();

    //----- Initialize end -----
    //std::cout<<"determineHessian\n";
    
    for(iCoord=iHessian=0; iHessian<sizeHessian_; iHessian++, iCoord++){
        for(int k=0; k<nCoord; k++)
            posDisplace[k] = 0;
        
        // Getting the index for the next coordinate to be accounted for
        while(!coordinatesToAccountFor_[iCoord]) 
            iCoord++;
        
        // Displacing one coordinate
        posDisplace[iCoord] = dr;
        
        subtract(posTemp,pos,posDisplace,nCoord);
        matterTemp.setPositions(posTemp);
        matterTemp.getForces(force1);
        
        add(posTemp,pos,posDisplace,nCoord);
        matterTemp.setPositions(posTemp);
        matterTemp.getForces(force2);
        
        for(jCoord=jHessian=0; jHessian<sizeHessian_; jHessian++,jCoord++){
            
            // Getting the index for the next coordinate to be accounted for
            while(!coordinatesToAccountFor_[jCoord]) 
                jCoord++;
            hessian[jHessian][iHessian] = -(force2[jCoord]-force1[jCoord])/tdr;
        }
    }
    
    forceCallsAtomsTemp = matterTemp.getForceCalls()-forceCallsAtomsTemp;
    
    totalForceCalls += forceCallsAtomsTemp;

    delete [] pos;
    delete [] posTemp;
    delete [] posDisplace;
    delete [] force1;
    delete [] force2;
    return;
}


bool Prefactors::getEigenValues(){
    assert(sizeHessian_>0);
    bool result = true;
    double **hessian;
    hessian = new double*[sizeHessian_];
    for(int i=0; i<sizeHessian_; i++)
        hessian[i] = new double[sizeHessian_];
    long negModesInSaddle = 0;
    //----- Initialize end -----
    //std::cout<<"getEigenValues\n";
    
    // Eigenvalues for minima1
    determineHessian(hessian, min1_);
    massScaleHessian(hessian);   
    
    eigenValues(sizeHessian_, eigenValMin1_, hessian);
    for(int i=0; i<sizeHessian_; i++){
        if(eigenValMin1_[i]<=0){
            result = false;
            break;
        }
    }
    if(result){
        // Eigenvalues for minima2
        determineHessian(hessian, min2_);
        massScaleHessian(hessian);
        
        eigenValues(sizeHessian_, eigenValMin2_, hessian);
        for(int i=0; i<sizeHessian_; i++){
            if(eigenValMin2_[i]<=0){
                result = false;
                break;
            }
        }
    }    
    if(result){
        // Eigenvalues for saddle
        determineHessian(hessian, saddle_);
        massScaleHessian(hessian);
        
        eigenValues(sizeHessian_, eigenValSaddle_, hessian);
        for(int i=0; i<sizeHessian_; i++){
            if(eigenValSaddle_[i]<0){
                negModesInSaddle = negModesInSaddle+1;
            }
        }
        if(negModesInSaddle!=1)
            result = false;
    }
    for(int i=0; i<sizeHessian_; i++)
        delete [] hessian[i];
    delete hessian;
    
    return(result); 
}


void Prefactors::massScaleHessian(double **hessian){
    long iCoord, jCoord;
    long iHessian, jHessian;
    double effMass;
    
    for(iCoord=iHessian=0; iHessian<sizeHessian_; iHessian++, iCoord++){
        
        // Getting the index for the next atom to be accounted for
        while(!coordinatesToAccountFor_[iCoord]) 
            iCoord++;
        
        for(jCoord=jHessian=0; jHessian<sizeHessian_; jHessian++,jCoord++){
            
            // Getting the index for the next atom to be accounted for
            while(!coordinatesToAccountFor_[jCoord]) 
                jCoord++;
            
            // Remember the omega = sqrt(k/m)
            effMass=sqrt(saddle_->getMass(jCoord/3)*saddle_->getMass(iCoord/3));
//            effMass = effMass/SU::G_PER_MOL;
            
            hessian[jHessian][iHessian] = hessian[jHessian][iHessian]/effMass;
        }
    }
    return;
}
