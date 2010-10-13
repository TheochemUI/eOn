/*
 *===============================================
 *  EON Prefactors.cpp
 *===============================================
 */
#include "Prefactors.h"

using namespace helper_functions;


Prefactors::Prefactors(){
    eigenValMin1.setZero();
    eigenValMin2.setZero();
    eigenValSaddle.setZero();
    coordinatesToAccountFor = 0;
    totalForceCalls = 0;
}


Prefactors::Prefactors(const Matter *saddle, const Matter *min1, 
                       const Matter *min2, Parameters *parameters){
    eigenValMin1.setZero();
    eigenValMin2.setZero();
    eigenValSaddle.setZero();
    coordinatesToAccountFor = 0;
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
    if(coordinatesToAccountFor != 0){
        delete [] coordinatesToAccountFor;
        coordinatesToAccountFor = 0;
    }
    
    return;
}


void Prefactors::initialize(const Matter *saddle_passed, const Matter *min1_passed, 
                            const Matter *min2_passed, Parameters *parameters_passed){
    clean();
    saddle = saddle_passed;
    min1 = min1_passed;
    min2 = min2_passed;
    parameters = parameters_passed;
    
    nAtoms = saddle->numberOfAtoms();
    coordinatesToAccountFor = new bool[3*nAtoms];
    return;
}


bool Prefactors::compute(double *prefactors){
    bool good = false;
    
    long forceCallsSaddle;
    long forceCallsMin1;
    long forceCallsMin2;

    sizeHessian = atomsToAccountForInHessian();
    fprintf(stdout, "Hessian size %li\n",sizeHessian);
    if (sizeHessian == 0) {
        ///XXX: This should not exit with a status of 1
        //      due to showing up as a boinc error.
        //      Users do not get credit for non-zero
        //      exit statuses.
        fprintf(stderr, "Error: size of hessian is zero. "
                "Try with smaller min_Displacement_Hessian\n");
        exit(1);
    };    
    eigenValMin1.resize(sizeHessian);
    eigenValMin2.resize(sizeHessian);
    eigenValSaddle.resize(sizeHessian);    
    
    forceCallsSaddle = saddle->getForceCalls();
    forceCallsMin1 = min1->getForceCalls();
    forceCallsMin2 = min2->getForceCalls();    
    
    good = getEigenValues();
    
    if(good){
        // Calculating the prefactor for the forward and backward process
        prefactors[0] = 1;
        prefactors[1] = 1;
        
        for(int i=0; i<sizeHessian; i++){
            prefactors[0] *= eigenValMin1[i];
            prefactors[1] *= eigenValMin2[i];
            if(eigenValSaddle[i]>0){
                prefactors[0] /= eigenValSaddle[i];
                prefactors[1] /= eigenValSaddle[i];
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
	if (parameters->hessianMaxSize == 0){
		minDisp = parameters->hessianMinDisplacement;
		sizeHessian = atomsMovedMoreThan(minDisp);
	}
	// Will ensure that there is not accounted for more
	// than a specified number of atoms
	else{
		int loop = 0;
		do {
			minDisp = parameters->hessianMinDisplacement + loop * 0.1;		
			sizeHessian = atomsMovedMoreThan(minDisp);
			loop = loop + 1;
		} while (parameters->hessianMaxSize < sizeHessian);
	}
	return(sizeHessian);
}

long Prefactors::atomsMovedMoreThan(double minDisplacement){
    long sizeHessian;
    double diffR1, diffR2, diffRSaddle;
    for(int i=0; i<3*nAtoms; i++)
        coordinatesToAccountFor[i] = false;
    //----- Initialize end -----
    //std::cout<<"determineActiveAtoms\n";
    
    // Picking out all atoms that are displaced
    for(int i=0; i<nAtoms; i++){
        diffR1 = saddle->distance(*min1, i);
        diffR2 = saddle->distance(*min2, i);
        if(((minDisplacement<diffR1) || 
			(minDisplacement<diffR2)) &&
		   !saddle->getFixed(i)){
            coordinatesToAccountFor[ 3*i ] = true;
            coordinatesToAccountFor[3*i+1] = true;
            coordinatesToAccountFor[3*i+2] = true;
            
            // Picking out free atoms in the vicinity of a displaced atom
            for(int j=0; j<nAtoms; j++){
                diffRSaddle = saddle->distance(i,j);
                
                if(diffRSaddle<parameters->hessianWithinRadiusDisplaced 
                   && (!saddle->getFixed(j))){
                    coordinatesToAccountFor[ 3*j ] = true;
                    coordinatesToAccountFor[3*j+1] = true;
                    coordinatesToAccountFor[3*j+2] = true;
                }
            }
        }
    }
    // Counting all the atoms to be accounted for in the Hessian
    sizeHessian = 0;
    for(int i=0; i<3*nAtoms; i++){
        if(coordinatesToAccountFor[i] == true)
            sizeHessian = sizeHessian+1;
    }
    return(sizeHessian);
}
	

void Prefactors::determineHessian(Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian, const Matter *matter){
    long iCoord, jCoord;
    long iHessian, jHessian;
    long forceCallsAtomsTemp;
    
    Matter matterTemp(parameters);
    matterTemp = *matter;
    
    double dr = 1e-6; // value used by graeme 1e-4;
    double tdr = 2*dr;
    
    Matrix<double, Eigen::Dynamic, 3> pos = matter->getPositions();
    Matrix<double, Eigen::Dynamic, 3> posDisplace(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> posTemp(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> force1(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> force2(nAtoms, 3);

    forceCallsAtomsTemp = matterTemp.getForceCalls();

    //----- Initialize end -----
    //std::cout<<"determineHessian\n";
   
    for(iCoord=iHessian=0; iHessian<sizeHessian; iHessian++, iCoord++){
        
        posDisplace.setZero(); 
        // Getting the index for the next coordinate to be accounted for
        while(!coordinatesToAccountFor[iCoord]) 
            iCoord++;
        
        // Displacing one coordinate
        posDisplace(iCoord/3,iCoord%3)  = dr;
        
        posTemp = pos - posDisplace;
        matterTemp.setPositions(posTemp);
        force1 = matterTemp.getForces();
        
        posTemp = pos + posDisplace; 
        matterTemp.setPositions(posTemp);
        force1 = matterTemp.getForces();
        
        for(jCoord=jHessian=0; jHessian<sizeHessian; jHessian++,jCoord++){
            
            // Getting the index for the next coordinate to be accounted for
            while(!coordinatesToAccountFor[jCoord]) 
                jCoord++;
            hessian(jHessian,iHessian) = -(force2(jCoord/3, jCoord%3)-force1(jCoord/3, jCoord%3))/tdr;
        }
    }
    
    forceCallsAtomsTemp = matterTemp.getForceCalls()-forceCallsAtomsTemp;
    
    totalForceCalls += forceCallsAtomsTemp;

    return;
}


bool Prefactors::getEigenValues(){
    assert(sizeHessian>0);
    bool result = true;

    Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian((int)sizeHessian, (int)sizeHessian);

    long negModesInSaddle = 0;
    //----- Initialize end -----
    //std::cout<<"getEigenValues\n";
    
    // Eigenvalues for minima1
    determineHessian(hessian, min1);
    massScaleHessian(hessian);   
    
    eigenValMin1 = hessian.eigenvalues();
    for(int i=0; i<sizeHessian; i++){
        if(eigenValMin1[i]<=0){
            result = false;
            break;
        }
    }
    if(result){
        // Eigenvalues for minima2
        determineHessian(hessian, min2);
        massScaleHessian(hessian);
        
        eigenValMin2 = hessian.eigenvalues();
        for(int i=0; i<sizeHessian; i++){
            if(eigenValMin2[i]<=0){
                result = false;
                break;
            }
        }
    }    
    if(result){
        // Eigenvalues for saddle
        determineHessian(hessian, saddle);
        massScaleHessian(hessian);
        
        eigenValSaddle = hessian.eigenvalues();
        for(int i=0; i<sizeHessian; i++){
            if(eigenValSaddle[i]<0){
                negModesInSaddle = negModesInSaddle+1;
            }
        }
        if(negModesInSaddle!=1)
            result = false;
    }
    
    return(result); 
}


void Prefactors::massScaleHessian(Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian){
    long iCoord, jCoord;
    long iHessian, jHessian;
    double effMass;
    
    for(iCoord=iHessian=0; iHessian<sizeHessian; iHessian++, iCoord++){
        
        // Getting the index for the next atom to be accounted for
        while(!coordinatesToAccountFor[iCoord]) 
            iCoord++;
        
        for(jCoord=jHessian=0; jHessian<sizeHessian; jHessian++,jCoord++){
            
            // Getting the index for the next atom to be accounted for
            while(!coordinatesToAccountFor[jCoord]) 
                jCoord++;
            
            // Remember the omega = sqrt(k/m)
            effMass=sqrt(saddle->getMass(jCoord/3)*saddle->getMass(iCoord/3));
//            effMass = effMass/SU::G_PER_MOL;
            
            hessian(jHessian, iHessian) /= effMass;
        }
    }
    return;
}
