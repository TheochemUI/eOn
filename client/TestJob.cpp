//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------
#include "TestJob.h"
#include "Matter.h"
#include "Constants.h"
#include "Potential.h"
#include "SaddleSearch.h"

#include <stdlib.h>

TestJob::TestJob(Parameters *params)
{
    tolerance = 0.01;
    parameters = params;
}

TestJob::~TestJob(){ }

void TestJob::run(int bundleNumber)
{
    checkPotentials();
    checkFullSearch();
}

void TestJob::checkFullSearch(void){
    printf("\n---Beginning tests of saddle point search---\n");
    printf("Checks the potential energies of located configurations.\n");
    printf("Reported as OK if within a tolerance of: %f\n",tolerance);

    long status;
    bool ok=1;
    double diffM1, diffM2, diffSP;
    
    Matter *initial;      
    Matter *saddle;
    Matter *displacement; 
    Matter *min1;
    Matter *min2; 
    Matter *matterTemp; 
    SaddlePoint *saddlePoint;     
    
    string reactant_passed("reactant_test.con");
    string displacement_passed("displacement_test.con");
    string mode_passed("mode_test.dat");
    
    parameters->potential = "emt";
    
    initial = new Matter(parameters);
    displacement = new Matter(parameters);
    saddle = new Matter(parameters);
    min1 = new Matter(parameters);
    min2 = new Matter(parameters);
    matterTemp = new Matter(parameters);

    saddle->con2matter(displacement_passed);
    initial->con2matter(reactant_passed);
    *min1 = *min2 = *initial;

    printf("\n---Output for saddle point search start---\n");    
    saddlePoint = new SaddlePoint();
    saddlePoint->initialize(initial, saddle, parameters);
    saddlePoint->loadMode(mode_passed);
    status = saddlePoint->locate();
	printf("---Output for saddle point search end---\n\n");
	
	printf("---Output relax from saddle point search start---\n");    
	// relax from the saddle point located
	
	Matrix<double, Eigen::Dynamic, 3> posSaddle = saddlePoint->getSaddlePositions();
	Matrix<double, Eigen::Dynamic, 3> displacedPos;
	
	*min1 = *saddle;
	//XXX: the distance displaced from the saddle should be a parameter
	displacedPos = posSaddle - saddlePoint->getEigenMode() * 0.2;
	min1->setPositions(displacedPos);
	ConjugateGradients cgMin1(min1, parameters);
	cgMin1.fullRelax();
//	fCallsMin += cgMin1.totalForceCalls;
	
	*min2 = *saddle;
	displacedPos = posSaddle + saddlePoint->getEigenMode() * 0.2;
	min2->setPositions(displacedPos);
	ConjugateGradients cgMin2(min2, parameters);  
	cgMin2.fullRelax();
//	fCallsMin += cgMin2.totalForceCalls;
	
    // If min2 corresponds to initial state swap min1 && min2
    if(!(*initial==*min1) && ((*initial==*min2))){
        *matterTemp = *min1;
        *min1 = *min2;
        *min2 = *matterTemp;
    }
	printf("---Output relax from saddle point search end---\n");
    
    // checking the energies of the obtained configurations
    diffM1 = abs(min1->getPotentialEnergy()-45.737426);
    diffM2 = abs(min2->getPotentialEnergy()-45.737433);
    diffSP = abs(saddle->getPotentialEnergy()-46.284511);
    
    if ((diffM1 < tolerance) and
        (diffM2 < tolerance) and
        (diffSP < tolerance)){
        ok *= 1; 
        printf("OK: Saddle search structural energies\n");
    }
    else{
        if (tolerance < diffSP){
            ok *= 0;
            printf("WARNING: Saddle point not within energy tolerance: %f\n", diffSP);
        }
        if (tolerance < diffM2){
            ok *= 0;
            printf("WARNING: Minimum 2 not within energy tolerance: %f\n", diffM2);
        }
        if (tolerance < diffM1){
            ok *= 0;
            printf("WARNING: Minimum 1 not within energy tolerance: %f\n", diffM2);
        }
    }
    
    // checking the structures of the obtained configurations
    diffM1 = abs((min1->getPositions()).row(384).norm()-19.123375);
    diffM2 = abs((min2->getPositions()).row(384).norm()-19.527995);
    diffSP = abs((saddle->getPositions()).row(384).norm()-19.026709);

    if ((diffM1 < tolerance) and
        (diffM2 < tolerance) and
        (diffSP < tolerance)){
        ok *= 1; 
        printf("OK: Saddle search, adatom positions\n");
    }
    else{
        if (tolerance < diffSP){
            ok *= 0;
            printf("WARNING: Saddle point, adatom not within position tolerance: %f\n", diffSP);
        }
        if (tolerance < diffM2){
            ok *= 0;
            printf("WARNING: Minimum 2, adatom not within position tolerance: %f\n", diffM2);
        }
        if (tolerance < diffM1){
            ok *= 0;
            printf("WARNING: Minimum 1, adatom not within position tolerance: %f\n", diffM1);
        }
    }
    

    if (ok){
        printf("Saddle search tests all good\n");
    }
    else{
        printf("Saddle search tests there were WARNINGS\n");
    }
    printf("SP done\n");
    
    return;
    
}
/*
void TestJob::checkMode(void){

    string reactant_passed("reactant_test.con");
    parameters->potential = 1;  // always LJ in test
    Matter *saddle = new Matter(parameters);
    saddle->con2matter(reactant_passed);
    
    LowestEigenmodeInterface *lowestEigenmode = new Dimer(saddle, parameters);
    
}
*/
void TestJob::checkPotentials(void)
{
    double energyDiff;
    double forceDiff;
    
    printf("\n---Beginning tests of potentials---\n");
    printf("Checks the potential energy and the max force.\n");
    printf("Reported as OK if within a tolerance of: %f\n\n",tolerance);
    
    energyDiff = getEnergyDiff("lj", -1475.984331);
    if (abs(energyDiff) > tolerance){
        printf("WARNING: LJ energy difference: %f\n", energyDiff);
    }else{
        forceDiff = getForceDiff("lj", 2.007213);
        if (abs(forceDiff) > tolerance){
            printf("WARNING: LJ force difference: %f\n", forceDiff);
        }else{
            printf("OK: LJ\n");
        }
    }    
    
    energyDiff = getEnergyDiff("emt", 46.086312);
    if (abs(energyDiff) > tolerance){
        printf("WARNING: EMT energy difference: %f\n", energyDiff);
    }else{
        forceDiff = getForceDiff("emt", 0.357493);
        if (abs(forceDiff) > tolerance){
            printf("WARNING: EMT force difference: %f\n", forceDiff);
        }else{
            printf("OK: EMT\n");
        }
    } 
    
    energyDiff = getEnergyDiff("edip" ,-1033.250950);
    if (abs(energyDiff) > tolerance){
        printf("WARNING: EDIP energy difference: %f\n", energyDiff);
    }else{
        forceDiff = getForceDiff("edip", 7.080115);
        if (abs(forceDiff) > tolerance){
            printf("WARNING: EDIP force difference: %f\n", forceDiff);
        }else{
            printf("OK: EDIP\n");
        }
    } 

    energyDiff = getEnergyDiff("tersoff_si",-1035.809985 );
    if (abs(energyDiff) > tolerance){
        printf("WARNING: Tersoff energy difference: %f\n", energyDiff);
    }else{
        forceDiff = getForceDiff("tersoff_si", 11.145002);
        if (abs(forceDiff) > tolerance){
            printf("WARNING: Tersoff force difference: %f\n", forceDiff);
        }else{
            printf("OK: Tersoff\n");
        }
    } 

    energyDiff = getEnergyDiff("sw_si", -1449.795645);
    if (abs(energyDiff) > tolerance){
        printf("WARNING: SW energy difference: %f\n", energyDiff);
    }else{
        forceDiff = getForceDiff("sw_si", 2.530904);
        if (abs(forceDiff) > tolerance){
            printf("WARNING: SW force difference: %f\n", forceDiff);
        }else{
            printf("OK: SW\n");
        }
    } 

    energyDiff = getEnergyDiff("lenosky_si", -1410.679106);
    if (abs(energyDiff) > tolerance){
        printf("Lenosky energy difference: %f\n", energyDiff);
    }else{
        forceDiff = getForceDiff("lenosky_si", 2.320168);
        if (abs(forceDiff) > tolerance){
            printf("Lenosky force difference: %f\n", forceDiff);
        }else{
            printf("OK: Lenosky\n");
        }
    } 

//    energyDiff = getEnergyDiff(Potential::POT_LJBINARY,0);
//    if (abs(energyDiff) > tolerance){
//        printf("WARNING: LJBinary energy difference: %f\n", energyDiff);
//    }else{
//        printf("OK: LJBinary\n");
//    } 

    energyDiff = getEnergyDiff("eam_al", -1206.825825);
    if (abs(energyDiff) > tolerance){
        printf("WARNING: Aluminum energy difference: %f\n", energyDiff);
    }else{
        forceDiff = getForceDiff("eam_al", 0.000246);
        if (abs(forceDiff) > tolerance){
            printf("WARNING: Aluminum force difference: %f\n", forceDiff);
        }else{
            printf("OK: Aluminum\n");
        }
    } 
    
//    energyDiff = getEnergyDiff(Potential::POT_EAM,0);
//    if (abs(energyDiff) > tolerance){
//        printf("WARNING: EAM energy difference: %f\n", energyDiff);
//    }else{
//        printf("OK: EAM\n");
//    } 

    energyDiff = getEnergyDiff("qsc", -1232.806318);
    if (abs(energyDiff) > tolerance){
        printf("WARNING: QSC energy difference: %f\n", energyDiff);
    }else{
        forceDiff = getForceDiff("qsc", 0.673444);
        if (abs(forceDiff) > tolerance){
            printf("WARNING: QSC force difference: %f\n", forceDiff);
        }else{
            printf("OK: QSC\n");
        }
    } 

//    energyDiff = getEnergyDiff(Potential::POT_ZPICE, 0);
//    if (abs(energyDiff) > tolerance){
//        printf("WARNING: ZPICE energy difference: %f\n", energyDiff);
//    }else{
//        printf("OK: ZPICE\n");
//    } 

    energyDiff = getEnergyDiff("tip4p", 4063.865115);
    if (abs(energyDiff) > tolerance){
        printf("WARNING: TIP4P energy difference: %f\n", energyDiff);
    }else{
        forceDiff = getForceDiff("tip4p", 73.655248);
        if (abs(forceDiff) > tolerance){
            printf("WARNING: TIP4P force difference: %f\n", forceDiff);
        }else{
            printf("OK: TIP4P\n");
        }
    } 

//    energyDiff = getEnergyDiff(Potential::POT_BOPFOX,0);
//    if (abs(energyDiff) > tolerance){
//        printf("WARNING: BOPFOX energy difference: %f\n", energyDiff);
//    }else{
//        printf("OK: BOPFOX\n");
//    } 

//    energyDiff = getEnergyDiff(Potential::POT_BOP,0);
//    if (abs(energyDiff) > tolerance){
//        printf("WARNING: BOP energy difference: %f\n", energyDiff);
//    }else{
//        printf("OK: BOP\n");
//    } 
}

double TestJob::getEnergyDiff(string pot, double refEnergy)
{
    string reactant_passed("reactant_test.con");
    parameters->potential = pot;
    Matter *reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);
//    printf("Energy: %f\n", reactant->getPotentialEnergy());
    return reactant->getPotentialEnergy()-refEnergy;
}


double TestJob::getForceDiff(string pot, double refForce)
{
    string reactant_passed("reactant_test.con");
    parameters->potential = pot;
    Matter *reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);
//    printf("Force: %f\n", reactant->maxForce());
    return reactant->maxForce()-refForce;
}


