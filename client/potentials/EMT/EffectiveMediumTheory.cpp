#include "EffectiveMediumTheory.h"
#include <string.h>

// General Functions
EffectiveMediumTheory::EffectiveMediumTheory(void){
	
	// dummy variables 
	AtomsObj = 0;
	EMTObj = 0;
	SuperCellObj = 0;
    numberOfAtoms = 0;

    // should have periodic boundary conditions in all directions
    periodicity[0] = true;      
    periodicity[1] = true;
    periodicity[2] = true;    
}

void EffectiveMediumTheory::cleanMemory(void){
	if(EMTObj != 0){
		delete EMTObj;
        EMTObj = 0;
    }
	if(SuperCellObj != 0){
		delete SuperCellObj;
        SuperCellObj = 0;
    }
	if(AtomsObj != 0){
        delete AtomsObj;
        AtomsObj = 0;
    }
    return;
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void EffectiveMediumTheory::force(long N, const double *R, const int *atomicNrs,
                                  double *F, double *U, const double *box)
{
    int i, j;
    double *pos;
    
    pos = new double[3*N];
        
    for(i=0; i<3*N; i++)
        pos[i] = R[i];

    
	// an atom has been deposited
	if(numberOfAtoms != N)
	{
        numberOfAtoms = N;
        cleanMemory();
        int *atomicNrsTemp;
        atomicNrsTemp = new int[N];
				
		// create new atoms / emt potential with N atoms
		Vec tempBasisX(box[0], 0.0, 0.0);
		Vec tempBasisY(0.0, box[1], 0.0);
		Vec tempBasisZ(0.0, 0.0, box[2]);
		Vec tempBasis[3];
		tempBasis[0] = tempBasisX;
		tempBasis[1] = tempBasisY;
		tempBasis[2] = tempBasisZ;
        
		SuperCellObj = new SuperCell(tempBasis, periodicity);
		AtomsObj = new Atoms((Vec *) pos, N, SuperCellObj);

        for(j=0; j<N; j++)
            atomicNrsTemp[j] = int(atomicNrs[j]);
        
		AtomsObj->SetAtomicNumbers(atomicNrsTemp);

		EMTObj = new EMT(NULL);
		AtomsObj->SetCalculator(EMTObj);        
        delete [] atomicNrsTemp;
	}
	AtomsObj->SetCartesianPositions((Vec *) pos);
	// update the box
	Vec tempBasisX(box[0], 0.0, 0.0);
	Vec tempBasisY(0.0, box[1], 0.0);
	Vec tempBasisZ(0.0, 0.0, box[2]);
	Vec tempBasis[3];
	tempBasis[0] = tempBasisX;
	tempBasis[1] = tempBasisY;
	tempBasis[2] = tempBasisZ;
    AtomsObj->SetUnitCell(tempBasis, true);	
    

	*U = EMTObj->GetPotentialEnergy();
	
	// converts data from EMT to suite EON
	const Vec *tempF = EMTObj->GetCartesianForces();
	memcpy(F, tempF, N*sizeof(Vec));

    delete [] pos;
    return;
}
