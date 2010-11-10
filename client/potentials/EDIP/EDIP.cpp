#include"EDIP.h"

EDIP::EDIP(void){
    return;
}

void EDIP::initialize(void){
    return;
}

void EDIP::cleanMemory(void){
    return;
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void EDIP::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box){
    edip_(&N, R, F, U, &box[0], &box[4], &box[8]);    
    return;
}
