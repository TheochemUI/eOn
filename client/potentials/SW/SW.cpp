#include"SW.h"

SW::SW(void){
    return;
}

void SW::initialize(void){
    return;
}

void SW::cleanMemory(void){
    return;
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void SW::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box){
    sw_(&N, R, F, U, &box[0], &box[4], &box[8]);    
    return;
}
