#include"Tersoff.h"

Tersoff::Tersoff(void){
    return;
}


void Tersoff::initialize(void){
    return;
}

void Tersoff::cleanMemory(void){
    return;
}
// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void Tersoff::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box){
    tersoff_(&N, R, F, U, &box[0], &box[1], &box[2]);    
    return;
}
