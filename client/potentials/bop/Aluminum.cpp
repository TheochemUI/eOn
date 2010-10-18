#include "Aluminum.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

Aluminum::Aluminum(void)
{
    return;
}

void Aluminum::initialize(void)
{
    potinit_();
    return;
}

void Aluminum::cleanMemory(void)
{
    return;
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void Aluminum::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box)
{
    assert((int)N > 1);
    force_(&N, R, F, U, &box[0], &box[1], &box[2]);    
    return;
}
