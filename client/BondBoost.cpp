#include "BondBoost.h"
#include <math.h>

using namespace helper_functions;

BondBoost::BondBoost(Matter *matt,Parameters *params)
{
    matter = matt;    
    parameters = params;
    nAtoms = matter->numberOfAtoms();

};


BondBoost::~BondBoost()
{
    return;
};


void BondBoost::boost()
{
  printf("hello,let's boost!\n");
} 



