
#ifndef PREFACTOR_H
#define PREFACTOR_H

#include "Eigen.h"
#include "Matter.h"

#include "Parameters.h"

namespace Prefactor {
int getPrefactors(Parameters *parameters, Matter *min1, Matter *saddle,
                  Matter *min2, double &pref1, double &pref2);
Vector<int> movedAtoms(Parameters *parameters, Matter *min1, Matter *saddle,
                       Matter *min2);
Vector<int> movedAtomsPct(Parameters *parameters, Matter *min1, Matter *saddle,
                          Matter *min2);
Vector<int> allFreeAtoms(Matter *matter);
VectorType removeZeroFreqs(Parameters *parameters, VectorType freqs);
void logFreqs(VectorType freqs, char *name);
} // namespace Prefactor
#endif
