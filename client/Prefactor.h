
#ifndef PREFACTOR_H
#define PREFACTOR_H

#include "Eigen.h"
#include "Matter.h"

#include "Prefactor.h"
#include "Parameters.h"

namespace Prefactor
{
    const std::string RATE_HTST =            "htst"s;
    const std::string RATE_QQHTST =          "qqhtst"s;
    const std::string FILTER_CUTOFF =        "cutoff"s;
    const std::string FILTER_FRACTION =      "fraction"s;

    int getPrefactors(Parameters *parameters, Matter *min1, Matter *saddle, Matter *min2, double &pref1, double &pref2);
    VectorXi movedAtoms(Parameters* parameters, Matter *min1, Matter *saddle, Matter *min2);
    VectorXi movedAtomsPct(Parameters* parameters, Matter *min1, Matter *saddle, Matter *min2);
    VectorXi allFreeAtoms(Matter *matter);
    VectorXd removeZeroFreqs(Parameters *parameters, VectorXd freqs);
    void logFreqs(VectorXd freqs, char *name);
}
#endif
