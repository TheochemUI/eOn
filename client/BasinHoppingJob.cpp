
#include <math.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>

#include "BasinHoppingJob.h"
#include "Dynamics.h"
#include "HelperFunctions.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Potential.h"
#include "Log.h"

using namespace std;
using namespace helper_functions;

BasinHoppingJob::BasinHoppingJob(Parameters *params)
{
    parameters = params;
    current = new Matter(parameters);
    trial = new Matter(parameters);
    fcalls = Potential::fcalls;
    atomListsCached = false;
}

BasinHoppingJob::~BasinHoppingJob()
{
    delete current;
    delete trial;

    for (unsigned int i = 0; i < uniqueStructures.size(); i++)
    {
        delete uniqueStructures[i];
    }
}

std::vector<std::string> BasinHoppingJob::run(void)
{
    bool swapMove;
    double swap_accept = 0.0;
    jump_count = 0; // count of jump movies
    swap_count = 0; // count of swap moves
    disp_count = 0; // count of displacement moves
    int consecutive_rejected_trials = 0;
    double totalAccept = 0.0;
    Matter *minTrial = new Matter(parameters);
    Matter *swapTrial = new Matter(parameters);

    string conFilename = getRelevantFile(parameters->conFilename);
    current->con2matter(conFilename);

    // Sanity Check
    vector<long> Elements;
    Elements = getElements(current);
    if (parameters->basinHoppingSwapProbability > 0 && Elements.size() == 1)
    {
        char msg[] = "error: [Basin Hopping] swap move probability "
                     "must be zero if there is only one element type\n";
        log(msg);
        exit(1);
    }

    double randomProb = parameters->basinHoppingInitialRandomStructureProbability;
    if (randomProb > 0.0)
    {
        log("generating random structure with probability %.4f\n", randomProb);
    }
    double u = helper_functions::random();
    if (u < parameters->basinHoppingInitialRandomStructureProbability)
    {
        AtomMatrix randomPositions = current->getPositionsFree();
        for (int i = 0; i < current->numberOfFreeAtoms(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                randomPositions(i, j) = helper_functions::random();
            }
        }
        randomPositions *= current->getCell();
        current->setPositionsFree(randomPositions);

        pushApart(current, parameters->basinHoppingPushApartDistance);
    }

    *trial = *current;
    *minTrial = *current;

    current->relax(true);

    double currentEnergy = current->getPotentialEnergy();
    double minimumEnergy = currentEnergy;

    Matter *minimumEnergyStructure = new Matter(parameters);
    *minimumEnergyStructure = *current;
    int nsteps = parameters->basinHoppingSteps + parameters->basinHoppingQuenchingSteps;
    long totalfc;
    FILE *pFile;
    pFile = fopen("bh.dat", "w");

    log("[Basin Hopping] %4s %12s %12s %12s %4s %5s %5s\n",
        "step", "current", "trial", "global min", "fc", "ar", "md");
    log("[Basin Hopping] %4s %12s %12s %12s %4s %5s %5s\n",
        "----", "-------", "-----", "----------", "--", "--", "--");

    int recentAccept = 0;
    double curDisplacement = parameters->basinHoppingDisplacement;

    // Initialize atom selection lists if needed
    if (parameters->basinHoppingDebug) {
        log("[Basin Hopping] DEBUG: Initializing atom selection\n");
    }
    std::vector<int> initialSelection = resolveDisplacementTargets();
    if (parameters->basinHoppingDebug) {
        log("[Basin Hopping] DEBUG: Atom selection initialized: selectedAtoms.size()=%zu\n", initialSelection.size());
    }

    if (!initialSelection.empty())
    {
        log("[Basin Hopping] Selected atoms: ");
        for (size_t i = 0; i < initialSelection.size(); i++)
        {
            log("%d ", initialSelection[i]);
        }
        log("\n");
    }

    for (int step = 0; step < nsteps; step++)
    {

        // Swap or displace
        if (randomDouble(1.0) < parameters->basinHoppingSwapProbability &&
            step < parameters->basinHoppingSteps)
        {
            *swapTrial = *current;
            randomSwap(swapTrial);
            swapMove = true;
            *minTrial = *swapTrial;
        }
        else
        {
            std::vector<int> targets = resolveDisplacementTargets();
            if (parameters->basinHoppingDebug) log("[Basin Hopping] DEBUG: Step %d resolved displacement targets: %zu atoms\n", step + 1, targets.size());

            AtomMatrix displacement = computeDisplacement(curDisplacement, targets);

            trial->setPositions(current->getPositions() + displacement);
            swapMove = false;
            pushApart(trial, parameters->basinHoppingPushApartDistance);

            // Debug: Log coordinates of selected atoms before relax
            if (!targets.empty() && parameters->basinHoppingDebug)
            {
                log("[Basin Hopping] DEBUG: Coordinates of target atoms BEFORE relax:\n");
                for (size_t i = 0; i < targets.size(); i++)
                {
                    int idx = targets[i];
                    if (idx >= 0 && idx < trial->numberOfAtoms())
                    {
                        log("[Basin Hopping] DEBUG: Atom #%d (idx0=%d): [%.6f, %.6f, %.6f]\n",
                            idx, idx, trial->getPosition(idx, 0), trial->getPosition(idx, 1), trial->getPosition(idx, 2));
                    }
                }
            }

            *minTrial = *trial;
        }

        if (parameters->writeMovies == true)
        {
            trial->matter2con("trials", true);
        }

        Potential::fcalls = 0;
        minTrial->relax(true);
        int minfcalls = Potential::fcalls;

        // Debug: Log coordinates of selected atoms after relax
        std::vector<int> targetsAfter = resolveDisplacementTargets();
        if (!targetsAfter.empty() && parameters->basinHoppingDebug)
        {
            log("[Basin Hopping] DEBUG: Coordinates of target atoms AFTER relax:\n");
            for (size_t i = 0; i < targetsAfter.size(); i++)
            {
                int idx = targetsAfter[i];
                if (idx >= 0 && idx < minTrial->numberOfAtoms())
                {
                    log("[Basin Hopping] DEBUG: Atom #%d (idx0=%d): [%.6f, %.6f, %.6f]\n",
                        idx, idx, minTrial->getPosition(idx, 0), minTrial->getPosition(idx, 1), minTrial->getPosition(idx, 2));
                }
            }
        }

        double deltaE = minTrial->getPotentialEnergy() - currentEnergy;
        double p = 0.0;
        if (step >= parameters->basinHoppingSteps)
        {
            if (deltaE <= 0.0)
            {
                p = 1.0;
            }
        }
        else
        {
            if (deltaE <= 0.0)
            {
                p = 1.0;
            }
            else
            {
                p = exp(-deltaE / (parameters->temperature * 8.6173324e-5));
            }
        }

        bool accepted = false;
        if (randomDouble(1.0) < p)
        {
            accepted = true;
            if (parameters->basinHoppingSignificantStructure)
            {
                *current = *minTrial;
            }
            else
            {
                *current = *trial;
            }
            if (swapMove)
            {
                swap_accept += 1;
            }
            if (step < parameters->basinHoppingSteps)
            {
                totalAccept += 1;
                recentAccept += 1;
            }

            currentEnergy = minTrial->getPotentialEnergy();

            if (currentEnergy < minimumEnergy)
            {
                minimumEnergy = currentEnergy;
                *minimumEnergyStructure = *minTrial;
                minimumEnergyStructure->matter2con("min.con");
            }

            consecutive_rejected_trials = 0;

            if (parameters->basinHoppingWriteUnique)
            {
                bool newStructure = true;
                for (unsigned int i = 0; i < uniqueEnergies.size(); i++)
                {
                    // if minTrial has a different energy or a different structure
                    // it is new, otherwise it is old
                    if (fabs(currentEnergy - uniqueEnergies[i]) <
                        parameters->energyDifference)
                    {
                        if (current->compare(uniqueStructures[i],
                                             parameters->indistinguishableAtoms) == true)
                        {
                            newStructure = false;
                        }
                    }
                }

                if (newStructure)
                {
                    uniqueEnergies.push_back(currentEnergy);
                    Matter *currentCopy = new Matter(parameters);
                    *currentCopy = *current;
                    uniqueStructures.push_back(currentCopy);

                    char fname[128];
                    snprintf(fname, 128, "min_%.5i.con", step + 1);
                    current->matter2con(fname);
                    returnFiles.push_back(fname);

                    snprintf(fname, 128, "energy_%.5i.dat", step + 1);
                    returnFiles.push_back(fname);
                    FILE *fh = fopen(fname, "w");
                    fprintf(fh, "%.10e\n", currentEnergy);
                    fclose(fh);
                }
            }
        }
        else
        {
            consecutive_rejected_trials++;
        }

        if (parameters->writeMovies == true)
        {
            minTrial->matter2con("movie", true);
        }

        totalfc = Potential::fcallsTotal;
        char acceptReject[2];
        acceptReject[1] = '\0';
        if (accepted)
        {
            acceptReject[0] = 'A';
        }
        else
        {
            acceptReject[0] = 'R';
        }
        log("[Basin Hopping] %5i %12.3f %12.3f %12.3f %4i %5.3f %5.3f %1s\n",
            step + 1, currentEnergy, minTrial->getPotentialEnergy(), minimumEnergy,
            minfcalls, totalAccept / ((double)step + 1), curDisplacement, acceptReject);
        fprintf(pFile, "%6i %9ld %12.4e %12.4e\n", step + 1, totalfc, currentEnergy,
                minTrial->getPotentialEnergy());

        if (minimumEnergy < parameters->basinHoppingStopEnergy)
        {
            break;
        }

        if (consecutive_rejected_trials == parameters->basinHoppingJumpMax &&
            step < parameters->basinHoppingSteps)
        {
            consecutive_rejected_trials = 0;
            AtomMatrix jump;
            for (int j = 0; j < parameters->basinHoppingJumpSteps; j++)
            {
                jump_count++;
                jump = computeDisplacement(curDisplacement, resolveDisplacementTargets());
                current->setPositions(current->getPositions() + jump);
                if (parameters->basinHoppingSignificantStructure)
                {
                    pushApart(current, parameters->basinHoppingPushApartDistance);
                    current->relax(true);
                }
                currentEnergy = current->getPotentialEnergy();
                if (currentEnergy < minimumEnergy)
                {
                    minimumEnergy = currentEnergy;
                    *minimumEnergyStructure = *current;
                }
            }
        }

        int nadjust = parameters->basinHoppingAdjustPeriod;
        double adjustFraction = parameters->basinHoppingAdjustFraction;
        if ((step + 1) % nadjust == 0 &&
            parameters->basinHoppingAdjustDisplacement == true)
        {
            double recentRatio = ((double)recentAccept) / ((double)nadjust);
            if (recentRatio > parameters->basinHoppingTargetRatio)
            {
                curDisplacement *= 1.0 + adjustFraction;
            }
            else
            {
                curDisplacement *= 1.0 - adjustFraction;
            }

            // log("recentRatio %.3f md: %.3f\n", recentRatio, curDisplacement);
            recentAccept = 0;
        }
    }
    fclose(pFile);

    /* Save Results */

    FILE *fileResults, *fileProduct;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    if (parameters->writeMovies == true)
    {
        std::string movieFilename("movie.xyz");
        returnFiles.push_back(movieFilename);
    }

    fprintf(fileResults, "%d termination_reason\n", 0);
    fprintf(fileResults, "%.6f minimum_energy\n", minimumEnergy);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%.3f acceptance_ratio\n",
            totalAccept / parameters->basinHoppingSteps);
    if (parameters->basinHoppingSwapProbability > 0)
    {
        fprintf(fileResults, "%.3f swap_acceptance_ratio\n",
                swap_accept / double(swap_count));
    }
    fprintf(fileResults, "%ld total_normal_displacement_steps\n",
            disp_count - jump_count - parameters->basinHoppingQuenchingSteps);
    fprintf(fileResults, "%d total_jump_steps\n", jump_count);
    fprintf(fileResults, "%d total_swap_steps\n", swap_count);
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
    fclose(fileResults);

    std::string productFilename("min.con");
    returnFiles.push_back(productFilename);
    fileProduct = fopen(productFilename.c_str(), "wb");
    minimumEnergyStructure->matter2con(fileProduct);
    fclose(fileProduct);

    std::string bhFilename("bh.dat");
    returnFiles.push_back(bhFilename);

    delete minTrial;
    delete minimumEnergyStructure;
    delete swapTrial;
    return returnFiles;
}



void BasinHoppingJob::randomSwap(Matter *matter)
{
    swap_count++;
    vector<long> Elements;
    Elements = getElements(matter);

    long ela;
    long elb;
    long ia = randomInt(0, Elements.size() - 1);
    ela = Elements.at(ia);
    Elements.erase(Elements.begin() + ia);

    long ib = randomInt(0, Elements.size() - 1);
    elb = Elements.at(ib);

    int changera = 0;
    int changerb = 0;

    changera = randomInt(0, matter->numberOfAtoms() - 1);
    while (matter->getAtomicNr(changera) != ela)
    {
        changera = randomInt(0, matter->numberOfAtoms() - 1);
    }

    changerb = randomInt(0, matter->numberOfAtoms() - 1);
    while (matter->getAtomicNr(changerb) != elb)
    {
        changerb = randomInt(0, matter->numberOfAtoms() - 1);
    }

    double posax = matter->getPosition(changera, 0);
    double posay = matter->getPosition(changera, 1);
    double posaz = matter->getPosition(changera, 2);

    matter->setPosition(changera, 0, matter->getPosition(changerb, 0));
    matter->setPosition(changera, 1, matter->getPosition(changerb, 1));
    matter->setPosition(changera, 2, matter->getPosition(changerb, 2));

    matter->setPosition(changerb, 0, posax);
    matter->setPosition(changerb, 1, posay);
    matter->setPosition(changerb, 2, posaz);
}

vector<long> BasinHoppingJob::getElements(Matter *matter)
{
    int allElements[118] = {0};
    vector<long> Elements;

    for (long y = 0; y < matter->numberOfAtoms(); y++)
    {
        if (!matter->getFixed(y))
        {
            int index = matter->getAtomicNr(y);
            allElements[index] = 1;
        }
    }

    for (int i = 0; i < 118; i++)
    {
        if (allElements[i] != 0)
        {
            Elements.push_back(i);
        }
    }

    return Elements;
}

VectorXd BasinHoppingJob::calculateDistanceFromCenter(Matter *matter)
{
    AtomMatrix pos = matter->getPositions();
    Vector3d cen(0, 0, 0);
    int num = matter->numberOfAtoms();

    cen = pos.colwise().sum() / (double)num;

    VectorXd dist(num);

    for (int n = 0; n < num; n++)
    {
        pos.row(n) -= cen;
        dist(n) = pos.row(n).norm();
    }

    return dist;
}

/**
 * Select atoms for displacement based on configuration parameters.
 * Priority: displace_atom_list > displace_type_list > all atoms (default).
 * Caches the atom lists on first call for efficiency.
 *
 * @return Vector of atom indices to use as displacement epicenters
 */
std::vector<int> BasinHoppingJob::resolveDisplacementTargets()
{
    // Cache explicit atom lists on first call
    if (!atomListsCached)
    {
        int nAtoms = current->numberOfAtoms();
        if (parameters->basinHoppingDebug)
        {
            log("[Basin Hopping] DEBUG: Caching atom lists (nAtoms=%d)\n", nAtoms);
        }

        if (!parameters->basinHoppingDisplaceAtomList.empty())
        {
            if (parameters->basinHoppingDebug)
                log("[Basin Hopping] DEBUG: Parsing displace_atom_list: '%s'\n",
                    parameters->basinHoppingDisplaceAtomList.c_str());
            cachedAtomList = helper_functions::parseAtomList(parameters->basinHoppingDisplaceAtomList, nAtoms, parameters->basinHoppingDebug);
            if (parameters->basinHoppingDebug) log("[Basin Hopping] DEBUG: Parsed %zu atoms from displace_atom_list\n", cachedAtomList.size());
        }

        if (!parameters->basinHoppingDisplaceTypeList.empty())
        {
            if (parameters->basinHoppingDebug)
                log("[Basin Hopping] DEBUG: Parsing displace_type_list: '%s'\n",
                    parameters->basinHoppingDisplaceTypeList.c_str());
            cachedTypeList = helper_functions::getAtomsByType(parameters->basinHoppingDisplaceTypeList, current, parameters->basinHoppingDebug);
            if (parameters->basinHoppingDebug) log("[Basin Hopping] DEBUG: Found %zu atoms matching displace_type_list\n", cachedTypeList.size());
        }

        atomListsCached = true;
    }

    std::vector<int> targets;

    // Step 1: Determine Candidate Pool
    // ---------------------------------
    if (!cachedAtomList.empty())
    {
        targets = cachedAtomList;
        if (parameters->basinHoppingDebug) {
            log("[Basin Hopping] DEBUG: Candidate pool: %zu atoms from displace_atom_list\n", targets.size());
        }
    }
    else if (!cachedTypeList.empty())
    {
        targets = cachedTypeList;
        if (parameters->basinHoppingDebug) log("[Basin Hopping] DEBUG: Candidate pool: %zu atoms from displace_type_list\n", targets.size());
    }
    else
    {
        // Default pool: All atoms
        targets.reserve(current->numberOfAtoms());
        for(int i=0; i<current->numberOfAtoms(); ++i) {
            targets.push_back(i);
        }
        if (parameters->basinHoppingDebug) log("[Basin Hopping] DEBUG: Candidate pool: All %d atoms\n", (int)targets.size());
    }

    // Step 2: Apply Single Atom Filter
    // ---------------------------------
    if (parameters->basinHoppingSingleAtomDisplace) 
    {
        if (!targets.empty())
        {
            int randIndex = randomInt(0, targets.size() - 1);
            int selectedAtom = targets[randIndex];
            
            // Log old pool size before clearing
            size_t poolSize = targets.size();
            targets.clear();
            targets.push_back(selectedAtom);
            
            if (parameters->basinHoppingDebug) 
                log("[Basin Hopping] DEBUG: Filtered down to single random atom #%d (from pool of %zu)\n", selectedAtom, poolSize);
        }
        else
        {
            if (parameters->basinHoppingDebug) log("[Basin Hopping] DEBUG: WARNING: Candidate pool is empty, cannot select single atom.\n");
        }
    }

    return targets;
}

/**
 * Apply random displacement to selected atoms only (no neighbors).
 * Only the atoms specified in selectedAtoms will be displaced.
 * Uses the same displacement algorithm and distribution as displaceRandom().
 *
 * @param maxDisplacement Maximum displacement magnitude
 * @param selectedAtoms Vector of atom indices to displace
 * @return Displacement matrix (zero for atoms not selected or fixed)
 */
AtomMatrix BasinHoppingJob::computeDisplacement(double maxDisplacement, const std::vector<int> &targets)
{
    disp_count++;
    AtomMatrix displacement;
    displacement.resize(trial->numberOfAtoms(), 3);
    displacement.setZero();

    if (targets.empty())
    {
        if (parameters->basinHoppingDebug) log("[Basin Hopping] DEBUG: No targets for displacement, returning zero matrix\n");
        return displacement;
    }

    VectorXd distvec = calculateDistanceFromCenter(current);
    double maxDist = distvec.maxCoeff(); // Calculate once based on current structure

    // Apply displacement to selected atoms only
    for (size_t i = 0; i < targets.size(); i++)
    {
        int idx = targets[i];
        if (idx < 0 || idx >= trial->numberOfAtoms() || trial->getFixed(idx))
        {
            continue;
        }

        double dist = distvec(idx);
        double disp = 0.0;

        // Calculate displacement magnitude
        if (parameters->basinHoppingDisplacementAlgorithm == "standard")
        {
            disp = maxDisplacement;
        }
        else if (parameters->basinHoppingDisplacementAlgorithm == "linear")
        {
            double Cs = maxDisplacement / maxDist;
            disp = Cs * dist;
        }
        else if (parameters->basinHoppingDisplacementAlgorithm == "quadratic")
        {
            double Cq = maxDisplacement / (maxDist * maxDist);
            disp = Cq * dist * dist;
        }
        else
        {
            log("Unknown displacement_algorithm\n");
            std::exit(1);
        }

        // Apply random displacement
        for (int j = 0; j < 3; j++)
        {
            if (parameters->basinHoppingDisplacementDistribution == "uniform")
            {
                displacement(idx, j) = randomDouble(2 * disp) - disp;
            }
            else if (parameters->basinHoppingDisplacementDistribution == "gaussian")
            {
                displacement(idx, j) = gaussRandom(0.0, disp);
            }
            else
            {
                log("Unknown displacement_distribution\n");
                std::exit(1);
            }
        }
    }

    // Debug: Log displacement vectors for verification
    if (parameters->basinHoppingDebug) {
        log("[Basin Hopping] DEBUG: Displacement vectors applied:\n");
        for (size_t i = 0; i < targets.size(); i++)
        {
            int idx = targets[i];
            if (idx >= 0 && idx < trial->numberOfAtoms())
            {
                double dx = displacement(idx, 0);
                double dy = displacement(idx, 1);
                double dz = displacement(idx, 2);
                double magnitude = sqrt(dx * dx + dy * dy + dz * dz);
                log("[Basin Hopping] DEBUG: Atom #%d (idx0=%d): displacement=[%.6f, %.6f, %.6f], magnitude=%.6f Å\n",
                    idx, idx, dx, dy, dz, magnitude);
            }
        }
        
        // Detailed Verification: Check all atoms to ensure STRICT adherence to target list
        int nonZeroCount = 0;
        for (int i = 0; i < trial->numberOfAtoms(); i++)
        {
            // O(N*M) verification is fine for debug mode
            bool isTarget = false;
            for(size_t j=0; j<targets.size(); ++j) {
                if(targets[j] == i) { isTarget=true; break; }
            }

            if (!isTarget && !trial->getFixed(i))
            {
                double magnitude = displacement.row(i).norm();
                if (magnitude > 1e-10)
                {
                    log("[Basin Hopping] DEBUG: WARNING: Atom #%d (not in target list) has non-zero displacement: %.10f Å\n", i, magnitude);
                    nonZeroCount++;
                }
            }
        }
        if (nonZeroCount == 0) log("[Basin Hopping] DEBUG: Verification passed: Only target atoms displaced\n");
        else log("[Basin Hopping] DEBUG: ERROR: %d non-target atoms moved!\n", nonZeroCount);
    }

    return displacement;
}
