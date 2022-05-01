#ifndef GPR_AIE_NEBJOB_H_
#define GPR_AIE_NEBJOB_H_

#include "Matter.h"
#include "Parameters.h"
#include "NudgedElasticBand.h"
#include "Job.h"
#include "GPRHelpers.h"
#include "GPRNEB.h"

class GPR_AIE_NEBJob : public Job {

    public:

        GPR_AIE_NEBJob(Parameters *parametersPassed);
        ~GPR_AIE_NEBJob(void);
        std::vector<std::string> run(void);

    private:

        // functions
        void printEndState(int status);
        void saveData(int status, GPRNEB *gpneb);
        void retrainGPR(std::vector<Matter>& newpath);
        void runGPRNEB(GPRNEB& gprneb);
        void runOuterLoop();
        void runRelaxationLoop(std::vector<Matter>& curpath);
        void checkConvergence(double curTrueEnergy);

        // variables
        std::vector<std::string> returnFiles;
        Parameters *eonp;
        size_t fCallsNEB;
        size_t fCallsGPR;
        std::vector<Matter> evaluatedIntermediates, linearMatter, matvec;
        std::vector<GPRMatter> linearPath;
        string reactantFilename, productFilename;
        Matter reactant, product;
        bool stoppedEarly, converged, mustUpdate, isWithin;
        std::shared_ptr<GPRobj> gpf;

};

#endif // GPR_AIE_NEBJOB_H_
