#ifndef NEBJOB_H
#define NEBJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "NudgedElasticBand.h"
#include "Job.h"

class NudgedElasticBandJob : public Job {

    public:

        NudgedElasticBandJob(Parameters *parametersPassed);
        ~NudgedElasticBandJob(void);
        std::vector<std::string> run(void);

    private:

        // functions
        void printEndState(int status);
        void saveData(int status, NudgedElasticBand *neb);

        // variables
        std::vector<std::string> returnFiles;
        Parameters *parameters;
        int fCallsNEB;

};

#endif
