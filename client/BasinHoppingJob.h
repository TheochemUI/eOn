
#ifndef BASINHOPPINGJOB_H
#define BASINHOPPINGJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "Job.h"

class BasinHoppingJob : public Job {
public:
  BasinHoppingJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)), current{new Matter(params)},
        trial{new Matter(params)}, fcalls{0} {}
  ~BasinHoppingJob(void) {
    delete current;
    delete trial;

    for (unsigned int i = 0; i < uniqueStructures.size(); i++) {
      delete uniqueStructures[i];
    }
  };

        std::vector<std::string> run(void) override;
    private:
        VectorXd calculateDistanceFromCenter(Matter *matter);
        AtomMatrix displaceRandom(double maxDisplacement);
        void randomSwap(Matter *matter);
        Matter *current;
        Matter *trial; // initial configuration
        vector<long> getElements(Matter *matter);
        std::vector<std::string> returnFiles;
        int jump_count; // count of jump moves
        int disp_count; // count of displacement moves
        int swap_count; // count of swap moves
        int fcalls;

        std::vector<Matter *> uniqueStructures;
        std::vector<double>   uniqueEnergies;
};

#endif
