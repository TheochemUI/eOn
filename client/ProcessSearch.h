#include "Matter.h"
#include "Parameters.h"
#include "SaddlePoint.h"
#include "Prefactors.h"
#include "Job.h"

class ProcessSearch : public Job {
    public:
        ProcessSearch(Parameters *params);
        ~ProcessSearch(void);
        void run(int bundleNumber);

    private:
        int  doProcessSearch(void);
        void printEndState(int status);
        void saveData(int status, int bundleNumber);

        Prefactors *prefactors;
        Parameters *parameters;
        SaddlePoint *saddlePoint; 
        Matter *initial;      // initial configuration.
        Matter *saddle;       // configuration used during the saddle point search.
        Matter *displacement; // configuration used during the saddle point search.
        Matter *min1;         // first minimum from the saddle
        Matter *min2;         // second minimum from the saddle

        double barriersValues[2];
        double prefactorsValues[2];
};
