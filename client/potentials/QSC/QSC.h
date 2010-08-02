#include <utility>
#include <map>

#ifndef QSC_STANDALONE
#include "../../PotentialsInterface.h"
#endif
    
class QSC
#ifndef QSC_STANDALONE
: public PotentialsInterface
#endif
{    
    public:
        QSC(void);
        // To satify interface
        void initialize();
        void cleanMemory();
        void force(long N, const double *R, const long *atomicNrs,
                   double *F, double *U, const double *box);
    private:
        bool init;
        struct qsc_parameters {
            int Z;
            double n;
            double m;
            double epsilon;
            double c;
            double a;
        };
        static const qsc_parameters qsc_element_params[];
        
        qsc_parameters get_qsc_parameters(int a, int b);
        double distance(const double *box, const double *R, int i, int j);
        double pair_potential(double r, double a, double n);
};
