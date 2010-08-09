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

        struct distance {
            double d[3];
            double r;
        };
        distance **distances;

        int **vlist;
        int  *nlist;
        double cutoff;
        double verlet_skin;
        double *oldR;

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
        double pair_potential(double r, double a, double n);
        void new_vlist(long N, const double *R, const double *box);
        void update_vlist(long N, const double *R, const double *box);
        void calc_distance(const double *box, const double *R, int i, int j, struct distance *d);
};
