#ifndef QSC_STANDALONE
#include "../../Potential.h"
#endif

class QSC
#ifndef QSC_STANDALONE
: public Potential
#endif
{    
    public:
        QSC(void);
        ~QSC(void);
        // To satify interface
        void initialize();
        void initialize(long N, const double *R, const int *atomicNrs,
                        const double *box);
        void cleanMemory();
        void force(long N, const double *R, const int *atomicNrs,
                   double *F, double *U, const double *box);
        void energy(long N, const double *R, const int *atomicNrs,
                    double *U, const double *box);
    private:
        bool init;

        struct distance {
            double d[3];
            double r;
        };
        distance **distances;

        int natomstoclear;
        int **vlist;
        int  *nlist;
        double cutoff;
        double verlet_skin;
        double *oldR;
        double *rho;
        double *sqrtrho;
        double **V;
        double **phi;

        struct qsc_parameters {
            int Z;
            double n;
            double m;
            double epsilon;
            double c;
            double a;
        };
        static const qsc_parameters qsc_element_params[];
        int unique_elements[9];
        int largest_element_num;
        qsc_parameters **qsc_param_cache;
        
        qsc_parameters get_qsc_parameters(int a, int b);
        double pair_potential(double r, double a, double n);
        void new_vlist(long N, const double *R, const double *box);
        void update_vlist(long N, const double *R, const double *box);
        void calc_distance(const double *box, const double *R, int i, 
                           int j, struct distance *d);
        int int_comp(const void* a,const void* b);
};
