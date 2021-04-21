#ifndef GLOBALOPTIMIZATION_HJOB_H
#define GLOBALOPTIMIZATION_HJOB_H

#include "Job.h"
#include "Parameters.h"
#include "Matter.h"

class GlobalOptimizationJob: public Job {
    public:
        GlobalOptimizationJob(Parameters *params);
        ~GlobalOptimizationJob(void);
        void hoppingStep(long,Matter *,Matter *);
        void decisionStep(Matter *,Matter *);
        void report(Matter *);
        void acceptRejectNPEW(Matter *,Matter *);
        void acceptRejectBoltzmann(Matter *, Matter *);
        void analyze(Matter *,Matter *);
        void examineEscape(Matter *,Matter *);
        void applyMoveFeedbackMD(void);
        void applyDecisionFeedback(void);
        void mdescape(Matter *);
        void randomMove(Matter *);
        void insert(Matter *);
        size_t hunt(double);
        void velopt(Matter *);
        std::vector<std::string> run(void);
		double beta1;
		double beta2;
		double beta3;
		double alpha1;
		double alpha2;
		size_t mdmin;
    private:
        Parameters *parameters;
		size_t nlmin;
		long fcallsMove;
		long fcallsRelax;
		double ediff;
		double ekin;
		bool firstStep;
		std::vector<double> earr;
		string escapeResult;
		//string moveFeedbackMethod;
		//string decisionMethod;
		string decisionResult;
		string hoppingResult;
};

#endif
