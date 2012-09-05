//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef GLOBALOPTIMIZATIONJOB_H
#define GLOBALOPTIMIZATIONJOB_H

#include "Job.h"
#include "Parameters.h"
#include "Matter.h"

class GlobalOptimizationJob: public Job {
    public:
        GlobalOptimizationJob(Parameters *params);
        ~GlobalOptimizationJob(void);
        void move_step(Matter *);
        void accept_reject_step(Matter *,Matter *);
        void accept_reject_minhopp(Matter *,Matter *);
        void examine_escape(Matter *,Matter *);
        void apply_move_feedback_p1(Matter *);
        void apply_move_feedback_p2(Matter *);
        void mdescape(Matter *);
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
		double etoler;
		double ediff;
		double ekin;
		std::vector<double> earr;
		string escaped;
		string move_type;
		string move_feedback;
		string acc_rej_type;
		string acc_rej_decision;
		string trial_minimum;
};

#endif
