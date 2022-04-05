/*
 * SCGMatterTest.cpp
 *
 *  Created on: 1 April 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"
#include "../potentials/GPRPotential/GPRPotential.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "SCGMatterTest.h"
#include "../subprojects/gprdimer/gpr/ml/GaussianProcessRegression.h"

namespace tests {

SCGMatterTest::SCGMatterTest() {
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

SCGMatterTest::~SCGMatterTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(SCGMatterTest, TestMatter) {
  string confile{"pos.con"};
  Parameters *parameters = new Parameters;
  parameters->potential = "morse_pt";
  parameters->potential = "morse_pt";
  log_init(parameters, (char *)"test.log");
  Matter *matter = new Matter(parameters);
  matter->con2matter("pos.con");
  // Setup the GPR variables
  gpr::OptimizationAlgorithmSettings settings;
  helper_functions::MatterHolder mh;
  mh.mrr = matter;
  // Set up solver
  settings.check_derivative = false;
  settings.lambda_limit = 1e20;
  settings.max_iter = 400;
  settings.tolerance_func = 1e-4;
  settings.tolerance_sol = 1e-4;
  settings.report_level = 0;
  scg.setAlgorithmSettings(settings);
  Eigen::VectorXd x_ind, y, w, x;
  x_ind.resize(4);
  w.resize(3);
  y.resize(4);
  x.resize(3*mh.mrr->numberOfAtoms());
  x = mh.mrr->getPositionsV();
  x_ind(0) = 0;
  x_ind(1) = 1;
  x_ind(2) = 2;
  x_ind(3) = 3;

    w(0) = -18.7861445569979e+000;
    w(1) = -7.02546595438817e+000;
    w(2) = -7.02546595438817e+000;

    y(0) = 0.00000000000000e+000;
    y(1) = 370.498489473903e-006;
    y(2) = 22.2174753484639e-003;
    y(3) = 21.4452456338877e-003;

  scg.optimize(x, x_ind, y, w,
                 &helper_functions::MatterHolder::getEnergyGradient,
                 mh);

  delete parameters;
  delete matter;
}

} /* namespace tests */
