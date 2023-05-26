/** @file
Parameters for potential ZhuPhilpott.
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
@see zhu_philpott.hpp
*/
#include "zhu_philpott_parameters.hpp"
using namespace forcefields::zhu_philpott_parameters;

double const Standard::sigmaO_ = 2.86;                     // Angstrom
double const Standard::epsilonO_ = 0.0023734176137013181;  // eV
double const Standard::sigmaH_ = 2.56;                     // Angstrom
double const Standard::epsilonH_ = 0.00087578073518673099; // eV
double const Standard::C10_O_ = 1.28;
double const Standard::C10_H_ = 1.2;
double const Standard::alpha_ = 0.8;
// Calculated using Smith and Kong mixing rules.
double const Standard::sigmaHPt_ = 2.730249677569295;      // Angstrom
double const Standard::sigmaOPt_ = 2.7735458150747108;     // Angstrom
double const Standard::epsilonHPt_ = 0.016217645873043762; // eV
double const Standard::epsilonOPt_ = 0.03387330127291549;  // eV

double const Iceland::sigmaO_ = 2.86;                     // Angstrom
double const Iceland::epsilonO_ = 0.0023734176137013181;  // eV
double const Iceland::sigmaH_ = 2.56;                     // Angstrom
double const Iceland::epsilonH_ = 0.00087578073518673099; // eV
double const Iceland::C10_O_ = 1.28;
double const Iceland::C10_H_ = 1.2;
double const Iceland::alpha_ = 0.8;
// Calculated using Smith and Kong mixing rules.
double const Iceland::sigmaHPt_ = 2.730249677569295;       // Angstrom
double const Iceland::sigmaOPt_ = 2.7735458150747108;      // Angstrom
double const Iceland::epsilonHPt_ = 0.0097305875238262573; // eV
double const Iceland::epsilonOPt_ = 0.020323980763749295;  // eV
