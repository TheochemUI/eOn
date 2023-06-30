#ifndef DIMER_H
#define DIMER_H

#include "HelperFunctions.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"

// dimer method to find the lowest curvature mode
class Dimer : public LowestEigenmode {

public:
  Dimer(std::shared_ptr<Matter> matter, std::shared_ptr<Parameters> params,
        std::shared_ptr<Potential> pot);
  ~Dimer() = default;

  void initialize(Matter *matter, AtomMatrix); // initialize the dimer
  void compute(std::shared_ptr<Matter> matter,
               AtomMatrix initialDirection); // compute the lowest eigenmode
  double getEigenvalue();                    // return the current eigenvalue
  AtomMatrix getEigenvector();               // return the current eigenvector

private:
  std::shared_ptr<spdlog::logger> log;
  std::shared_ptr<Matter> matterCenter; // center of the dimer
  std::shared_ptr<Matter> matterDimer;  // one configuration of the dimer
  AtomMatrix direction;                 // direction along the dimer
  AtomMatrix rotationalPlane; // direction normal to the plane of dimer rotation
  double eigenvalue;          // current curvature along the dimer
  int nAtoms;

  // The rotational plane that is going to be used is determined with the
  // conjugate gradient method
  void determineRotationalPlane(AtomMatrix rotationalForce,
                                AtomMatrix &rotationalForceOld,
                                AtomMatrix rotationalPlaneNormOld,
                                double *lengthRotationalForceOld);

  void
  rotate(double rotationAngle); // rotate the dimer by rotationAngle (radians)
  double calcRotationalForceReturnCurvature(
      AtomMatrix &forceDiffOrthogonalToDimer); // determine the rotational force
                                               // on the dimer
};

#endif
