#ifndef SADDLESEARCHMETHOD_H
#define SADDLESEARCHMETHOD_H

#include "Parameters.h"
#include "Potential.h"

class SaddleSearchMethod {
protected:
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Parameters> params;

public:
  SaddleSearchMethod(std::shared_ptr<Potential> potPassed,
                     std::shared_ptr<Parameters> paramsPassed)
      : pot{potPassed}, params{paramsPassed} {};
  virtual ~SaddleSearchMethod(){};
  virtual int run() = 0;
  virtual double getEigenvalue() = 0;
  virtual AtomMatrix getEigenvector() = 0;

  int status;
};

#endif
