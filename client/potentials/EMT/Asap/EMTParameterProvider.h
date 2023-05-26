// EMTParameterProvider.h  --  Abstract class for objects providing EMT
// parameters.

#ifndef _EMTPARAMETERPROVIDER_H
#define _EMTPARAMETERPROVIDER_H

#include "TinyMatrix.h"
#include <string>
using std::string;

// The stucture used to for store EMT parameters.
struct emt_parameters {
  double e0, seq, neq, V0, eta2, kappa, lambda, mass, invmass;
  double gamma1, gamma2;
  // double pairA, pairD;
  double
      lengthscale; // Not an EMT parameter.  1/sqrt(2) of the lattice constant.
  int Z;
  std::string name;
  int index; // An index into various arrays.  Counts the elements in the
             // simulation.
};

// A geometric constant.
// static const double Beta = 1.80939979; // ((16*pi/3)^(1/3))/sqrt(2)
static const double Beta = 1.809; // Preserve same rounding as in ARTwork.

// An EMTParameterProvider provides the EMT parameters to the
// potential.  The GetParameters() method returns the parameters for a
// given element.  If it is called multiple times with the same
// element it must return a pointer to the same emt_parameters struct.
// The Provider owns the parameters, they are deleted when the
// provider is deleted.  After getting _all_ the elements,
// CalcGammaEtc() should be called to calculate the cutoff radius, the
// gammas and other quantities that cannot be calculated before it is
// known which elements are in the simulation.  The gammas and the
// cutoff are determined by the Provider since changing the way they
// are determined corresponds to changing the EMT potential.
class EMTParameterProvider {
public:
  virtual ~EMTParameterProvider(){};
  virtual const emt_parameters *GetParameters(int element) = 0;
  virtual void CalcGammaEtc() = 0;
  virtual double GetCutoffDistance() = 0; // Can it be made element-dependent?
  virtual double GetCutoffSlope() = 0;
  virtual double GetListCutoffDistance() = 0; // Cutoff for neighbor list.
  virtual double
  GetLengthScale() = 0; // The potential delegates this to the provider.
  virtual const TinyDoubleMatrix *GetChi() = 0;
  virtual int GetNumberOfElements() = 0;
};

#endif // ! _EMTPARAMETERPROVIDER_H
