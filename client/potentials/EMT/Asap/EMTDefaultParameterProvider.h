// EMTDefaultParameterProvider.h  --  provides the default EMT parameters.

#ifndef _EMTDEFAULTPARAMETERPROVIDER_H
#define _EMTDEFAULTPARAMETERPROVIDER_H

#include "EMTParameterProvider.h"
#include <vector>
using std::vector;

class EMTDefaultParameterProvider : public EMTParameterProvider {
public:
  EMTDefaultParameterProvider();
  virtual ~EMTDefaultParameterProvider();
  virtual const emt_parameters *GetParameters(int element);
  virtual void CalcGammaEtc();
  // Can the cutoff be made element-dependent?
  virtual double GetCutoffDistance() { return cutoff; }
  virtual double GetListCutoffDistance() { return cutoff * listcutofffactor; }
  virtual double GetCutoffSlope() { return cutslope; }
  // The potential delegates GetLengthScale to the provider.
  // What should be done for multiple elements?
  virtual double GetLengthScale() { return params[0]->lengthscale; }
  virtual int GetNumberOfElements() { return params.size(); }
  virtual const TinyDoubleMatrix *GetChi() { return chi; }
  virtual void Debug();

#ifndef SWIG // SWIG cannot parse all of this
protected:
  virtual emt_parameters *GetNewParameters(int element);
  virtual void calc_cutoff();
  virtual void calc_gammas();
  virtual void calc_chi();

  std::vector<emt_parameters *> params;
  TinyDoubleMatrix *chi;
  double maxseq;
  double cutoff;
  double cutslope;
  double listcutofffactor;

  // Two constants: the last shell included in interaction range, and in nb list
  static const int shell0;
  static const int shell1;
#endif
};

#endif // ! _EMTDEFAULTPARAMETERPROVIDER_H
