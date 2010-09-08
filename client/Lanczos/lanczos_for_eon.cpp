/** @file
      Interface Lanczos for EON.     
      @author Jean Claude C. Berthet
      @date 2007
      University of Iceland
      */
#include <blitz/array.h>
#include "../Constants.h"
#include "../Parameters.h"
#include "../Matter.h"

#include "lanczos_for_eon.hpp"
#include "tools.hpp"
#include "lanczos.hpp"

using namespace blitz;
namespace gs=gradient_scanning;

namespace {
      class Matter2GradientObject : public gradient_scanning::GradientObject {
      public:
            Matter2GradientObject(Matter & matter);
            void compute(blitz::Array<double, 1>& coordinates, blitz::Array<double, 1>& gradient);
      private:
            Matter2GradientObject();
            Matter* matter_;
      };

      Matter2GradientObject::Matter2GradientObject(Matter & matter)
      {
            matter_=&matter;
      }

      void Matter2GradientObject::compute(blitz::Array<double, 1>& coordinates, blitz::Array<double, 1>& gradient)
      {
            matter_->updateForces(coordinates.data(), 0, gradient.data());
            gradient*=-1;
      }
}// end anonymous

Lanczos::Lanczos(Matter *const, Parameters *parameters) :
      matter_(parameters)
{
      setConvergenceLimit(parameters->getConvergenceLimit_Lanczos());
      setIterationLimit(parameters->getIterationLimit_Lanczos());
      setFiniteDifference(1e-5);//Angstrom
      setInitial(PREVIOUS);
}

void Lanczos::startNewSearchAndCompute(Matter const *matter, double *)
{
      eigenvector_.free();
      moveAndCompute(matter);
}

void Lanczos::moveAndCompute(Matter const *matter)
{
      Array<double, 1> coordinates(3*matter->numberOfFreeAtoms()), gradient;
      matter->getFreePositions(coordinates.data());
      matter_=*matter;
      Matter2GradientObject forcefield(matter_); 
      minimumMode(static_cast<gs::GradientObject&>(forcefield), coordinates, eigenvalue_, eigenvector_, gradient);
}

double Lanczos::returnLowestEigenmode(double *result)
{
      const int n=eigenvector_.size();
      Array<double, 1> _result(result, n, neverDeleteData);
      _result=eigenvector_;
      return eigenvalue_;
}
