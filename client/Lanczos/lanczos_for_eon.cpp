/** @file
Interface Lanczos for EON.     
@author Jean Claude C. Berthet
@date 2007
University of Iceland
*/
/** @page Lanczos Lanczos Itertive Method
The Lanczos Iterative method is an alternative to the Dimer method. It has been included in eOn. To compile Lanczos two extra libraries are required GSL (GNU Scientific Library) and Blitz++.
These can be downloaded at: http://www.gnu.org/software/gsl/ and http://www.oonumerics.org/blitz/download/ (or check out at
@verbatim
export CVSROOT=:pserver:anonymous@blitz.cvs.sourceforge.net:/cvsroot/blitz
cvs login
cvs -z3 checkout blitz </cc>
@endverbatim
). Because of the inconvenience, Lanczos is not compiled by default. To include it, the client must be with the environment variable WITH_LANCZOS defined. In bash:
@code
export WITH_LANCZOS=1
@endcode
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
