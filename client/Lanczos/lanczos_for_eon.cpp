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
#include "oldlanczos.hpp"

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
        int const n=coordinates.size();
        for (int i=0; i<n; ++i) {
            matter_->setPosition(i/3, i%3, coordinates(i));
        };
        Matrix<double, Eigen::Dynamic, 3> forces=matter_->getForces();
        for (int i=0; i<n; ++i) {
            coordinates(i)=matter_->getPosition(i/3, i%3);
            gradient(i)= -forces(i/3, i%3);
        };
    }
} // end anonymous

OldLanczos::OldLanczos(Matter *const, Parameters *parameters) : matter_(parameters)
{
    setConvergenceLimit(parameters->lanczosTolerance);
    setIterationLimit(parameters->lanczosMaxIterations);
    setFiniteDifference(parameters->lanczosFiniteDiff); //Angstrom
    setInitial(USER);
    assert(getIterationLimit() > 0);
}

void OldLanczos::compute(Matter const *matter, AtomMatrix direction)
{
    int const nfree=matter->numberOfFreeAtoms();
    eigenvector_.resize(nfree*3);
    int i, j;
    for (i=0,j=0;i<nfree;i++) {
        if (!matter->getFixed(i)) {
            eigenvector_(j) = direction(i, 0);
            eigenvector_(j+1) = direction(i, 1);
            eigenvector_(j+2) = direction(i, 2);
            j+=3;
        }
    }
    Array<double, 1> coordinates(3*nfree), gradient;
    double * data=coordinates.data();
    j=0;
    for (i=0; i<nfree; ++i) {
        if (not matter->getFixed(i)) {
            for (int a=0; a<3; ++a) {    
                data[j]=matter->getPosition(i, a);
                ++j;
            };
        };
    };
    matter_=*matter;
    Matter2GradientObject forcefield(matter_); 
    minimumMode(static_cast<gs::GradientObject&>(forcefield), coordinates, eigenvalue_, eigenvector_, gradient);
}

double OldLanczos::getEigenvalue()
{
    return eigenvalue_;
}

Matrix<double, Eigen::Dynamic, 3> OldLanczos::getEigenvector()
{
    int const n=matter_.numberOfAtoms();
    Matrix<double, Eigen::Dynamic, 3> result(n, 3);
    int j=0;
    for (int i=0; i < n; ++i) {
        if (not matter_.getFixed(i)) {
            for (int a=0; a < 3; ++a) {
                result(i, a)=eigenvector_(j);
                ++j;
            };
        }
        else {
            for (int a=0; a < 3; ++a) {
                result(i, a)=0.0;
            };
        };
    };
    return result;
}

void OldLanczos::setEigenvector(Matrix<double, Eigen::Dynamic, 3> const eigenvector)
{
    int const n=matter_.numberOfAtoms();
    int j=0;
    for (int i=0; i<n; ++i) {
        if (not matter_.getFixed(i)) {
            for (int a=0; a<3; ++a) {    
                eigenvector_(j)=eigenvector(i, a);
                ++j;
            };
        };
    };
}
