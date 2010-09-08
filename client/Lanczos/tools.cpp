/** @file
      Tools to handle functions or object computing the gradient
      @author Jean Claude C. Berthet
      @date 2006-2007
      University of Iceland
      */

/** @namespace gradient_scanning
      @brief ConjugateGradient, SaddlePoint searches.
      Tools to scan a surface. The namespace contains a minimiser: class ConjugateGradient, and a class to search for saddle points: SaddlePoint. The namespace contains also other classes and functions needed by the two latter. \n
      Class Lanczos performs the Lanczos iterative method and is used by SaddlePoint. Files gradient_scanning_tools.h and gradient_scanning_tools.cpp contains basic tools used by all major classes. In particular, it contains types (#GradientFunction) and classes (GradientObject ...) used to manage the functionss computing the gradient.\n
      Tools in this namespace requires library blitz++ which can be downloaded at http://www.oonumerics.org/blitz . Class Lanczos and SaddlePoint requires library GSL (GNU Scientific Library) to diagonalise matrices. GSL can be downloaded at http://www.gnu.org/software/gsl .
      @see ConjugateGradient, SaddlePoint.
      @note Requires library blitz++ and GSL.
      */
      
#include "tools.hpp"
#include <iostream>
using namespace std;
using namespace blitz;
using namespace gradient_scanning;

/** @class gradient_scanning::GradientObject
      @brief To compute gradient.
      To compute gradient from a function (#GradientFunction), or from a member function (GradientTemplate) 
      @see #GradientFunction, GradientTemplate
      */
      
GradientObject::GradientObject(GradientFunction function) :
      function_(function)
{}

/** Compute gradient.
      */
void GradientObject::compute(
      blitz::Array<double, 1>& coordinates, 
      blitz::Array<double, 1>& gradient
) {
      function_(coordinates, gradient);
}

/** @namespace gradient_scanning::random
      @brief random number generator.
      The class is just a wrapper for random number generators defined in Blitz++ (http://www.oonumerics.org/blitz/) library
      @see random/uniform.h and random/normal.h 
      */

/// Seed generator
unsigned int gradient_scanning::random::seedWithTime()
{
      unsigned int seed=time(0);
      return random::seed(seed);
}

unsigned int gradient_scanning::random::seed(unsigned int seed)
{
      normal.seed(seed);
      uniform.seed(seed);
      return seed;
}

ranlib::Normal<double> gradient_scanning::random::normal(0.0, 1.0);
/// Uniform distribution in interval [0;1]
ranlib::Uniform<double> gradient_scanning::random::uniform;

blitz::Array<double, 1> gradient_scanning::wrapToBlitz(gsl_vector * const vgsl)
{
      assert(vgsl->stride == 1);
      blitz::Array<double, 1> vblitz(vgsl->data, vgsl->size, neverDeleteData);
      return vblitz;
}

blitz::Array<double, 2> gradient_scanning::wrapToBlitz(gsl_matrix * const mgsl)
{
      assert(mgsl->tda == mgsl->size1);
      TinyVector<double, 2> size=shape(mgsl->size1, mgsl->size2);
      blitz::Array<double, 2> mblitz(mgsl->data, size, neverDeleteData);
      return mblitz;
}

gsl_vector gradient_scanning::wrapToGsl(blitz::Array<double, 1> & vblitz)
{
      assert(vblitz.isStorageContiguous());
      gsl_vector const vgsl={vblitz.size(), 1, vblitz.data(), 0, 0};
      return vgsl;
}

gsl_matrix gradient_scanning::wrapToGsl(blitz::Array<double, 2> & mblitz)
{
      assert(mblitz.isStorageContiguous());
      gsl_matrix mgsl={mblitz.rows(), mblitz.cols(), mblitz.cols(), mblitz.data(), 0, 0};
      return mgsl;
}

gsl_vector const gradient_scanning::fastWrapToGsl(blitz::Array<double, 1> const & vblitz)
{
      assert(vblitz.isStorageContiguous());
      gsl_vector vgsl={vblitz.size(), 1, const_cast<double*>(vblitz.data()), 0, 0};
      return vgsl;
}

gsl_matrix const gradient_scanning::fastWrapToGsl(blitz::Array<double, 2> const & mblitz)
{
      assert(mblitz.isStorageContiguous());
      gsl_matrix mgsl={mblitz.rows(), mblitz.cols(), mblitz.cols(), const_cast<double*>(mblitz.data()), 0, 0};
      return mgsl;
}
