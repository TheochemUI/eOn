
/** @file
      Class Lanczos.      
      @author Jean Claude C. Berthet
      @date 2007-2008
      University of Iceland
      */
#define NDEBUG
#ifdef NDEBUG
      #define HAVE_INLINE
      #define GSL_RANGE_CHECK_OFF
#else
      #warning NDEBUG not defined. Compute with NDEBUG defined for release.
      #define BZ_DEBUG
#endif
#include <cfloat>
#include <cassert>
#include <iostream>

#include <blitz/array.h>
#include <random/normal.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_eigen.h>
#include "gsl_eigen.h"
#include "lanczos.hpp"

using namespace blitz;
using namespace gradient_scanning;
//-------------------------------------------------Declaration variable and functions-----------------------------
namespace {
/// Parameters to be used with blitz::Array @{
blitz::firstIndex I;
blitz::secondIndex J;
blitz::thirdIndex K;
blitz::Range const All=Range::all();
/// @}

/** Diagonalise a tridiagonal matrix.
      @pre Size of @a diagonal @em n and size of @a subDiagonal @em n-1.
      eigenvectors are in columns (e.g. first eigenvector: @em eigenvectors(All,0);)
      */
void diagonalise(Array<double, 1> const & diagonal, Array<double, 1> const & subDiagonal, Array<double, 1> & eigenvalues, Array<double, 2> & eigenvectors)
{
      int const n=diagonal.size();
      assert(n == subDiagonal.size()+1);
      gsl_vector const diag=fastWrapToGsl(diagonal);
      gsl_vector const sub_diag=fastWrapToGsl(subDiagonal);      
      
      gsl_vector eval=wrapToGsl(eigenvalues);
      gsl_matrix_view evec=
            gsl_matrix_view_array_with_tda(
                  eigenvectors.data(),
                  eigenvectors.rows(),
                  eigenvectors.cols(),
                  eigenvectors.stride(0)
            ); 
      gsl_eigen_tridiagv_workspace * const w=gsl_eigen_tridiagv_alloc(n);
      gsl_eigen_tridiagv(&diag, &sub_diag, &eval, &evec.matrix, w);
      gsl_eigen_tridiagv_free(w);
      gsl_eigen_symmv_sort(&eval, &evec.matrix, GSL_EIGEN_SORT_VAL_ASC);
}

BZ_DECLARE_FUNCTION_RET(isnan, bool)

} // end anonymous namespace

//----------------------------------------------Definition functions in namespace lanczos------------------------
/** @class gradient_scanning::Lanczos
      @brief Compute the lowest eigenvalue and eigenvector.
      @section hessian Hessian Minimum Mode
            Compute the hessian lowest eigenvalue from function computing gradient by using finite differences.
            @see minimumMode(), gradient_scanning#GradientFunction.
      @section reference Reference
      <em> Comparison of methods fo finding the saddle points without knowledge of the final sate.</em> R.A. Olsen, G.J. Kroes et al. J. Chem. Physics, (2004), Vol. @b 121, n 20.
      */

Lanczos::Lanczos()
{
      initialiseOperatingMembers();
      initialiseParametersMembers();
      initialiseResultsMembers();
}

Lanczos::Lanczos(Lanczos const & lanczos)
{
      operator=(lanczos);
}

/// Delete data about past Lanczos Iterative Methods performed.
void Lanczos::clear()
{
      initialiseOperatingMembers();
      initialiseResultsMembers();
}

/// Convergence reached by last Lanczos.
double Lanczos::convergence() const
{
      return convergence_;
}

/** Parameter controling when the iteration process stops.
      The iteration process stops when @f$ \left| \lambda_i-\lambda_{i-1} \right| < tolerance\_ @f$ where @f$ \lambda @f$ is the smallest eigenvalue.
      @see getMaximumIterations() and setMaximumIterations()
      */
/// @{
double Lanczos::getConvergenceLimit() const
{
      return convergenceLimit_;
}

void Lanczos::setConvergenceLimit(double convergenceLimit)
{
      convergenceLimit_=convergenceLimit;
}
/// @}

/// Lowest eigenvalue for last call to minimumMode().
double Lanczos::eigenvalue() const
{
      return eigenvalue_;
}

/// Lowest mode's eigenvector for last call to minimumMode(). 
blitz::Array<double, 1> Lanczos::eigenvector() const
{
      return eigenvector_;
}

/** Step to calculate derivatives.
      To calculate the second derivatives of the gradient, a finite difference approximation is used: @f$ \mathbf g'=\frac{\mathbf g(\mathbf x+\Delta \mathbf x)-\mathbf g(\mathbf x)}{|\Delta \mathbf x|} @f$. The functions getFiniteDifference() and setFiniteDifference() allows to read and set the step @f$ |\Delta \mathbf x| @f$. 
      @{ */
double Lanczos::getFiniteDifference() const
{
      return finiteDifference_;
}

void Lanczos::setFiniteDifference(double difference)
{
      finiteDifference_=difference;
}
/// @}

/* Generate random vector from Gaussian distribution.
      Average is 0 and Stdv 1.
      @param[out]  vector      Vector to randomise.
      */
void Lanczos::gaussian(blitz::Array<double, 1> & vector)
{
      for (Array<double, 1>::iterator i=vector.begin(); i != vector.end(); i++)
            *i=random::normal.random();
}

/** Tell how to initialise the Lanczos iterative algorithm.
       Lanczos needs to initial vector to start the iterative process.
      The convergence of the process is tested by comparing the lowest eigenvalues found by two successive iterations.
      #Initial controls how the initial vector is chosen and to what value is compared the first eigenvalue.
      */
/// @{
Lanczos::Initial Lanczos::getInitial() const
{
      return initial_;
}

void Lanczos::setInitial(Initial initial)
{
      initial_=initial;
}
/// @}

/* Random isotropic vector of norm one.
      @param[out]  vector      Vector to randomise.
      */
void Lanczos::isotropic(blitz::Array<double, 1> & vector)
{
      double norm;
      do {
            gaussian(vector);
            norm=sum(vector*vector) ;
      } while (norm == 0.0);
      vector/=sqrt(norm);
}

/// Number of iterations used by last Lanczos.
int Lanczos::iteration() const
{
      return iteration_;
}

/** Maximum number of iterations.
      @{ */
int Lanczos::getIterationLimit() const
{
      return iterationLimit_;
}

void Lanczos::setIterationLimit(int iterationLimit)
{
      iterationLimit_=iterationLimit;
}
/// @}

/** Hessian lowest eigenvalue and corresponding eigenvector.
      @param[in] function Function to comput e the gradient.
      @param[in] coordinates Coordinates at which to compute the minimum mode.
      @param[in,out] eigenvalue  Output lowest eigenvalue. For input depends on getInitial().
      @param[in,out] eigenvector Output eigenvector of the lowest eigenvalue. For input depends on getInitial().
      @param[out] gradient Gradient at @a coordinates.
      @return     The status (see status()).
      @note Input values in @a eigenvalue and @a eigenvector are used only if getInitial() is #USER.\n
            When returned status is #UNCOMPLETE, the function returns with the best estimation  it has of the lowest eigenvalue and eigenvector.
            When the returned status is #UNCOMPLETE because the initial vector was an eigenvector, the function returns the eigenvector (i.e. initial vector) and its eigenvalue.
      @pre  If getInitial() is #USER the norm of @a eigenvector must not be zero.
      */
Lanczos::Status Lanczos::minimumMode(
            GradientFunction function,
            blitz::Array<double, 1> const & coordinates,
            double & eigenvalue,
            blitz::Array<double, 1> & eigenvector,
            blitz::Array<double, 1> & gradient
) {
      GradientObject object(function);
      return minimumMode(object, coordinates, eigenvalue, eigenvector, gradient);
}

/** Hessian lowest eigenvalue and corresponding eigenvector.
      @see minimumMode().
      */            
Lanczos::Status Lanczos::minimumMode(
      GradientObject & object,
      blitz::Array<double, 1> const & coordinates,
      double & eigenvalue,
      blitz::Array<double, 1> & eigenvector,
      blitz::Array<double, 1> & gradient
) {
      gradientObject_=&object;
      const unsigned ncoord=coordinates.size();
      gradient.resize(ncoord);
      Array<double, 1> coord(ncoord);
      coord=coordinates;
      coordinates_=&coord;
      gradient_=&gradient;
      setStatus(STARTING);
      object.compute(coord, gradient);
      switch (initial_) {
            case RANDOM:
                  eigenvalue=numeric_limits<double>::infinity();
                  eigenvector.resize(ncoord);
                  isotropic(eigenvector);
                  break;
            case GRADIENT:
                  eigenvalue=numeric_limits<double>::infinity();
                  eigenvector.resize(ncoord);
                  eigenvector=gradient;
            case USER:
                  assert(sum(eigenvector*eigenvector) != 0.0);
                  break;
            case PREVIOUS:
                  eigenvalue=eigenvalue_;
                  eigenvector.resize(ncoord);
                  swap(eigenvector_, eigenvector);
                  break;
            default:
                  std::cerr << __FILE__ << ':' << __FUNCTION__ << ':' << __LINE__ << std::endl;
                  exit(EXIT_FAILURE);                  
      };
      if (eigenvector.size() != ncoord) {
            eigenvector.resize(ncoord);
            isotropic(eigenvector);
            eigenvalue=numeric_limits<double>::infinity();
      };
      Status status=lanczos(&Lanczos::finiteH, eigenvalue, eigenvector);
      gradientObject_=0;
      coordinates_=0;
      gradient_=0;
      if (report_)
            report();
      assert(eigenvector.size() == ncoord);
      return status;
}

/// Tell if a final report is to be printed.
bool Lanczos::getReport() const
{
      return report_;
}

/// Tell if a final report is to be printed.
void Lanczos::setReport(bool turn)
{
      report_=turn;
}

/** Status of the Lanczos process.
      When the lanczos is not running, status indicates either NOT_RUNNING (before any call to minimumMode()) or why the last minimumMode() execution has ended (e.g. CONVERGED). 
      Function status() may also be used while minimumMode() is running, during a call back to the function computing the gradient
      (see parameter @a gradient in
            @link
                  minimumMode(GradientFunction, blitz::Array<double, 1> const &, double &, blitz::Array<double, 1> &, blitz::Array<double, 1> &)
                  minimumMode(GradientFunction, ... )
            @endlink
            and gradient in
            @link
                  minimumMode(GradientObject&, blitz::Array<double, 1> const &, double &, blitz::Array<double, 1> &, blitz::Array<double, 1> &)
                  minimumMode(GradientObject&, ... )
            @endlink
      ).
      In this case status() indicates where the minimisation is currently.
      @see #Status.
      */
Lanczos::Status Lanczos::status() const
{
      return status_;
}

/** Warning when a certain event occurs.
      Here you may specify a certain number of events, you want to be warned about when they occurs.
      The default value is #NOTHING which means that no warning will occur.
      Events may be combined with operator @c bitor.
      Example:
      @code
            warnWhen(EXCEED_ITERATIONS | EIGENVECTOR);
      @endcode
      @see status().
      @{ */
Lanczos::Status Lanczos::getWarnWhen() const
{
      return warnWhen_;
}

void Lanczos::setWarnWhen(Status status)
{
      warnWhen_ = status;
}

void Lanczos::setWarnWhen(int status)
{
      warnWhen_ = static_cast<Status>(status);
}

// @}
/** Copy parameters.
      Copy parameters controling how the Lanczos Iterative method is performed. Note that it does does copy the results displayed by @a lanczos. In other word, results produced by iteration(), eigenvalue(), etc ... will remain different.
      */
Lanczos const & Lanczos::operator=(Lanczos const & lanczos)
{
      initialiseResultsMembers();
      initialiseOperatingMembers();
      convergenceLimit_=lanczos.convergenceLimit_;
      iterationLimit_=lanczos.iterationLimit_;
      finiteDifference_=lanczos.finiteDifference_;
      initial_=lanczos.initial_;
      report_=lanczos.report_;
      warnWhen_=lanczos.warnWhen_;
      return *this;
}

/** Compute product H.q.
      Calling the finite difference @f$ \delta @f$. The product @f$ \mathrm H. \mathbf q @f$ is calculated from the gradient:
       @f$  \mathrm H. \mathbf q =\frac{ \mathbf G(\mathbf x+\delta \mathbf q)-\mathbf G(\mathbf x)}{\delta}  @f$.
       The gradient is computed by function compute();
       @see #finiteDifference_, realH(), and compute().
       @ingroup finiteH
       */
blitz::Array<double, 1> Lanczos::finiteH(blitz::Array<double, 1> const & q)
{
      const int N=q.size();
      Array<double, 1> xk(N), gk(N), result(N);
      xk=*coordinates_+q*getFiniteDifference();
      gradientObject_->compute(xk, gk);
      result=(gk - *gradient_)/getFiniteDifference();
      return result;
}

/// Clear and (re-)initialise members used while running the lanczos iterative method.
void Lanczos::initialiseOperatingMembers()
{
      gradientObject_=0;
      hessian_.free();
      coordinates_=0;
      gradient_=0;
      eigenvalue_=numeric_limits<double>::infinity();
}

/// Clear and (re-)initialise members which controls how the lanczos iterative method is performed.
void Lanczos::initialiseParametersMembers()
{
      convergenceLimit_=0;
      iterationLimit_=20;
      finiteDifference_=1e-5;
      initial_=RANDOM;
      report_=false;
      warnWhen_=NOTHING;
}

/// Clear and (re-)initialise members which contains data about the last lanczos iterative method executed.
void Lanczos::initialiseResultsMembers()
{
      iteration_=0;
      convergence_=0.0;
      eigenvalue_=0.0;
      eigenvector_.free();
      status_=NOTHING;
}

/** Lanczos Iterative Method.
      The function tridiagonalise a symmetric matrix using Lanczos Iterative Method.
      An approximation of the solution may be provided by the user through @a lowestEigenvalue and @a eigenvector.
      If the user does not wish to provide an approximation @a lowestEigenvalue must be set to INFINITY (to make first convergence test fail) and eigenvector must be set to a non-zero random vector.
      The function returns with the latest estimation of @a lowestEigenvalue and @a eigenvector 
      @param[in] H Function computing the product @f$ H.\mathbf{q} @f$.
      @param[in, out] lowestEigenvalue Lowest eigenvalue. 
      @param[in, out] eigenvector Eigenvector corresponding to lowest eigenvalue.
      @return     Why the iterative process stops.
      @see Status minimumMode().
      @note The returned values are copied into #eigenvalue_ and #eigenvector_.
      */
Lanczos::Status Lanczos::lanczos(
            blitz::Array<double, 1> (Lanczos::*H)(const blitz::Array<double, 1>&),
            double & lowestEigenvalue,
            blitz::Array<double, 1> & eigenvector
) {
      int const ncoord=eigenvector.size();
      assert(ncoord > 0);
      assert(iterationLimit_ > 0);
      // Equations and notations are based on the description of the algorithm in RA Olsen, GJ Kroes et al, J. chem. phys. vol 121, no 20 p 9776
      int maxIteration;
      if (iterationLimit_ >= ncoord)
            maxIteration=ncoord;
      else
            maxIteration=iterationLimit_;
      Array<double, 1> eval(maxIteration);// eigenvalues
      Array<double, 2> evec(maxIteration, maxIteration);// eigenvectors
      Array<double, 1> alpha(maxIteration+1), beta(maxIteration+1);// diagonal and sub-diagonal
      Array<double, 1> u(ncoord), r(ncoord);
      Array<double, 2> q(ncoord, maxIteration+1);// From tridiagonal matrix space to canonical
      q(All, 0)=0.0;
      Array<double, 2> vT;// tridiagonal matrix  eigenvectors
      Array<double, 1> a, b, lambda;
      r=eigenvector;
      beta(0)=sqrt(sum(r*r));
      assert(beta(0) > 0.0);
      alpha(0)=0.0;
      setStatus(RUNNING);
      double oldEigenvalue=lowestEigenvalue;
      double convergence=0;
      int k=0;
      Status status=RUNNING;
      while (status == RUNNING) {
            k++;
            q(All, k)=r/beta(k-1);
            u=(this->*H)( q(All, k) );
            r = u - q(All, k-1)*beta(k-1);
            alpha(k)=sum(q(All, k)*r);
            r -= q(All, k)*alpha(k);
            beta(k)=sqrt(sum(r*r));
            if (beta(k) <= 1e-10*fabs(alpha(k))) {
                  -- k;
                  status=UNCOMPLETE;           
                  // if  beta is zero, then wrong vector, so <=
                  break;
            };
            Range range(0, k-1);
            lambda.reference(eval(range));
            vT.reference(evec(range, range));
            if (k > 1) {
                  a.reference(alpha(Range(1, k)));
                  b.reference(beta(Range(1, k - 1)));
                  diagonalise(a, b, lambda, vT);
            }
            else {
                  lambda(0)=alpha(1);
                  vT(0)=1.0;
            };
            convergence=fabs(eval(0) - oldEigenvalue);
            oldEigenvalue=eval(0);
            if (k >= maxIteration) 
                  status=EXCEED_ITERATION;
            if (convergence < convergenceLimit_)
                  status=CONVERGED;
            if (k == ncoord)
                  status=FULL_MATRIX;
      };
      setStatus(status);
      // Memorise the final state of the lanczos
      iteration_ = k;
      convergence_=convergence;
      // It may happen that the number of iteration is 0, lanczos exited with EIGENVECTOR. In this case keep data provided by user.
      if (k > 0) {
            lowestEigenvalue=lambda(0);
            Array<double, 2> Q=q(All, Range(1, k));
            Array<double, 1> v=vT(All, 0);
            eigenvector.resize(ncoord);
            eigenvector=sum(Q(I, J)*v(J), J);
      }
      else {
            assert(status == UNCOMPLETE);
            // case where initial vector is an eigenvector, calculate eigenvalue.
            lowestEigenvalue=sum(u*q(All, 1));// q is normal.
      };
      eigenvalue_=lowestEigenvalue;
      eigenvector_=eigenvector;
      return status_;
}

/** Compute product @f$ \mathrm H. \mathbf q @f$.
      compute the product @f$ \mathrm H. \mathbf q @f$ using the matrix stored in #hessian_ .
      @see #hessian_ and finiteH()
      @ingroup realH
      */
blitz::Array<double, 1> Lanczos::realH(const blitz::Array<double, 1>& q)
{
      const int N=q.size();
      Array<double, 1> u(N);
      u=sum(  hessian_(I, J)*q(J), J  );
      return u;
}

///< Print final report.
void Lanczos::report()
{
      cout << "Lanczos::\n";
      string a, b;
      if (convergence() <= getConvergenceLimit()) {
            a="Has converged: ";
            b=" <= ";
      }
      else {
            a="Has not converged: ";
            b=" > ";
      };
      cout << a << convergence() << b << getConvergenceLimit() << '\n';
      if (iteration() >= getIterationLimit()) {
            a="Has exceeded iteration limit: ";
            b=" >= ";
      }
      else {
            a="Has not exceeded iteration limit: ";
            b=" < ";
      };
      cout << a << iteration() << b << getIterationLimit() << '\n';
      if (iteration() == eigenvector_.size() ) {
            a="Has computed full Hesse matrix: ";
            b=" == ";
      }
      else {
            a="Has not computed full Hesse matrix: ";
            b=" < ";
      };
      cout << a << iteration() << b << eigenvector_.size() << '\n';
      cout << "Status is ";
      switch (status_) {
            case NOTHING:
                  cout << "NOTHING\n";
                  break;
            case CONVERGED:
                  cout << "CONVERGED\n";
                  break;
            case EXCEED_ITERATION:
                  cout << "EXCEED_ITERATION\n";
                  break;
            case FULL_MATRIX:
                  cout << "FULL_MATRIX\n";
                  break;
            case UNCOMPLETE:
                  cout << "UNCOMPLETE\n";
                  break;
            case STARTING:
                  cout << "STARTING\n";
                  break;
            case RUNNING:
                  cout << "RUNNING\n";
                  break;
            default:
                  cerr << "Lanczos: status unknown\n";
                  exit(EXIT_FAILURE);
      };
      cout << "The lowest eigenvalue is " << eigenvalue() << '\n';
      cout << flush;
}

/** Set the new status of the Lanczos.
      Warning may be displayed depending on @a status.
      @see #warnWhen_, warnWhen(), ....
      */
void Lanczos::setStatus(Status status)
{
      status_=status;
      if (status_ & warnWhen_ & CONVERGED)
            cerr << "Lanczos::CONVERGED\n";
      if (status_ & warnWhen_ & EXCEED_ITERATION)
            cerr << "Lanczos::EXCEED_ITERATION\n";
      if (status_ & warnWhen_ & FULL_MATRIX)
            cerr << "Lanczos::FULL_MATRIX\n";
      if (status_ & warnWhen_ & UNCOMPLETE)
            cerr << "Lanczos::UNCOMPLETE\n";
      if (status_ & warnWhen_ & RUNNING)
            cerr << "Lanczos::RUNNING\n";
      if (status_ & warnWhen_ & EXCEED_ITERATION)
            cerr << "Lanczos::STARTING\n";
}
