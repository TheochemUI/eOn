#pragma once
#ifndef GRADIENT_SCANNING_LANCZOS_HPP
#define GRADIENT_SCANNING_LANCZOS_HPP

/** @file
      Class Lanczos.      
      @author Jean Claude C. Berthet
      @date 2007-2008
      University of Iceland
      */

#include <blitz/array.h>
#include "tools.hpp"
namespace gradient_scanning {
      class Climber;
      class Lanczos {
      public:
            /** Status of the Lanczos process.
                  May represent the state of the process both while running or when stopped.
                  @note Status #UNCOMPLETE is not necessarily an error.
                  @see status()
                  */
            enum Status {
                  NOTHING=0, ///< Before anything is run.
                  CONVERGED=1, ///< Lowest eigenvalue has converged (see setTolerance()).
                  EXCEED_ITERATION=2, ///< Maximum iterations reached (see setIterationLimit()).
                  FULL_MATRIX=4,///< The full Matrix has been computed.
                  /** Could not complete iterative process.
                        This happens when the initial vector is normal to some eigenvectors.
                        */
                  UNCOMPLETE=8,
                  RUNNING =32, ///< Building the lanczos matrix.
                  STARTING =64 ///< First call the gradient function.
            };
            Lanczos();
            Lanczos(Lanczos const & lanczos);
            virtual ~Lanczos(){};
            void clear();
            double convergence() const;
            void setConvergenceLimit(double convergenceLimit);
            double getConvergenceLimit() const;
            double eigenvalue() const;
            blitz::Array<double, 1> eigenvector() const;
            void setFiniteDifference(double difference);
            double getFiniteDifference() const;
            static void gaussian(blitz::Array<double, 1> & vector);        
            /** Specify how to start the Lanczos iterative algorithm.
                  @see setInitial() and minimumMode()
                  */
            enum Initial {
                  RANDOM, ///< Start with a random vector and perform at leat two iterations.
                  GRADIENT,///< Start with gradient vector and perform at leat two iterations.
                  USER,///< Use @a eigenvalue and @a eigenvector provided by user.
                  PREVIOUS///< Use eigenvalue and eigenvector returned by last call to minimumMode(). If not available default to #RANDOM.
            };
            Initial getInitial() const;
            void setInitial(Initial initial);
            static void isotropic(blitz::Array<double, 1> & vector);
            int iteration() const;
            int getIterationLimit() const;
            void setIterationLimit(int iterationLimit);
            Status minimumMode(
                  GradientFunction function,
                  blitz::Array<double, 1> const & coordinates,
                  double& eigenvalue,
                  blitz::Array<double, 1> & eigenvector,
                  blitz::Array<double, 1> & gradient
            );
            Status minimumMode(
                  GradientObject & object,
                  blitz::Array<double, 1> const & coordinates,
                  double & eigenvalue,
                  blitz::Array<double, 1> & eigenvector,
                  blitz::Array<double, 1> & gradient
            );
            template <class C>
            Status minimumMode(
                  C & object,
                  void (C::*member)(blitz::Array<double, 1>&, blitz::Array<double, 1>&),
                  blitz::Array<double, 1> const & coordinates,
                  double & eigenvalue,
                  blitz::Array<double, 1> & eigenvector,
                  blitz::Array<double, 1> & gradient
            );
            bool getReport() const;
            void setReport(bool turnon);
            Status status() const;
            Status getWarnWhen() const;
            void setWarnWhen(Status status);
            void setWarnWhen(int status);
            Lanczos const & operator=(Lanczos const & lanczos);
      private:
            blitz::Array<double, 1> finiteH(const blitz::Array<double, 1>& q);
            GradientObject * gradientObject_;///< Used by finiteH to compute gradient.
            void initialiseOperatingMembers();
            void initialiseParametersMembers();
            void initialiseResultsMembers();
            Status lanczos(
                  blitz::Array<double, 1> (Lanczos::*H)(const blitz::Array<double, 1>&),
                  double & lowestEigenvalue,
                  blitz::Array<double, 1> & eigenvector
            );
            blitz::Array<double, 1> realH(const blitz::Array<double, 1>& q);
            ///< Symmetric matrix used by realH() to compute the product @f$ \mathrm H. \mathbf q @f$ .
            blitz::Array<double, 2> hessian_;
            bool report_;///< Tell to print a report at the end of lanczos.
            void report();
            void setStatus(Status status);
            Status status_;///< Status of the lanczos process. @see status() and #Status
            Status warnWhen_;///< Tell what warnings to display when executing setStatus().
            double convergence_;///< Convergence reached by last Lanczos computation.
            /** Parameter controlling when the iteration process stops.
                  The iteration process stops when @f$ \left| \lambda_i - \lambda_{i-1} \right| < tolerance\_ @f$ where @f$ \lambda @f$ is the smallest eigenvalue.
                  @see #iterationLimit_.*/
            double convergenceLimit_;
            /** Back up the coordinates and forces. For finiteH().*/
            blitz::Array<double, 1> * coordinates_, * gradient_;
            /// Lowest eigenvalue for last search.
            double eigenvalue_;
            /** Last lowest mode eigenvector. 
                  This vector is also the starting vector for the next lanczos calculation. It is renamed @a r in function lanczos().
                  */
            blitz::Array<double, 1> eigenvector_;
            double finiteDifference_;///< Step for derivatives. @see setFiniteDifference()
            Initial initial_;///< @copydoc getInitial().
            int iteration_;///< Number of iterations used by the last Lanczos computation.
            /// Maximum number of iterations (i.e. maximum size of the tridiagonal matrix built by Lanczos). @see setIterationLimit()
            int iterationLimit_;
            friend class Climber;
      };
      template <class C>
      Lanczos::Status Lanczos::minimumMode(
            C & object,
            void (C::*member)(blitz::Array<double, 1>&, blitz::Array<double, 1>&),
            blitz::Array<double, 1> const & coordinates,
            double & eigenvalue,
            blitz::Array<double, 1> & eigenvector,
            blitz::Array<double, 1> & gradient
      ) {
            GradientTemplate<C> gt(object, member);
            return minimumMode(gt, coordinates, eigenvalue, eigenvector, gradient);
      }
}

inline
gradient_scanning::Lanczos::Status operator|=(gradient_scanning::Lanczos::Status a, gradient_scanning::Lanczos::Status b)
{return reinterpret_cast<gradient_scanning::Lanczos::Status&>(reinterpret_cast<int&>(a) |= reinterpret_cast<int&>(b));}

#endif
