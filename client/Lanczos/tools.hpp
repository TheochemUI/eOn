#pragma once
#ifndef GRADIENT_SCANNING_TOOLS_HPP
#define GRADIENT_SCANNING_TOOLS_HPP
#include <boost/format.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <blitz/array.h>
#include <random/uniform.h>
#include <random/normal.h>

/** @file
      Tools to handle functions or object computing the gradient
      @author Jean Claude C. Berthet
      @date 2006-2007
      University of Iceland
      */

namespace blitz {
      template <class C, int N> class Array;
}

namespace gradient_scanning {
      /** Function to compute the gradient.
            The user enters an array of @a coordinates. The function returns an array containing the @a gradient. The parameter @a coordinates is not constant as it may be changed by the function to apply constraints (bond length constraints, periodic boundaries, etc...)
            @param[in,out] coordinates Coordinates. Return coordinates with applied constraints is any.
            @param[out] gradient Gradient at @a coordinates.
            */
      typedef void (*GradientFunction)(blitz::Array<double, 1>& coordinates, blitz::Array<double, 1>& gradient);
      class GradientObject {
      public:
            GradientObject(GradientFunction function);
            virtual ~GradientObject() {};
            virtual void compute(blitz::Array<double, 1> & coordinates, blitz::Array<double, 1> & gradient);
      protected:
            GradientObject() {};
      private:
            GradientObject(GradientObject const &);
            void operator=(GradientObject const &);
            GradientFunction function_;///< Pointer to function to be excuted by compute().
      };
      /// When gradient is computed by member function. 
      template <class T>
      class GradientTemplate : public GradientObject {
      public:
            typedef void (T::*GradientMember)(blitz::Array<double, 1>& coordinates, blitz::Array<double, 1>& gradient);
            GradientTemplate(T & object, GradientMember member);
            void compute(blitz::Array<double, 1>& coordinates, blitz::Array<double, 1>& gradient);
      private:
            GradientTemplate();
            GradientTemplate(GradientTemplate<T> const &);
            void operator=(GradientTemplate<T>);
            T * object_;
            GradientMember member_;
      };

      template <class T>
      GradientTemplate<T>::GradientTemplate(T& object, GradientMember member) :
            object_(&object),
            member_(member)
      {}
      
      template <class T>
      void GradientTemplate<T>::compute(blitz::Array<double, 1>& coordinates, blitz::Array<double, 1>& gradient)
      {
            (object_->*member_)(coordinates, gradient);
      }
         
      namespace random {
            unsigned int seedWithTime();
            unsigned int seed(unsigned int seed);
            /// Gaussian distribution of mean zero and variance 1.
            extern ranlib::Normal<double> normal;
            /// Uniform distribution in interval [0;1]
            extern ranlib::Uniform<double> uniform;
      };
      
      /** @defgroup gslblitz  Wrappers between blitz and gsl.
            Make gsl_vector and gsl_matrix visible to blitz and vice versa.
            @section fastWrapping Fast Wrapping for const Array.
            When wrapping, the memory in blitz::Array is shared with the gsl_vector (or gsl_matrix). Though marked as constant it is still possible change the data through pointer gsl_vector.data. Doing this would affect the blitz::Array declared as constant. It is an error. A safer solution is to copy the blitz beforehand and use the standard wrapToBlitz function. 
            @{ */
      /// Wrap gsl_vector into blitz::Array
      blitz::Array<double, 1> wrapToBlitz(gsl_vector * const vgsl);
      /// Wrap gsl_matrix into blitz::Array
      blitz::Array<double, 2> wrapToBlitz(gsl_matrix * const mgsl);
      /// Wrap blitz::Array into gsl_vector.
      gsl_vector wrapToGsl(blitz::Array<double, 1> & vblitz);
      /// Wrap blitz::Array into gsl_matrix.
      gsl_matrix wrapToGsl(blitz::Array<double, 2> & mblitz);
      
      /// Wrap blitz::Array into gsl_vector.
      /// @see @ref fastWrapping
      gsl_vector const fastWrapToGsl(blitz::Array<double, 1> const & vblitz);
      /// Wrap blitz::Array into gsl_matrix.
      /// @see @ref fastWrapping
      gsl_matrix const fastWrapToGsl(blitz::Array<double, 2> const & mblitz);
      /// @}
}
#endif
