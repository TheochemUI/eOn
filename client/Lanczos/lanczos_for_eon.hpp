/** @file
      Interface Lanczos for EON.     
      @author Jean Claude C. Berthet
      @date 2007
      University of Iceland
      */
#ifndef LANCZOS_FOR_EON_HPP
#define LANCZOS_FOR_EON_HPP
#include "../LowestEigenmodeInterface.h"
#include "oldlanczos.hpp"

/** Lanczos class compatible with EON.
      */
class OldLanczos : private gradient_scanning::Lanczos, public LowestEigenmodeInterface {
public:
    /** Constructor.
    @param[in]      *parameters                Pointer to the Parameter object containing the runtime parameters.*/    
    OldLanczos(Matter *const, Parameters *parameters);
    ~OldLanczos() {};///< Destructor.

    void compute(Matter const *matter, AtomMatrix direction);

    double getEigenvalue();
    /// Return eigenvector.
    Matrix<double, Eigen::Dynamic, 3> getEigenvector();
        /** Set initial direction manually.*/
    void setEigenvector(Matrix<double, Eigen::Dynamic, 3> const eigenvector);
private:
    OldLanczos();
    blitz::Array<double, 1> eigenvector_;///< Back up eigenvector for last Lanczos.
    double eigenvalue_;
    Matter matter_;
};
#endif
