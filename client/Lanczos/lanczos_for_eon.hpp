/** @file
      Interface Lanczos for EON.     
      @author Jean Claude C. Berthet
      @date 2007
      University of Iceland
      */
#ifndef LANCZOS_FOR_EON_HPP
#define LANCZOS_FOR_EON_HPP
#include "../LowestEigenmodeInterface.h"
#include "lanczos.hpp"
/** Lanczos class compatible with EON.
      */
class Lanczos : private gradient_scanning::Lanczos, public LowestEigenmodeInterface {
public:
    /** Constructor.
    @param[in]      *parameters                Pointer to the Parameter object containing the runtime parameters.*/    
    Lanczos(Matter *const, Parameters *parameters);
    
    ~Lanczos() {};///< Destructor.
    void virtual startNewSearchAndCompute(Matter const *matter, Matrix<double, Eigen::Dynamic, 3> matrix); 
    void virtual moveAndCompute(Matter const *matter);  
    double virtual getEigenvalue();
    /// Return eigenvector.
    virtual Matrix<double, Eigen::Dynamic, 3> getEigenvector();
        /** Set initial direction manually.*/
    virtual void setEigenvector(Matrix<double, Eigen::Dynamic, 3> const eigenvector);
private:
    Lanczos();
    blitz::Array<double, 1> eigenvector_;///< Back up eigenvector for last Lanczos.
    double eigenvalue_;
    Matter matter_;
};
#endif
