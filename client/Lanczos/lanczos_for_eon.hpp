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
    
    /// Execute Lanczos. It uses gradient to start lanczos. Use at the beggining of a saddle point search.
    virtual void startNewSearchAndCompute(Matter const *matter, double *);
    /// Execute Lanczos. It uses last eigenvector to start Lanczos. Use in subsequent step of a saddle point search
    virtual void moveAndCompute(Matter const *matter);
    /** Lowest eigenvalue and corresponding eigenvector.
        @param[in] result Pointer to array to store the eigenvector. The length of the vector is @f$ 3\times \textrm{Number of Movable Atoms} @f$.
        @return Eigenvalue
        */
    virtual double returnLowestEigenmode(double *result);
private:
    Lanczos();
    blitz::Array<double, 1> eigenvector_;///< Back up eigenvector for last Lanczos.
    double eigenvalue_;
    Matter matter_;
};
#endif
