
#ifndef BONDBOOST_H
#define BONDBOOST_H

#include "HelperFunctions.h"
#include "Matter.h"
#include "Parameters.h"

#include "Eigen.h"

/** Functionality relying on the conjugate gradients algorithm. The object is
 * capable of minimizing an Matter object or modified forces being passed in.*/
class BondBoost {

public:
  /** Constructor to be used when a structure is minimized.
  @param[in]   *matter        Pointer to the Matter object to be relaxed.
  @param[in]   *parameters    Pointer to the Parameter object containing the
  runtime parameters.*/
  BondBoost(Matter *matt, Parameters *params);
  ~BondBoost(); ///< Destructor.

  void initialize();
  double boost();

private:
  Matrix<double, Eigen::Dynamic, 1> Rmdsteps();
  long BondSelect();
  double Booststeps();
  long nAtoms;    ///< Number of free coordinates.
  Matter *matter; ///< Pointer to atom object \b outside the scope of the class.
  Parameters *parameters; ///< Pointer to a structure outside the scope of the
                          ///< class containing runtime parameters.
  long *BAList;
  long *RAList;
  long *TABAList;
  long *BBAList;
  double *Epsr_Q;
  Matrix<double, Eigen::Dynamic, 1>
      TABLList; // EquilibriumTaggedAtomInvolvedBondLengthList;
  Matrix<double, Eigen::Dynamic, 1> EBBLList; // EquilibriumBoostBondLengthList
  Matrix<double, Eigen::Dynamic, 1> CBBLList; // CurrentBoostBondLengthList
  long nBAs;
  long nRAs;
  long nTABs;
  long nReg;
  long nBBs;
  std::shared_ptr<spdlog::logger> log;
};

class Hyperdynamics {
public:
  //            static const string NONE;
  //            static const string BOND_BOOST;
  static const char NONE[];
  static const char BOND_BOOST[];
};

#endif
