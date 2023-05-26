/// \file EMT.h  --  Implements the EMT potential for multiple elements.

#ifndef _EMT_H
#define _EMT_H

#include "EMTParameterProvider.h"
#include "Potential.h"
#include "TinyMatrix.h"
#include <vector>
using std::vector;

class AtomsBase;
class Atoms;
class GhostAtoms;
class NeighborList;
struct emt_parameters;

/// The Effective Medium Theory (EMT) potential.

/// The Effective Medium Theory (EMT) as originally implemented in Per
/// Stoltze's ARTwork program, and documented in
///   - K. W. Jacobsen, P. Stoltze and J. K. Nørskov,
///     Surf. Sci. vol. 366, p. 394-402 (1996).
class EMT : public AsapPotential {
public:
  /// Create an EMT potential optionally using a parameter provider.
  EMT(EMTParameterProvider *prov = 0);

  /// Delete the EMT potential.
  virtual ~EMT();
  void SetSubtractE0(bool subtractE0) { this->subtractE0 = subtractE0; }

  virtual void SetAtoms(Atoms *atoms);
  virtual const Vec *GetCartesianForces();
  virtual const double *GetPotentialEnergies();
  virtual double GetPotentialEnergy();
  virtual const symTensor *GetStresses(const Vec *momenta = 0);
  virtual void GetStress(double stress[6], const Vec *momenta = 0);
  /// Check that neighbor lists are up to date, update them if not.
  virtual void CheckNeighborLists();

  virtual double GetCutoffRadius() const { return rNbCut; }
  virtual double GetLatticeConstant() const;
  virtual int GetNumberOfAtoms() const { return nAtoms; }

  /// \brief Specify the element in the continuum region of a
  /// QuasiContinuum simulation.
  virtual void SetContinuumElement(int z);
  virtual void UpdateSuperCell(const SuperCell *newSuperCell);

  // Quasicontinuum stuff
  virtual double CalculateLatticeEnergy(const Vec a[3]);
  virtual void CalculateDerivatives(const Vec a[3], double dEdaDota[6]);
  virtual double GetData() const { return dataSigma1; }

  /// Return a pointer to the EMT "density" sigma1.
  const double *GetSigma1(int n) { return &sigma1[n][0]; }
  /// Return a pointer to the EMT "density" sigma2.
  const double *GetSigma2(int n) { return &sigma2[n][0]; }

  /// Print the EMT parameters.
  void PrintParameters();

  /// Return a pointer to the neighbor list.
  NeighborList *GetNeighborList() const { return nblist; }

protected:
  /// Initialization of the EMT parameters.
  virtual void InitParameters();
  /// (Re)allocate storage for forces, energies and intermediate results.
  virtual void Allocate();
  /// (Re)allocate storage for stresses
  virtual void AllocateStress();
  /// Calculate type numbers from the atomic numbers.
  virtual void CalculateIDs();
  /// Calculate sigma1 and perhaps sigma2.
  virtual void CalculateSigmas(int calculatesigma2 = 1);
  /// Calculate energies once sigma1 is known

  /// If Epot is NULL, energies are not calculated, only derivatives
  /// needed for the force.
  virtual void CalculateEnergiesAfterSigmas(double *Epot = 0);

  /// Calculate forces in a multicomponent system.
  virtual void CalculateForcesAfterEnergies(Vec *forces,
                                            double stress[][6] = 0);

  /// Calculate forces in a system with only one element.
  virtual void CalculateForcesAfterEnergiesSingle(Vec *forces,
                                                  double stress[][6] = 0);

  /// Internal function used for the QuasiContinuum stuff.
  virtual void CalculateStuff(const Vec a[3], double &sigma1, double &sigma2,
                              double *diffs2, double *wght, double *dsigma1,
                              double *dsigma2);

private: // Methods.
  /// sigma_batch does the hard work in CalculateSigmas().

  ///
  /// Not virtual, do not reimplement without reimplementing
  /// CalculateSigmas().
  void sigma_batch(int *self, int *other, Vec rnb[], double *sq_dist, int zs,
                   int zo, int n, int calculatesigma2);

  /// force_batch does the hard work in CalculateForcesAfterEnergy().

  ///
  /// Not virtual, do not reimplement without reimplementing
  /// CalculateForcesAfterEnergy().
  void force_batch(int *self, int *other, Vec rnb[], double *sq_dist,
                   double dEdss[], double dEdso[], int zs, int zo, int n,
                   Vec *force, double (*stress)[6]);

protected:                        // Data
  Atoms *atoms;                   ///< The atoms we are working on
  GhostAtoms *ghostatoms;         ///< Non-NULL if atoms are GhostAtoms.
  int nAtoms;                     ///< The number of (real) atoms.
  int nSize;                      ///< Number of atoms including ghost atoms.
  NeighborList *nblist;           ///< The neighborlist object.
  EMTParameterProvider *provider; ///< The source of the EMT parameters
  bool ownProvider;               ///< May we delete the provider?

  bool subtractE0; ///< Whether we subtract E0 from atomic energies (defines the
                   ///< zero of potential energy; if we don't subtract, zero
                   ///< corresponds to infinite seperation. If we do, it is the
                   ///< bulk fcc crystal at the equilibrium lattice constant.)

  std::vector<const emt_parameters *> parameters; ///< The EMT parameters
  const emt_parameters *continuumelement;
  ///< The element in the continuum region
  const TinyDoubleMatrix *chi; ///< The Chi matrix of EMT.
  int nelements; ///< The number of different elements in the simulation
  /// Cutoff parameters (from EMTParameterProvider).
  double rFermi, rNbCut;
  /// Cutoff slope
  double cutoffslope;
  //  int maxleaflen;         // Lenght of the longest leaf.
  int sigma2isvalid; ///< For consistency checks.

  /// Temporary data for the atoms

  /// Each atom has two sigmas for each possible neighboring element.
  /// So these are arrays (nelements long).
  vector<vector<double>> sigma1;
  vector<vector<double>> sigma2;

  //@{
  /// Each atom has a single Ec, Eas, radius etc.
  vector<double> Ec;
  vector<double> Eas;
  vector<double> potentialenergy;
  vector<double> radius;
  vector<double> dEds;
  vector<Vec> force;
  /// The stresses.  Always remember to use six per atoms!

  /// Should really be declared as vector<double[6]>, but that breaks
  /// the vector library when attempting to resize!
  vector<double> stress;
  //@}

  int nAtomsRes, nSizeRes;
  ///< If there are ghostatoms, some extra space is reserved for the arrays

  /// The atomic numbers are translated into IDs, integers in [0, nelements-1]
  vector<int> id;

  /* The buffers for batch processing */
  // static const int BUFLEN;  // The batch--buffer size.  Avoid powers of 2
  // static const int NMAXELEMENTS;

  double unnormalizedstress[6]; ///< The total stress before division by volume
  double totalvolume;           ///< Total volume during stresscalculation

  int nHalfNeighbors; ///< For the quasicontinuum stuff
  int (*coef)[6];     ///< For the quasicontinuum stuff
  double dataSigma1;  ///< For the quasicontinuum stuff

  /// A structure of counters to check if recalculations are necessary.
  struct {
    int ids;
    int nblist;
    int sigma1;
    int sigma2;
    int beforeforces;
    int energies;
    int forces;
    int stresses;
    int fullstresses;
  } counters;
};

#endif // ! _EMT_H
