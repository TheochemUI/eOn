// Emacs: this is -*- C++ -*-

#include "Image.h"
#include <vector>
using std::vector;

class SuperCell;
class Atoms;
class Vec;
// class Image;

/// "Half" neighbor lists for atoms.

/// This implements a <em>"half neighbor list"</em>, i.e. if atom a is
/// on atom b's neighbor list, atom b will NOT be on atom a's neighbor
/// list.
///
/// GetNeighbors(...) returns "half" of the neighbors inside a sphere
/// of radius rCut. CheckAndUpdateNeighborLists() must be called every
/// time the positions change.  If an atom has moved more than
/// "drift", CheckAndUpdateNeighborLists() will call MakeList(), and a
/// new list is generated.

class NeighborList {
public:
  /// Generate a neighbor list for atoms a with cutoff rCut.

  /// The neighbor list will contain all neighbors within the distance
  /// rCut.  The neighborlist can be reused until an atom has moved
  /// more than rCut*driftfactor.
  NeighborList(Atoms *a, double rCut, double driftfactor = 0.05,
               int defaultNeighborEstimate = 50);
  ~NeighborList();

  /// Check if the neighbor list can still be reused, update if not.
  bool CheckAndUpdateNeighborList();

  /// Get all neighbors of atom n.  The most important method :-)

  /// Input values: n is the number of the atom.  r (optional) is a
  /// cutoff, must be less than rCut in the constructor (not
  /// checked!).
  ///
  /// In-out values: size contains the maximum space in the arrays.
  /// It is decremented by the number of neighbors placed in the
  /// arrays.  It is an error to call GetNeighbors with too small a
  /// value of size.
  ///
  /// Out values: neighbors[] contains the numbers of the atoms,
  /// diffs[] contains the \em relative positions of the atoms,
  /// diffs2[] contains the norms of the diffs vectors.
  ///
  /// Return value: The number of neighbors.
  int GetNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2, int &size,
                   double r = -1.0) const;

  /// Inform the neighbor list that the supercell has changed.
  void UpdateSuperCell(const SuperCell *superCell);

  /// \em Low-level interface to the neighbor list.

  /// Returns a const pointer to the internal list of neighbors of
  /// atom a, the length of the list is placed in n.  The list will
  /// contain all neighbors within the cutoff radius, and some
  /// neighbors outside the cutoff.  Some neighbors will be image
  /// atoms instead of real atoms (if they are neighbors thanks to the
  /// periodic boundaries), use GetImage to get the number of the real
  /// atom.
  const int *GetList(int a, int &n);

  /// Get the real atom corresponding to an image.

  /// The input argument i is the number of the image counted from the
  /// first image.  It is thus the number returned by GetList
  /// <em>minus the number of atoms</em>.  The return value is a const
  /// reference to an Image.
  const Image &GetImage(int i) { return images[i]; }

  /// Return the guaranteed maximal length of a single atom's NB list.

  /// Call this before using GetNeighbors() to make sure the arrays
  /// are big enough.  The value may change when the neighbor list is
  /// updated.
  int MaxNeighborListLength() const { return maxlistlen; }

  /// Get the number of atoms in the corresponding list of atoms.
  int GetNumberOfAtoms() const { return nAtoms; } // Used from swig.
private:
  /// Generate a new neighbor list.
  void MakeList();
  /// Update the positions of the images when the real atoms may have moved.
  void UpdateImagePositions();
  /// (Re)allocate the memory for the neighbor list.
  void Allocate();

private:
  Atoms *atoms;   ///< A pointer to the atoms.
  int nAtoms;     ///< The number of atoms.
  int nSize;      ///< Number of atoms including ghost atoms.
  double rCut2;   ///< The square of the cutoff radius.
  double drift2;  ///< The square of the maximally allowed drift of an atom.
  double rMax;    ///< Cutoff when \em generating neighbor list.
  Vec *positions; ///< A pointer into the positions array of the atoms.
  vector<Vec> referencePositions; ///< The positions at the last update.
  vector<Image> images;           ///< The images of neighbors through PBC.
  int *neighborList;              ///< The actual neighbor list.
  int nMax;                       ///< The size of the neighborList array.
  int *indices;   ///< Offsets into neighborList giving each atom's neighbors.
  int maxlistlen; ///< The length of the longest neighbor list of an atom.
};
