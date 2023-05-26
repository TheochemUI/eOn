// Emacs: This is -*- C++ -*-

/// Represents atoms seen through the periodic boundaries.

/// The Image structure is used by NeighborList to represent atoms
/// seen through the periodic boundary conditions.
struct Image {
  /// Dummy constructor needed by some STL containers.
  Image(){};
  /// Construct an image of atom n translated by translation number nTrans.
  Image(int num, int nTrans) : number(num), nTranslation(nTrans) {}

  /// The number of the real atom.
  int number;

  /// The number of the translation to apply to the real atom.
  int nTranslation;
};
