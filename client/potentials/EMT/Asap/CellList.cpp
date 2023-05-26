#include "CellList.h"
#include "Atoms.h"
#include "Exception.h"
#include "GhostAtoms.h"
#include "Image.h"
#include "SuperCell.h"
#include <map>
#include <math.h>
using std::map;
#include <assert.h>
using std::cerr;
using std::endl;

CellList::CellList(Atoms *atoms, double r) : rMax(r), rMax2(r * r) {
  nAtoms = atoms->GetNumberOfRealAtoms();
  nAllAtoms = nAtoms;
  int i;
  // Do we have any ghosts?
  GhostAtoms *ghostAtoms;
  //  Atoms* tmpAtoms = atoms;

  //  try {
  //	ghostAtoms = dynamic_cast<GhostAtoms *>(tmpAtoms);

  //  }catch(...)
  //  {
  ghostAtoms = NULL;
  //  }
  if (ghostAtoms != 0) {
    nAllAtoms += ghostAtoms->GetNumberOfGhosts();
  }

  superCell = atoms->GetSuperCell();
  positions = atoms->GetPositionsPtr();
  periodic = superCell->periodic;
  superCellHeights = superCell->heights;
  for (i = 0; i < 3; i++)
    if (periodic[i] && superCellHeights[i] < 2 * rMax)
      throw Exception("The height of the cell (")
          << superCellHeights[i] << ") must be larger than " << 2 * rMax;

  superCell->TransformToScaledSpace(&positions[0], nAllAtoms);
  for (i = 0; i < 3; i++) {
    if (periodic[i]) {
      size[i] = 1.0;
      minimum[i] = 0.0;
#if 0
            double *p = &(positions[0][i]);
            for (int a = 0; a < nAllAtoms; a++)
            {
                *p -= floor(*p);
                p += 3;
            }
#endif
    } else {
      double *p = &(positions[0][i]);
      double min = *p;
      double max = *p;
      p += 3;
      for (int a = 1; a < nAllAtoms; a++) {
        double x = *p;
        if (x > max)
          max = x;
        else if (x < min)
          min = x;
        p += 3;
      }
      minimum[i] = min;
      size[i] = max - min;
    }
    nCells[i] = int(size[i] * superCellHeights[i] / rMax);
    if (nCells[i] == 0) {
      nCells[i] = 1;
      minimum[i] -= 0.5 * (rMax / superCellHeights[i] - size[i]);
      size[i] = rMax / superCellHeights[i];
    }
  }
  nTotalCells[0] = 1;
  for (i = 0; i < 3; i++)
    nTotalCells[i + 1] = nTotalCells[i] * nCells[i];
  cells.resize(nTotalCells[3]);
  const double *x = &(positions[0][0]);
  extern int verbose;
  if (verbose > 1) {
    cerr << minimum << " <---> " << minimum + size << endl
         << nAtoms << " atoms and " << nAllAtoms - nAtoms << " ghosts in "
         << nCells[0] << " * " << nCells[1] << " * " << nCells[2] << " cells"
         << periodic[0] << periodic[1] << periodic[2] << endl;
  }
  cellIndices.resize(nAllAtoms);
  for (int a = 0; a < nAllAtoms; a++) {
    int index = 0;
    for (i = 0; i < 3; i++) {
      int k = int(((*x++) - minimum[i]) / size[i] * nCells[i]);
      assert(k >= 0);
      if (k == nCells[i]) {
        k--;
      }
      assert(k < nCells[i]);
      index += nTotalCells[i] * k;
    }
    cells[index].push_back(a);
    cellIndices[a] = index;
  }
  superCell->TransformToRealSpace(&positions[0], nAllAtoms);
}

// Don't try to understand this method:
void CellList::MakeNeighborList(int neighborList[], int indices[], int &nMax,
                                vector<Image> &images)

{
  const Vec *translations = superCell->translations;
  images.resize(0);
  vector<map<int, int>> imageNumbers(27);
  int n = 0;
  indices[0] = 0;
  for (int a1 = 0; a1 < nAtoms; a1++) {
    Vec pos1 = positions[a1];
    int index1 = cellIndices[a1];
    int i0 = index1 % nCells[0];
    int i1 = (index1 / nTotalCells[1]) % nCells[1];
    int i2 = index1 / nTotalCells[2];
    for (int j2 = i2 - 1; j2 <= i2 + 1; j2++) {
      int k2 = j2;
      if (periodic[2]) {
        if (j2 < 0)
          k2 = nCells[2] - 1;
        else if (j2 >= nCells[2])
          k2 = 0;
      } else if (j2 < 0 || j2 >= nCells[2])
        continue;
      for (int j1 = i1 - 1; j1 <= i1 + 1; j1++) {
        int k1 = j1;
        if (periodic[1]) {
          if (j1 < 0)
            k1 = nCells[1] - 1;
          else if (j1 >= nCells[1])
            k1 = 0;
        } else if (j1 < 0 || j1 >= nCells[1])
          continue;
        for (int j0 = i0 - 1; j0 <= i0 + 1; j0++) {
          int k0 = j0;
          if (periodic[0]) {
            if (j0 < 0)
              k0 = nCells[0] - 1;
            else if (j0 >= nCells[0])
              k0 = 0;
          } else if (j0 < 0 || j0 >= nCells[0])
            continue;
          int index2 = k0 + nTotalCells[1] * k1 + nTotalCells[2] * k2;
          const vector<int> &cell2 = cells[index2];
          if (n >= nMax)
            throw Exception("Too many neighbor pairs: ") << n << " pairs";
          if (k0 == j0 && k1 == j1 && k2 == j2) {
            for (unsigned int n2 = 0; n2 < cell2.size(); n2++) {
              int a2 = cell2[n2];
              if (a2 <= a1)
                continue;
              if (Length2(positions[a2] - pos1) < rMax2)
                neighborList[n++] = a2;
              assert(n < nMax);
            }
          } else {
            int diff[3] = {(j0 - k0), (j1 - k1), (j2 - k2)};
            int dn = 1;
            int nTranslation = 1 + 3 + 3 * 3;
            for (int i = 0; i < 3; i++, dn *= 3)
              if (periodic[i])
                nTranslation += diff[i] / nCells[i] * dn;
            map<int, int> &iNumbers = imageNumbers[nTranslation];
            Vec tpos1 = pos1 - translations[nTranslation];
            for (unsigned int n2 = 0; n2 < cell2.size(); n2++) {
              int a2 = cell2[n2];
              if (nTranslation > 27 / 2 && a2 < nAtoms)
                continue;
              if (Length2(positions[a2] - tpos1) < rMax2) {
                int a3;
                if (iNumbers.find(a2) == iNumbers.end()) {
                  a3 = nAllAtoms + images.size();
                  images.push_back(Image(a2, nTranslation));
                  iNumbers[a2] = a3;
                } else
                  a3 = iNumbers[a2];
                neighborList[n++] = a3;
              }
            }
          }
        }
      }
    }
    indices[a1 + 1] = n;
  }
  extern int verbose;
  if (verbose > 1)
    cerr << images.size() << " images" << endl
         << n << " neighbor pairs (" << nMax << ')' << endl;
  nMax = n;
}

bool CellList::AnyNeighbors(Vec point) const {
  int neighbors[100];
  return (GetNeighbors(point, neighbors, 100) > 0);
}

int CellList::GetNeighbors(Vec point, int *neighbors, int space) const {
  int ii[3];
  int i;
  Vec center(0.0, 0.0, 0.0);
  for (i = 0; i < 3; i++) {
    double scaled = point * superCell->inverse[i];
    if (periodic[i]) {
      scaled -= floor(scaled);
      ii[i] = int(scaled * nCells[i]);
      if (ii[i] == nCells[i])
        ii[i] = nCells[i] - 1;
    } else
      ii[i] = int(floor((scaled - minimum[i]) / size[i] * nCells[i]));
    // center += scaled * superCell->vectors[i];
    Vaxpy(scaled, superCell->vectors[i], center);
  }
  int nNeighbors = 0;
  for (int j2 = ii[2] - 1; j2 <= ii[2] + 1; j2++) {
    int k2 = j2;
    if (periodic[2]) {
      if (j2 < 0)
        k2 = nCells[2] - 1;
      else if (j2 >= nCells[2])
        k2 = 0;
    } else if (j2 < 0 || j2 >= nCells[2])
      continue;
    for (int j1 = ii[1] - 1; j1 <= ii[1] + 1; j1++) {
      int k1 = j1;
      if (periodic[1]) {
        if (j1 < 0)
          k1 = nCells[1] - 1;
        else if (j1 >= nCells[1])
          k1 = 0;
      } else if (j1 < 0 || j1 >= nCells[1])
        continue;
      for (int j0 = ii[0] - 1; j0 <= ii[0] + 1; j0++) {
        int k0 = j0;
        if (periodic[0]) {
          if (j0 < 0)
            k0 = nCells[0] - 1;
          else if (j0 >= nCells[0])
            k0 = 0;
        } else if (j0 < 0 || j0 >= nCells[0])
          continue;
        int index = k0 + nTotalCells[1] * k1 + nTotalCells[2] * k2;
        const vector<int> &cell = cells[index];
        int diff[3] = {(j0 - k0), (j1 - k1), (j2 - k2)};
        Vec tCenter = center;
        for (i = 0; i < 3; i++)
          if (diff[i] != 0)
            tCenter -= superCell->vectors[i] * (diff[i] / nCells[i]);
        for (unsigned int n = 0; n < cell.size(); n++) {
          int a = cell[n];
          if (Length2(positions[a] - tCenter) < rMax2) {
            assert(space > 0);
            neighbors[nNeighbors++] = a;
            space--;
          }
        }
      }
    }
  }
  return nNeighbors;
}
