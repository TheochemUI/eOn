#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "EAM.h"
#include "Parameters.h"

void EAM::cleanMemory() {
  if (initialized) {
    delete[] celllist_old;
    delete[] celllist_new;
    delete[] neigh_list;
  }
}

// Calculate here long num_cells, long *num_axis, long *cell_length, //become
// global variables -long *celllist_old, long *celllist_new, long *neigh_list,
// long fcalled)
void EAM::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, double *variance, const double *fullbox) {
  variance = nullptr;
  double box[3];
  box[0] = fullbox[0];
  box[1] = fullbox[4];
  box[2] = fullbox[8];

  /* -- Code starting here only needs to be done once. -- */
  // num_axis contains the number of cell lengths on each axis
  long num_axis[3];
  // cell length contains dimensions of each cell
  long cell_length[3];

  for (long i = 0; i < 3; i++) {
    num_axis[i] =
        (long)(box[i] / rc[i]) +
        1; // sets num-axis according to box size and optimal dimensions of cell
  }

  long num_cells;
  for (long i = 0; i < 3; i++)
    cell_length[i] =
        (long)(box[i] /
               (num_axis[i] -
                1)); // from num_axis, find the best fit for dimensions of cell

  num_cells = num_axis[0] * num_axis[1] * num_axis[2];

  /* -- Code ending here only needs to be done once. -- */

  if (!initialized) // first time force is called
  {
    celllist_old = new long[num_cells * (N + 1)];
    celllist_new = new long[num_cells * (N + 1)];
    neigh_list = new long[N * (N + 1)];
  }

  *U = 0;

  // Are the passed in forces already initialized?
  for (long kk = 0; kk < 3 * N; kk++) {
    F[kk] = 0;
  }

  long xmin = LONG_MAX;
  long ymin = LONG_MAX;
  long zmin = LONG_MAX;

  // shift of the R values so that all coordinates are positive (makes creating
  // cell table easier)
  double *Rtemp = new double[3 * N];
  double *Rnew = new double[3 * N];
  double *Rold = new double[3 * N];

  // finds the minimum value on each axis. This is necessary because within the
  // method, the R values will be shifted so all are positive, making the
  // creation of a cell list really easy.
  for (long i = 0; i < 3 * N; i += 3) {
    Rtemp[i] = R[i];
    Rtemp[i + 1] = R[i + 1];
    Rtemp[i + 2] = R[i + 2];
    if (R[i] < xmin)
      xmin = R[i];
    if (R[i + 1] < ymin)
      ymin = R[i + 1];
    if (R[i + 2] < zmin)
      zmin = R[i + 2];
  }
  // if no coordinates are less than 0, no shift is necessary
  if (xmin > 0)
    xmin = 0;
  if (ymin > 0)
    ymin = 0;
  if (zmin > 0)
    zmin = 0;
  // shifting of R
  for (long i = 0; i < 3 * N; i += 3) {
    Rtemp[i] += abs(xmin);
    Rtemp[i + 1] += abs(ymin);
    Rtemp[i + 2] += abs(zmin);
  }

  // enforce periodic boundary conditions
  for (long i = 0; i < 3 * N; i++) {
    while (Rtemp[i] > box[i % 3]) {
      Rtemp[i] -= box[i % 3];
    }
  }

  for (long i = 0; i < 3 * N; i++) {
    Rnew[i] = Rtemp[i];
  }

  if (!initialized) {
    new_celllist(N, box, num_axis, cell_length, celllist_new, num_cells, Rnew);
    cell_to_neighbor(N, num_cells, num_axis, cell_length, celllist_new,
                     neigh_list);
  } else {
    // only updates cell and neighbor lists if an atom moved cells
    if (update_cell_list(N, num_cells, num_axis, cell_length, celllist_old,
                         Rnew) > 0) {
      new_celllist(N, box, num_axis, cell_length, celllist_new, num_cells,
                   Rnew);
      cell_to_neighbor(N, num_cells, num_axis, cell_length, celllist_new,
                       neigh_list);
    }
  }
  // does force calculations for all atoms
  calc_force(N, Rnew, atomicNrs, F, U, box);

  for (long i = 0; i < num_cells * (N + 1); i++) {
    celllist_old[i] = celllist_new[i];
  }

  // double diff = R[0] - R[3];
  // if(diff > box[0]/2)
  //{
  //     diff = box[0] - diff;
  // }
  delete[] Rtemp;
  delete[] Rnew;
  delete[] Rold;
  initialized = true;
}

EAM::element_parameters EAM::get_element_parameters(int atomic_number) {
  bool found = false;
  int i;
  for (i = 0; i < NPARAMS; i++) {
    if (el_params[i].Z == atomic_number) {
      found = true;
      break;
    }
  }
  if (found == false) {
    /* This sucks. We need to have a way to alert user that the
     * parameters are missing */
    throw 14324;
  }
  return el_params[i];
}

void EAM::calc_force(long N, double *R, const int *atomicNrs, double *F,
                     double *U, const double *box) {
  double *drho_dr = new double[3 * N];
  for (long i = 0; i < 3 * N; i++) {
    F[i] = 0;
  }

  *U = 0;

  for (long i = 0; i < N; i++) {
    for (long j = 0; j < 3 * N; j++) {
      drho_dr[j] = 0;
    }
    element_parameters params = get_element_parameters(atomicNrs[i]);

    double dens = 0; // sum of density for atom 1
    double phi_r = 0;

    // for (long j=0;j<neigh_list[i*(N+1)+N];j++) //loops through each atom in
    // atom i's neighbor list
    //{
    //     long neigh= neigh_list[i*(N+1)+j];
    for (long j = 0; j < N;
         j++) // loops through each atom in atom i's neighbor list
    {
      if (i == j) {
        continue;
      }
      phi_r = 0;
      long neigh = j;

      // finds distance on each axis between atoms i and j
      double disx = R[3 * i] - R[3 * neigh];
      double disy = R[3 * i + 1] - R[3 * neigh + 1];
      double disz = R[3 * i + 2] - R[3 * neigh + 2];
      double dirs[] = {disx, disy, disz};

      // finds smallest distance between the two atoms (periodic boundary
      // conditions)
      for (long u = 0; u < 3; u++) {
        if (dirs[u] > box[u] / 2) {
          dirs[u] = -box[u] + dirs[u];
        } else if (dirs[u] < -box[u] / 2) {
          dirs[u] = box[u] + dirs[u];
        }
      }

      // distance between atoms i and j
      double r = 0;

      for (long k = 0; k < 3; k++) {
        r += dirs[k] * dirs[k];
      }
      r = sqrt(r);
      if (r > params.r_cut) {
        continue;
      }
      dens += pow(r, 6) *
              (exp(-params.beta1 * r) + 512 * exp(-2 * params.beta2 * r));
      double mag_force_den =
          6 * pow(r, 5) *
              (exp(-params.beta1 * r) + 512 * exp(-2 * params.beta2 * r)) +
          pow(r, 6) * (-params.beta1 * exp(-params.beta1 * r) -
                       1024 * params.beta2 * exp(-2 * params.beta2 * r));
      for (long k = 0; k < 3; k++) {
        drho_dr[3 * i + k] += dirs[k] / r * -mag_force_den;
        drho_dr[3 * neigh + k] -= dirs[k] / r * -mag_force_den;
      }
      if (neigh > i) {
        // Morse potential portion of energy
        phi_r = params.Dm * pow(1 - exp(-params.alphaM * (r - params.Rm)), 2) -
                params.Dm;

        // magnitude of force from Morse potential
        double mag_force =
            2 * params.alphaM * params.Dm *
            (exp(params.alphaM * r) - exp(params.alphaM * params.Rm)) *
            (exp(params.alphaM * params.Rm - 2 * params.alphaM * r));

        // calculates forces on atom i and j because of each other
        for (long k = 0; k < 3; k++) {
          F[3 * i + k] -= dirs[k] / r * mag_force;
          F[3 * neigh + k] += dirs[k] / r * mag_force;
        }
        *U += phi_r;
      }
    }

    double f_of_rho = embedding_function(params.func_coeff, dens);
    *U += f_of_rho;
    double dF_drho = embedding_force(params.func_coeff, dens);
    for (int k = 0; k < 3 * N; k++) {
      F[k] -= dF_drho * drho_dr[k];
    }
  }
  delete[] drho_dr;
}

// Creates a new cell list for the first time force is called
void EAM::new_celllist(long N, const double *box, long *num_axis,
                       long *cell_length, long *celllist_new, long num_cells,
                       double *Rnew) {
  // the cell list is a 1d array - each N+1 indexes correspond to 1 cell. The
  // last index of each cell holds the number of atoms in that cell
  long *cell_list = new long[num_cells * (N + 1)];

  for (long i = 0; i < num_cells; i++)
    for (long j = 0; j < N + 1; j++)
      if (j % N == 0 && j > 0)
        cell_list[i * (N + 1) + j] = 0;
      else
        cell_list[i * (N + 1) + j] = -1;

  // loop to put each atom into appropriate cell
  for (long i = 0; i < N; i++) {
    long *coor = new long[3];
    for (long j = 0; j < 3; j++)
      // this calculates, for each axis, how many lengths away from 0,0,0 the
      // atom is
      coor[j] = (long)((Rnew[3 * i + j]) / cell_length[j]);

    // this calculates which cell number the atom is in
    int cell =
        coor[0] * num_axis[1] * num_axis[2] + coor[1] * num_axis[2] + coor[2];

    // puts atom in appropriate spot of appropriate cell
    cell_list[cell * (N + 1) + cell_list[cell * (N + 1) + N]] = i;
    // increases the value that gives number of atoms in that particular cells
    cell_list[cell * (N + 1) + N]++;

    delete[] coor;
  }

  for (long i = 0; i < num_cells; i++)
    for (long j = 0; j < N + 1; j++)
      celllist_new[i * (N + 1) + j] = cell_list[i * (N + 1) + j];
  delete[] cell_list;
}

// creates a neighbor list for each cell from given cell list
void EAM::cell_to_neighbor(long N, long num_of_cells, long *num_axis,
                           long *cell_length, long *celllist_new,
                           long *neigh_list) {
  num_of_cells = num_axis[0] * num_axis[1] * num_axis[2];

  for (long i2 = 0; i2 < N; i2++)
    for (long k = 0; k < N + 1; k++)
      neigh_list[i2 * (N + 1) + k] = -1;

  for (long j = 0; j < num_axis[0]; j++) {
    for (long j1 = 0; j1 < num_axis[1]; j1++) {
      for (long j2 = 0; j2 < num_axis[2]; j2++) {
        // last index contains number of atoms in array.
        long *neighbors = new long[N + 1];
        for (long i = 0; i < N + 1; i++) {
          neighbors[i] = -1;
        }

        long cur_index = j * num_axis[1] * num_axis[2] + j1 * num_axis[2] +
                         j2; // converts 3d location of cell to 1d location in
                             // use in cell_list

        // making a copy of cell list so that no cell is counted twice as a
        // neighbor of a dif cell in small atom systems
        long *cell_list_copy = new long[num_of_cells * (N + 1) + N + 1];
        for (long d1 = 0; d1 < num_of_cells; d1++) {
          for (long d2 = 0; d2 < N + 1; d2++) {
            cell_list_copy[d1 * (N + 1) + d2] = celllist_new[d1 * (N + 1) + d2];
          }
        }

        if (cell_list_copy[cur_index * (N + 1) + N] == 0) {
          delete[] cell_list_copy;
          delete[] neighbors;
          continue; // if cell size is 0, continues to next cell
        }

        for (long d1 = -1; d1 < 2; d1++) {
          for (long d2 = -1; d2 < 2; d2++) {
            for (long d3 = -1; d3 < 2; d3++) {
              long *pos = new long[3];
              pos[0] = j - d1;
              pos[1] = j1 - d2;
              pos[2] =
                  j2 -
                  d3; // 3d location of neighbor cell
                      // check if pos is within bounds, if not, pos wraps around
                      // to meet neighbors (periodic boundary conditions)
              for (long y = 0; y < 3; y++) {
                if (pos[y] < 0) {
                  pos[y] = num_axis[y] - 1;
                } else if (pos[y] >= num_axis[y]) {
                  pos[y] = 0;
                }
              }
              // 1d location of neighbor cell
              long neigh_index = pos[0] * num_axis[1] * num_axis[2] +
                                 pos[1] * num_axis[2] + pos[2];

              // number of atoms in neighbor cell
              long num = cell_list_copy[neigh_index * (N + 1) + N];
              for (long y = 0; y < num; y++) {
                // makes sure that atom has not already been added to neigbors
                long check = 0;
                for (long k = 0; k < N; k++) {
                  if (neighbors[check] ==
                      cell_list_copy[neigh_index * (N + 1) + y]) {
                    check = -1;
                  }
                }
                // if not, add atom to neighbors
                if (check != -1) {
                  neighbors[N]++;
                  neighbors[neighbors[N]] =
                      cell_list_copy[neigh_index * (N + 1) + y];
                }
              }

              delete[] pos;
            }
          }
        }
        // adds all atoms in neighbors to the array of each atom in current cell
        for (long i = 0; i < celllist_new[cur_index * (N + 1) + N]; i++) {
          // atom for which neighbors are being added to neigh_list
          long cur = celllist_new[cur_index * (N + 1) + i];
          neigh_list[cur * (N + 1) + N] = 0;
          for (long temp = 0; temp <= neighbors[N]; temp++) {
            if (cur != neighbors[temp]) // checks to make sure not same atom
            {

              neigh_list[cur * (N + 1) + neigh_list[cur * (N + 1) + N]] =
                  neighbors[temp];
              neigh_list[cur * (N + 1) + N]++;
            }
          }
        }
        delete[] neighbors;
        delete[] cell_list_copy;
      }
    }
  }
}

// Returns 0 if unchanged, >0 if changed - represents number of atoms that
// changed lists
int EAM::update_cell_list(long N, long num_cells, long *num_axis,
                          long *cell_length, long *celllist_old, double *Rnew) {
  int changed = 0;

  long *table =
      new long[N]; // holds the cell number for each atom in the old cell list

  for (long i = 0; i < num_cells; i++)
    for (long j = 0; j < celllist_old[i * (N + 1) + N]; j++)
      table[celllist_old[i * (N + 1) + j]] = i;

  for (long i = 0; i < N; i++) {
    long *coor = new long[3];
    for (long j = 0; j < 3; j++)
      // this calculates, for each axis, how many lengths away from 0,0,0 the
      // atom is
      coor[j] = (long)((Rnew[3 * i + j]) / cell_length[j]);
    // this calculates which cell number the atom is in
    long cell =
        coor[0] * num_axis[1] * num_axis[2] + coor[1] * num_axis[2] + coor[2];
    // if cell location is different, changed increments
    if (cell != table[i]) {
      changed++;
    }
    delete[] coor;
  }
  delete[] table;
  return changed;
}

double EAM::embedding_function(const double *func_coeff, double rho) {
  double result = func_coeff[8];
  // Evaluate the polynomial using Horner's method.
  for (int i = 7; i >= 0; i--) {
    result = result * rho + func_coeff[i];
  }
  return result;
}

double EAM::embedding_force(const double *func_coeff, double rho) {
  double result = func_coeff[8] * 8;
  // Evaluate the polynomial using Horner's method.
  for (int i = 7; i >= 1; i--) {
    result = result * rho + i * func_coeff[i];
  }
  return -result;
}
