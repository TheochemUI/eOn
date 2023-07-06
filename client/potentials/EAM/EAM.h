#include <iostream>
#include <math.h>

#ifndef EAM_STANDALONE
#include "../../Potential.h"
#endif

class EAM
#ifndef EAM_STANDALONE
    : public Potential
#endif
{
public:
  EAM(std::shared_ptr<Parameters> params)
      : Potential(PotType::EAM_AL, params) {
    celllist_new = 0;
    neigh_list = 0;
    initialized = false;

    rc = new double[3];
    // 6 is arbitrary number. rc represents the optimal size
    // for each cell in cell list.
    rc[0] = rc[1] = rc[2] = 6.0;
  };
  // To satify interface
  void cleanMemory();
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *fullbox) override;

private:
  struct element_parameters {
    const int Z;                // Atomic Number
    const double Dm;            // Morse potential well depth
    const double alphaM;        // Curvative at Morse minimum
    const double Rm;            // Position of Morse minimum
    const double beta1;         // Density parameter 1
    const double beta2;         // Density parameter 2
    const double r_cut;         // Cutoff distance
    const double func_coeff[9]; // 8th order poly for embedding function
  };
  static const element_parameters el_params[];
  // Variables
  long *celllist_old;
  long *celllist_new;
  long *neigh_list;
  bool initialized;
  double *rc;
  void calc_force(long N, double *R, const int *atomicNrs, double *F, double *U,
                  const double *box);
  void new_celllist(long N, const double *box, long *num_axis,
                    long *cell_length, long *celllist_new, long num_cells,
                    double *Rnew);
  void cell_to_neighbor(long N, long num_of_cells, long *num_axis,
                        long *cell_length, long *celllist_new,
                        long *neigh_list);
  // returns 0 if unchanged, >0 if changed - represents number of atoms that
  // changed lists
  int update_cell_list(long N, long num_cells, long *num_axis,
                       long *cell_length, long *celllist_old, double *Rnew);
  // calculates local density of single atom
  double density(long N, long atom, double *R, const long *atomicNrs,
                 const double *box);
  double embedding_function(const double *func_coeff, double rho);
  double embedding_force(const double *func_coeff, double rho);
  element_parameters get_element_parameters(int atomic_number);
};
