#include <math.h> 
#include <iostream>

#include "../../PotentialsInterface.h"
    
class EAM : public PotentialsInterface
{    
    public:
        EAM(void);
        // To satify interface
        void initialize();
        void cleanMemory();
        void force(long N, const double *R, const long *atomicNrs, double *F,
                   double *U, const double *box);
    private:
        // Variables
	    long *celllist_old;
	    long *celllist_new;
	    long *neigh_list;
	    bool initialized;
	    double *rc;
        void calc_force(long N, double *R, const long *atomicNrs, double *F,
                        double *U, const double *box);
        void new_celllist(long N, const double *box, long *num_axis,
                          long *cell_length, long *celllist_new, 
                          long num_cells, double *Rnew);
        void cell_to_neighbor(long N, long num_of_cells, long *num_axis,
                              long *cell_length, long *celllist_new,
                              long *neigh_list);
        //returns 0 if unchanged, >0 if changed - represents number of atoms that changed lists
        int update_cell_list(long N, long num_cells,long *num_axis,
                             long *cell_length, long *celllist_old, double *Rnew);
        //calculates local density of single atom
        double density (long N, long atom, double *R, const long *atomicNrs,
                        const double *box);
        double embedding_function(double *func_coeff, double rho);
        double embedding_force(double *func_coeff, double rho);
};
