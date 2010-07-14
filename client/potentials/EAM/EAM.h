#include <math.h> 
#include <iostream>

//#include "../../system_unit.h" // unit converters
#include "../../PotentialsInterface.h"
    
class EAM : public PotentialsInterface
{    
    public:
        EAM(void);
        // Variables
	    long *celllist_old;
	    long *celllist_new;
	    long *neigh_list;
	    bool initialized;
	    double *rc;
        // To satify interface
        void initialize();
        void cleanMemory();
        void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);
        void calc_force(long N, double *R, const long *atomicNrs, double *F, double *U, const double *box);
        void new_celllist(long N, const double *box, long *num_axis, long *cell_length, long *celllist_new, long num_cells, double *Rnew);
        void cell_to_neighbor(long N, long num_of_cells, long *num_axis, long *cell_length, long *celllist_new, long *neigh_list);
        int update_cell_list(long N, long num_cells,long *num_axis, long *cell_length, long *celllist_old, double *Rnew); //returns 0 if unchanged, >0 if changed - represents number of atoms that changed lists
        double density (long N, long atom, double *R, const long *atomicNrs, const double *box); //calculates local density of single atom
};
	

