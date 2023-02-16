// serves as an interface between LAMMPS potentials maintained by SANDIA

#ifndef LAMMPS
#define LAMMPS

#include "../../Potential.h"
#include "../../Parameters.h"

class lammps : public Potential {

    public:
        lammps(Parameters *p);
        ~lammps(void);
        void initialize() {};
        void cleanMemory(void);
        void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
        std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                             const VectorXi atmnrs,
                                             const Matrix3d box) override {
          double energy{std::numeric_limits<double>::infinity()};
          long nAtoms{pos.rows()};
          AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
          this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy,
                      box.data());
          return std::make_pair(energy, forces);
        };

    private:
        long numberOfAtoms;
        double oldBox[9];
        void *LAMMPSObj;
        void makeNewLAMMPS(long N, const double *R,  const int *atomicNrs, const double *box);
        Parameters *parameters;
        bool realunits;
};
#endif
