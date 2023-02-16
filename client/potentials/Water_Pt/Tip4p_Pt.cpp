#include "Tip4p_Pt.hpp"

std::pair<double, AtomMatrix> Tip4p_Pt::get_ef(const AtomMatrix pos,
                                               const VectorXi atmnrs,
                                               const Matrix3d m_box) {
  double energy{0};
  long N{pos.rows()};
  AtomMatrix forces{Eigen::MatrixXd::Zero(N, 3)};
  const double *R = pos.data();
  const double *box = m_box.data();
  const int *atomicNrs = atmnrs.data();
  double *F = forces.data();
  double *U = &energy;

  double diagbox[3];
  diagbox[0] = box[0];
  diagbox[1] = box[4];
  diagbox[2] = box[8];
  int i = 0;
  while (atomicNrs[i] == 1)
    i += 2;
  computeHH_O_Pt_(i / 2, N - i * 3 / 2, R, F, *U, diagbox, 0);
}
