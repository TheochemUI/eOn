#include <stdio.h>

#include "EAM.h"

int main() {
  EAM pot;
  pot.initialize();

  long N = 2;
  const double R[] = {0.0, 0.0, 0.0, 5.0, 5.0, 5.0};
  const long atomicNrs[] = {13, 13};
  double *F = new double[3 * N];
  double U = 0;
  const double box[] = {10.0, 10.0, 10.0};
  pot.force(N, R, atomicNrs, F, &U, box);
  printf("Energy: %10.4g\n", U);
  printf("Forces: ");
  for (int i = 0; i < 3 * N; i++) {
    if (i % 3 == 0)
      printf("\n");
    printf("%10.4g ", F[i]);
  }
  printf("\n");
  pot.cleanMemory();
  return 0;
}
