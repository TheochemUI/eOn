/*
 * SimQA: Copyright (C) Jakob Schiotz 1996, 1997
 *
 * IBM POWER2 optimized vector library:
 * Vectorized mathematica functions from the IBM mass library.
 */

#include <math.h>

inline void vsqrt(double res[], const double a[], const int *n) {
  int nn = *n;
  for (int i = 0; i < nn; i++)
    res[i] = sqrt(a[i]);
}

inline void vexp(double res[], const double a[], const int *n) {
  int nn = *n;
  for (int i = 0; i < nn; i++)
    res[i] = exp(a[i]);
}

inline void vrec(double res[], const double a[], const int *n) {
  int nn = *n;
  for (int i = 0; i < nn; i++)
    res[i] = 1.0 / a[i];
}

inline void vlog(double res[], const double a[], const int *n) {
  int nn = *n;
  for (int i = 0; i < nn; i++)
    res[i] = log(a[i]);
}
