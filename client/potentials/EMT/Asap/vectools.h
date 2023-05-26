/*
 * SimQA: Copyright (C) Jakob Schiotz 1996, 1997
 *
 * IBM POWER2 optimized vector library:
 * Header files for all the vectorized functions
 */

/* XXX Added extern"C" to everything? */

#ifndef VECTOOLS_H
#define VECTOOLS_H

#ifdef _MSC_VER
// XXX Need to disable under Microsoft console and _WIN32: what is the analog
// for _MICROSOFT
#pragma warning(disable : 4068) // turns off spurious Windows pragma warnings
#endif                          // _MSC_VER

/*   res[] = a[] * b + c */
inline void vec_mul_scal_add_scal(double res[], const double a[],
                                  const double b, const double c, int n) {
  for (int i = 0; i < n; i++)
    res[i] = a[i] * b + c;
}

/* r[] = a[] * b[offset[]] + c[offset[]] */
inline void vec_mul_indir_add_indir(double res[], const double a[],
                                    const double b[], const double c[],
                                    const int offset[], int n) {
  for (int i = 0; i < n; i++)
    res[i] = a[i] * b[offset[i]] + c[offset[i]];
}

/*   res[] = a[] * b[] */
inline void vec_mul(double res[], const double a[], const double b[], int n) {
  for (int i = 0; i < n; i++)
    res[i] = a[i] * b[i];
}

/*   res[] = a[] * b */
inline void vec_mul_scal(double res[], const double a[], double b, int n) {
  for (int i = 0; i < n; i++)
    res[i] = a[i] * b;
}

/* r[] = a[] * b[offset[]] */
inline void vec_mul_indir(double res[], const double a[], const double b[],
                          const int offset[], int n) {
  for (int i = 0; i < n; i++)
    res[i] = a[i] * b[offset[i]];
}

/*   res[] = res[] * a[] + b * c[] */
inline void vec_eq_self_mul_addscal_mul(double res[], const double a[],
                                        double b, const double c[], int n) {
  for (int i = 0; i < n; i++)
    res[i] = res[i] * a[i] + b * c[i];
}

/* r[] = r[] * a[] * b[offset[]] */
inline void vec_self_mul_mul_indir(double res[], const double a[],
                                   const double b[], const int offset[],
                                   int n) {
  for (int i = 0; i < n; i++)
    res[i] = res[i] * a[i] * b[offset[i]];
}

/*   res[] = (a[] - a[] * a[]) * b */
inline void vec_diffweight(double res[], const double a[], double b, int n) {
  for (int i = 0; i < n; i++)
    res[i] = (a[i] - a[i] * a[i]) * b;
}

/*   res[] = (a[] * b + c[]) * d[] */
inline void vec_mul_scal_add_mul(double res[], const double a[], double b,
                                 const double c[], const double d[], int n) {
  for (int i = 0; i < n; i++)
    res[i] = (a[i] * b + c[i]) * d[i];
}

/* r[] = (a[] * b[offset[]] + c[offset[]]) * d[] */
inline void vec_mul_indir_add_indir_mul(double res[], const double a[],
                                        const double b[], const double c[],
                                        const double d[], const int offset[],
                                        int n) {
  for (int i = 0; i < n; i++)
    res[i] = (a[i] * b[offset[i]] + c[offset[i]]) * d[i];
}

/* r[] = a[] * b[offset[]] + c[] * d[offset[]] */
inline void vec_dbl_mul_indir_add(double res[], const double a[],
                                  const double b[], const double c[],
                                  const double d[], const int offset[], int n) {
  for (int i = 0; i < n; i++)
    res[i] = a[i] * b[offset[i]] + c[i] * d[offset[i]];
}

/*  res[] = (d[offset[]] * a[] + e[offset[]]) * b[] + f[offset[]] * c[]  */
inline void vec_dEds(double res[], const double a[], const double b[],
                     const double c[], const double d[], const double e[],
                     const double f[], const int offset[], int n) {
  for (int i = 0; i < n; i++)
    res[i] = (d[offset[i]] * a[i] + e[offset[i]]) * b[i] + f[offset[i]] * c[i];
}

/*  r[] = a[] * (b[] * c[] * d + e[] * f + g[] * h[] * x + y[] * z)  */
/* Note that sometimes b==g and e==y !!!!! */
inline void vec_df(double *res, const double *a, const double *b,
                   const double *c, double d, const double *e, double f,
                   const double *g, const double *h, double x, const double *y,
                   double z, int n) {
  for (int i = 0; i < n; i++)
    res[i] = a[i] * (b[i] * c[i] * d + e[i] * f + g[i] * h[i] * x + y[i] * z);
}

/* r[] = a[] + b */
inline void vec_add_scal(double res[], const double a[], double b, int n) {
  for (int i = 0; i < n; i++)
    res[i] = a[i] + b;
}

#endif // VECTOOLS_H
