#ifndef __VEC_H__
#define __VEC_H__

#include <iostream>
using std::istream;
using std::ostream;
// class istream;
// class ostream;

/// A 3-vector useful for postions etc.

/// The only data is the three positions (and there are no virtual
/// functions), so the memory layout of an array of Vecs will be x0,
/// y0, z0, x1, y1, z1, x2, ...
///
/// Almost all operations are inline for speed.

class Vec {
public:
  /// Dummy constructor needed by STL containers.
  Vec(){};
  /// Construct a 3-vector from three doubles.
  Vec(double x0, double x1, double x2);
  /// Dot product
  double operator*(const Vec &v) const;
  /// Multiplication with scalar
  Vec operator*(const double &s) const;
  /// Add two Vecs
  Vec operator+(const Vec &v) const;
  /// Subtract two Vecs
  Vec operator-(const Vec &v) const;
  /// Unary minus
  Vec operator-() const;
  /// Add a Vec to this one.
  Vec &operator+=(const Vec &v);
  /// Subtract a vec from this one.
  Vec &operator-=(const Vec &v);
  /// Multiply this vec with a scalar
  Vec &operator*=(double s);
  /// Divide this Vec with a scalar.
  Vec &operator/=(double s);
  /// const indexing
  double operator[](int n) const;
  /// Non-const indexing
  double &operator[](int n);
  /// Cross product of two Vecs.
  friend Vec Cross(const Vec &v1, const Vec &v2);
  /// The length of a Vec
  friend double Length2(const Vec &v);
  /// Increment y with a times x.
  friend void Vaxpy(double a, const Vec &x, Vec &y);
  /// Print a Vec
  friend ostream &operator<<(ostream &out, const Vec &v);
  /// Read a Vec
  friend istream &operator>>(istream &in, Vec &v);

public:
  double x[3]; ///< The actual data.
};

inline Vec::Vec(double x0, double x1, double x2) {
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
}

inline double Vec::operator*(const Vec &v) const {
  return x[0] * v.x[0] + x[1] * v.x[1] + x[2] * v.x[2];
}

inline Vec Vec::operator*(const double &s) const {
  return Vec(s * x[0], s * x[1], s * x[2]);
}

inline Vec Vec::operator+(const Vec &v) const {
  return Vec(x[0] + v.x[0], x[1] + v.x[1], x[2] + v.x[2]);
}

inline Vec Vec::operator-(const Vec &v) const {
  return Vec(x[0] - v.x[0], x[1] - v.x[1], x[2] - v.x[2]);
}

inline Vec Vec::operator-() const { return Vec(-x[0], -x[1], -x[2]); }

inline Vec &Vec::operator+=(const Vec &v) {
  x[0] += v.x[0];
  x[1] += v.x[1];
  x[2] += v.x[2];
  return *this;
}

inline Vec &Vec::operator-=(const Vec &v) {
  x[0] -= v.x[0];
  x[1] -= v.x[1];
  x[2] -= v.x[2];
  return *this;
}

inline Vec &Vec::operator*=(double s) {
  x[0] *= s;
  x[1] *= s;
  x[2] *= s;
  return *this;
}

inline Vec &Vec::operator/=(double s) {
  x[0] /= s;
  x[1] /= s;
  x[2] /= s;
  return *this;
}

inline double Vec::operator[](int n) const { return x[n]; }

inline double &Vec::operator[](int n) { return x[n]; }

inline Vec Cross(const Vec &v1, const Vec &v2) {
  return Vec(v1.x[1] * v2.x[2] - v1.x[2] * v2.x[1],
             v1.x[2] * v2.x[0] - v1.x[0] * v2.x[2],
             v1.x[0] * v2.x[1] - v1.x[1] * v2.x[0]);
}

inline void Vaxpy(double a, const Vec &x, Vec &y) {
  y.x[0] += a * x.x[0];
  y.x[1] += a * x.x[1];
  y.x[2] += a * x.x[2];
}

inline double Length2(const Vec &v) {
  return v.x[0] * v.x[0] + v.x[1] * v.x[1] + v.x[2] * v.x[2];
}

#endif // __VEC_H__
