#include "Vec.h"
#include <iostream>

ostream &operator<<(ostream &out, const Vec &v) {
  out << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
  return out;
}

istream &operator>>(istream &in, Vec &v) {
  in >> v[0] >> v[1] >> v[2];
  return in;
}
