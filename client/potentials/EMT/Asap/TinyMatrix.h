// TinyMatrix.h  ---  Implements a rudimentary matrix class

#ifndef _TINYMATRIX_H
#define _TINYMATRIX_H

template <class Type> class TinyMatrix {
public:
  TinyMatrix() { data = 0; }
  TinyMatrix(int rows, int columns) { Allocate(rows, columns); }
  void Allocate(int rows, int columns) {
    r = rows;
    c = columns;
    data = new Type[r * c];
  }
  ~TinyMatrix() { delete[] data; }
  Type *operator[](int row) { return data + c * row; }
  const Type *operator[](int row) const { return data + c * row; }

private:
  int r, c;
  Type *data;
};

typedef TinyMatrix<double> TinyDoubleMatrix;

#endif // _TINYMATRIX_H
