// -*- c++ -*-

// TODO These functions will be removed in the near future and replaced wtih
// LGPL-compatible functions

#define COLVARS_OK 0
#define COLVARS_ERROR 1

namespace NR_Jacobi {

  /// Numerical recipes diagonalization
  int jacobi(double a[4][4], double d[4], double v[4][4], int *nrot);

  /// Eigenvector sort
  int eigsrt(double d[4], double v[4][4]);

  /// Transpose the matrix
  int transpose(double v[4][4]);

}

