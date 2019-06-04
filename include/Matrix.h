// Sam Svindland
//
// Verified to run on MacOS 10.14.2
//                    Apple LLVM version 10.0.0 (clang-1000.11.45.5)
//                    cmake 3.13.2
//                    VTK 8.1.2
#include <iostream>
#include <vtkPolyDataReader.h>

#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

#endif