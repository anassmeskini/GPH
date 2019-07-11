#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include <algorithm>
#include <cassert>
#include <cstring>
#include <vector>

struct SparseMatrix
{
   SparseMatrix() {}

   SparseMatrix(int, int);

   SparseMatrix(SparseMatrix&&);

   SparseMatrix& operator=(SparseMatrix&&);

   int ncols;
   int nrows;
   std::vector<double> coefficients;
   std::vector<int> indices;
   std::vector<int> rowStart;
};

#endif
