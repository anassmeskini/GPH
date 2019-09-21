#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include <algorithm>
#include <cassert>
#include <cstring>
#include <vector>

struct SparseMatrix
{
   SparseMatrix() = default;

   SparseMatrix(SparseMatrix&&) noexcept;

   SparseMatrix& operator=(SparseMatrix&&) noexcept;

#ifdef UNIT_TEST
   SparseMatrix& operator=(SparseMatrix&) = default;
#endif

   int ncols;
   int nrows;
   std::vector<double> coefficients;
   std::vector<int> indices;
   std::vector<int> rowStart;
};

#endif
