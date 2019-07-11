#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(int nnz, int nCols)
  : ncols(nCols)
  , nrows(0)
{
   coefficients.resize(nnz);
   indices.resize(nnz);
   rowStart.push_back(0);
}

SparseMatrix::SparseMatrix(SparseMatrix&& other)
{
   ncols = other.ncols;
   nrows = other.nrows;

   coefficients = std::move(other.coefficients);
   indices = std::move(other.indices);
   rowStart = std::move(other.rowStart);
}

SparseMatrix&
SparseMatrix::operator=(SparseMatrix&& other)
{
   ncols = other.ncols;
   nrows = other.nrows;

   coefficients = std::move(other.coefficients);
   indices = std::move(other.indices);
   rowStart = std::move(other.rowStart);

   return *this;
}
