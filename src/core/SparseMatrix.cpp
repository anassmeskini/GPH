#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(SparseMatrix&& other) noexcept
{
   ncols = other.ncols;
   nrows = other.nrows;

   coefficients = std::move(other.coefficients);
   indices = std::move(other.indices);
   rowStart = std::move(other.rowStart);
}

SparseMatrix&
SparseMatrix::operator=(SparseMatrix&& other) noexcept
{
   ncols = other.ncols;
   nrows = other.nrows;

   coefficients = std::move(other.coefficients);
   indices = std::move(other.indices);
   rowStart = std::move(other.rowStart);

   return *this;
}
