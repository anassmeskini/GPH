#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include <algorithm>
#include <cassert>
#include <cstring>
#include <vector>

template<typename REAL>
struct VectorView;

template<typename REAL>
struct SparseMatrix
{
   SparseMatrix() {}

   SparseMatrix(size_t, size_t);

   SparseMatrix(SparseMatrix<REAL>&&);

   SparseMatrix<REAL>& operator=(SparseMatrix<REAL>&&);

   void addRow(VectorView<REAL>);

   size_t ncols;
   size_t nrows;
   std::vector<REAL> coefficients;
   std::vector<size_t> indices;
   std::vector<size_t> rowStart;
};

template<typename REAL>
SparseMatrix<REAL>::SparseMatrix(size_t nnz, size_t nCols)
  : ncols(nCols)
  , nrows(0)
{
   coefficients.resize(nnz);
   indices.resize(nnz);
   rowStart.push_back(0);
}

template<typename REAL>
void
SparseMatrix<REAL>::addRow(VectorView<REAL> row)
{
   assert(all_of(row.indices, row.indices + row.size, [&](size_t id) {
      return id < ncols;
   }));

   size_t nnz = coefficients.size();
   rowStart.push_back(nnz);

   auto coefBegin = coefficients.data() + nnz;
   auto indBegin = indices.data() + nnz;

   coefficients.resize(nnz + row.size);
   indices.resize(nnz + row.size);

   memcpy(coefBegin, row.coefs, sizeof(REAL) * row.size);
   memcpy(indBegin, row.indices, sizeof(size_t) * row.size);

   ++nrows;
}

template<typename REAL>
SparseMatrix<REAL>::SparseMatrix(SparseMatrix<REAL>&& other)
{
   ncols = other.ncols;
   nrows = other.nrows;

   coefficients = std::move(other.coefficients);
   indices = std::move(other.indices);
   rowStart = std::move(other.rowStart);
}

template<typename REAL>
SparseMatrix<REAL>&
SparseMatrix<REAL>::operator=(SparseMatrix<REAL>&& other)
{
   ncols = other.ncols;
   nrows = other.nrows;

   coefficients = std::move(other.coefficients);
   indices = std::move(other.indices);
   rowStart = std::move(other.rowStart);

   return *this;
}

#endif
