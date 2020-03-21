#include "MIP.h"
#include "Common.h"
#include "Numerics.h"
#include "io/Message.h"

#include <algorithm>
#include <numeric>

MIP::MIP(const Rows& rows, const Cols& cols, std::vector<double>&& coefsT,
         std::vector<int>&& idxT, std::vector<int>&& rstartT,
         std::vector<double>&& _rhs, std::vector<double>&& _lhs,
         std::vector<double>&& _lbs, std::vector<double>&& _ubs,
         std::vector<double>&& _obj, const dynamic_bitset<>& integer,
         std::vector<int>& rowSize, std::vector<std::string>&& _varNames)
{
   int ncols = cols.size();
   int nrows = rows.size();

   assert(ncols == static_cast<int>(rstartT.size() - 1));
   assert(idxT.size() == coefsT.size());

   // move constraint sides
   lhs = std::move(_lhs);
   rhs = std::move(_rhs);

   // move column bounds
   lb = std::move(_lbs);
   ub = std::move(_ubs);

   // move the objective
   objective = std::move(_obj);

   // fill the column major matrix
   constMatrixT.coefficients = std::move(coefsT);
   constMatrixT.indices = std::move(idxT);
   constMatrixT.rowStart = std::move(rstartT);
   constMatrixT.ncols = nrows;
   constMatrixT.nrows = ncols;
   stats.ncols = ncols;
   stats.nrows = nrows;

   for (int col = 0; col < ncols; ++col)
   {
      if (constMatrixT.rowStart[col + 1] == constMatrixT.rowStart[col])
         Message::warn("column {} has zero support", col);
   }

   // get the row-major matrix by transposing
   constMatrix = transpose(constMatrixT, rowSize);

   // fill coefficients statistics
   stats.nnzmat = constMatrixT.coefficients.size();
   for (auto cost : objective)
   {
      if (cost != 0.0)
         ++stats.nnzobj;
   }

   // move constraint names
   consNames.resize(nrows);
   for (auto& elem : rows)
   {
      int row = elem.second.second;
      consNames[row] = std::move(elem.first);
   }

   // move variable names
   assert(ncols == static_cast<int>(_varNames.size()));
   varNames.resize(ncols);
   for (int i = 0; i < ncols; ++i)
      varNames[i] = std::move(_varNames[i]);

   // compute the locks and max/min row support
   downLocks.resize(ncols);
   upLocks.resize(ncols);
   for (int row = 0; row < nrows; ++row)
   {
      auto [rowcoefs, indices, rowsize] = getRow(row);

      bool lhsfinite = !Num::isMinusInf(lhs[row]);
      bool rhsfinite = !Num::isInf(rhs[row]);

      if (lhs[row] == rhs[row])
         ++stats.nequality;

      for (int id = 0; id < rowsize; ++id)
      {
         int col = indices[id];
         double coef = rowcoefs[id];
         assert(rowcoefs[id] != 0.0);

         if (lhsfinite)
         {
            if (coef > 0.0)
               ++downLocks[col];
            else
               ++upLocks[col];
         }
         if (rhsfinite)
         {
            if (coef > 0.0)
               ++upLocks[col];
            else
               ++downLocks[col];
         }
      }

      stats.avgRowSupport += rowsize;
      if (rowsize > stats.maxRowSupport)
         stats.maxRowSupport = rowsize;
      if (rowsize < stats.minRowSupport)
         stats.minRowSupport = rowsize;
   }

   for (int col = 0; col < ncols; ++col)
   {
      auto [colcoefs, colindices, size] = getCol(col);

      if (std::max(upLocks[col], downLocks[col]) > stats.maxLocks)
         stats.maxLocks = std::max(upLocks[col], downLocks[col]);
      if (std::min(upLocks[col], downLocks[col]) < stats.minLocks)
         stats.minLocks = std::min(upLocks[col], downLocks[col]);

      stats.avgColSupport += size;
      if (size > stats.maxColSupport)
         stats.maxColSupport = size;
      if (size < stats.minColSupport)
         stats.minColSupport = size;

      if (integer[col])
      {
         if (lb[col] == 0.0 && ub[col] == 1.0)
            ++stats.nbin;
         else
            ++stats.nint;
      }
      else
         ++stats.ncont;
   }

   stats.avgRowSupport /= nrows;
   stats.avgColSupport /= ncols;

   // reorder variables | binary | int | continuous
   // permutation maps new id -> old id
   std::vector<int> permutation(ncols);
   std::generate(std::begin(permutation), std::end(permutation), [](){
      static int i = 0;
      return i++;
   });

   std::sort(permutation.begin(), permutation.end(),
             [&integer, this](int left, int right) -> bool {
                bool leftbin =
                    integer[left] && lb[left] == 0.0 && ub[left] == 1.0;
                bool rightbin =
                    integer[right] && lb[right] == 0.0 && ub[right] == 1.0;

                bool leftint = integer[left];
                bool rightint = integer[right];

                return (leftbin && (!rightbin || !rightint)) ||
                       (leftint && !rightint);
             });

   // permute all information related to columns
   // permutation maps old id -> new id
   std::vector<int> mapping(ncols);
   dynamic_bitset<> int_buf(ncols);
   std::vector<double> lb_buf(ncols);
   std::vector<double> ub_buf(ncols);
   std::vector<int> dl_buf(ncols);
   std::vector<int> ul_buf(ncols);
   std::vector<double> obj_buf(ncols);
   std::vector<std::string> name_buf(ncols);
   for (int i = 0; i < ncols; ++i)
   {
      mapping[permutation[i]] = i;

      lb_buf[i] = lb[permutation[i]];
      ub_buf[i] = ub[permutation[i]];
      obj_buf[i] = objective[permutation[i]];
      name_buf[i] = std::move(varNames[permutation[i]]);
      dl_buf[i] = downLocks[permutation[i]];
      ul_buf[i] = upLocks[permutation[i]];
   }

   lb = std::move(lb_buf);
   ub = std::move(ub_buf);
   objective = std::move(obj_buf);
   varNames = std::move(name_buf);
   downLocks = std::move(dl_buf);
   upLocks = std::move(ul_buf);

   // change the old indices by the new ones in the row major matrix
   for (size_t i = 0; i < constMatrix.indices.size(); ++i)
   {
      constMatrix.indices[i] = mapping[constMatrix.indices[i]];
   }

   // permute the rows of the columns major matrix
   SparseMatrix new_transp;
   new_transp.ncols = nrows;
   new_transp.nrows = ncols;
   new_transp.coefficients.resize(stats.nnzmat);
   new_transp.indices.resize(stats.nnzmat);
   new_transp.rowStart.push_back(0);
   int nnz = 0;
   for (int i = 0; i < ncols; ++i)
   {
      int old_col = permutation[i];
      int size = constMatrixT.rowStart[old_col + 1] -
                 constMatrixT.rowStart[old_col];

      std::memcpy(new_transp.coefficients.data() + nnz,
                  constMatrixT.coefficients.data() +
                      constMatrixT.rowStart[old_col],
                  size * sizeof(double));

      std::memcpy(new_transp.indices.data() + nnz,
                  constMatrixT.indices.data() +
                      constMatrixT.rowStart[old_col],
                  size * sizeof(int));

      nnz += size;
      new_transp.rowStart.push_back(nnz);
   }

   assert(nnz == stats.nnzmat);
   constMatrixT = std::move(new_transp);

#ifndef NDEBUG
   printStats(stats);
#endif
}

MIP&
MIP::operator=(MIP&& other) noexcept
{
   objective = std::move(other.objective);

   lhs = std::move(other.lhs);
   rhs = std::move(other.rhs);

   ub = std::move(other.ub);
   lb = std::move(other.lb);

   varNames = std::move(other.varNames);
   consNames = std::move(other.consNames);

   constMatrix = std::move(other.constMatrix);
   constMatrixT = std::move(other.constMatrixT);

   upLocks = std::move(other.upLocks);
   downLocks = std::move(other.downLocks);

   stats = other.stats;

   return *this;
}

SparseMatrix
MIP::transpose(const SparseMatrix& matrix, const std::vector<int>& rowSize)
{
   int nnz = matrix.coefficients.size();
   int ncols = matrix.nrows;
   int nrows = matrix.ncols;

   assert(matrix.rowStart.size() == static_cast<size_t>(matrix.nrows + 1));
   assert(std::accumulate(rowSize.begin(), rowSize.end(), 0) == nnz);
   assert(std::all_of(matrix.coefficients.begin(),
                      matrix.coefficients.end(),
                      [](double coef) { return coef != 0.0; }));

   SparseMatrix transposed;
   transposed.nrows = nrows;
   transposed.ncols = ncols;
   transposed.coefficients.resize(nnz);
   transposed.indices.resize(nnz);

   auto& rowStart = transposed.rowStart;
   rowStart.reserve(nrows + 1);
   rowStart.push_back(0);

   std::vector<int> offset(nrows, 0);

   for (int i = 0; i < nrows; ++i)
      rowStart.push_back(rowSize[i] + rowStart[i]);

   for (int col = 0; col < ncols; ++col)
   {
      for (int rowid = matrix.rowStart[col];
           rowid < matrix.rowStart[col + 1]; ++rowid)
      {
         int row = matrix.indices[rowid];
         double coef = matrix.coefficients[rowid];
         int rowstart = rowStart[row];

         assert(coef != 0.0);
         transposed.coefficients[rowstart + offset[row]] = coef;
         transposed.indices[rowstart + offset[row]] = col;
         ++offset[row];
         assert(offset[row] <= rowSize[row]);
      }
   }

   assert(std::all_of(transposed.coefficients.begin(),
                      transposed.coefficients.end(),
                      [](double coef) { return coef != 0.0; }));

   return transposed;
}

VectorView
MIP::getRow(int row) const noexcept
{
   const double* coefBegin =
       constMatrix.coefficients.data() + constMatrix.rowStart[row];
   const int* indBegin =
       constMatrix.indices.data() + constMatrix.rowStart[row];
   const int size =
       constMatrix.rowStart[row + 1] - constMatrix.rowStart[row];

   return {coefBegin, indBegin, size};
}

VectorView
MIP::getCol(int col) const noexcept
{
   const double* coefBegin =
       constMatrixT.coefficients.data() + constMatrixT.rowStart[col];
   const int* indBegin =
       constMatrixT.indices.data() + constMatrixT.rowStart[col];
   const int size =
       constMatrixT.rowStart[col + 1] - constMatrixT.rowStart[col];

   return {coefBegin, indBegin, size};
}

int
MIP::getRowSize(int row) const noexcept
{
   return constMatrix.rowStart[row + 1] - constMatrix.rowStart[row];
}

int
MIP::getColSize(int col) const noexcept
{
   return constMatrixT.rowStart[col + 1] - constMatrixT.rowStart[col];
}

void
printStats(Statistics st)
{
   int ncols = st.ncols;
   int nrows = st.nrows;

   Message::print("\n==Statistics==");

   float percbin = 100.0 * (static_cast<double>(st.nbin) / ncols);
   float percint = 100.0 * (static_cast<double>(st.nint) / ncols);
   float perccont = 100.0 * (static_cast<double>(st.ncont) / ncols);
   Message::print("Columns: {}", st.ncols);
   Message::print(
       "\tbin: {} ({:0.1f}%)\n\tint: {} ({:0.1f}%)\n\tcont: {} ({:0.1f}%)",
       st.nbin, percbin, st.nint, percint, st.ncont, perccont);

   Message::print("\tavg support: {}\n\tmin support: {}\n\tmax support: "
                  "{}\n\tmin locks: {}\n\tmax locks: {}",
                  st.avgColSupport, st.minColSupport, st.maxColSupport,
                  st.minLocks, st.maxLocks);

   float perceq = 100.0 * (static_cast<double>(st.nequality) / nrows);
   Message::print("Rows: {}", st.nrows);
   Message::print("\tneq: {} ({:0.1f}%)\n\tavg support: {}\n\tmin support"
                  ": {}\n\tmax support {}",
                  st.nequality, perceq, st.avgRowSupport, st.minRowSupport,
                  st.maxRowSupport);

   float percmatnnz =
       100.0 * (static_cast<double>(st.nnzmat) / (ncols * nrows));
   float percobjnnz = 100.0 * (static_cast<double>(st.nnzobj) / ncols);
   Message::print("Coefficients:");
   Message::print(
       "\tmatrix nnz: {} ({:0.1f}%)\n\tobjective nnz: {} ({:0.1f}%)",
       st.nnzmat, percmatnnz, st.nnzobj, percobjnnz);
   Message::print("");
}
