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
         std::vector<double>&& _obj, dynamic_bitset<>&& _integer,
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

   // move integer bitset
   integer = std::move(_integer);

   // move the objective
   objective = std::move(_obj);

   // now before moving the transposed matrix, we remove
   // slack variables (non-integer unit columns), and remember these
   // changes to adjust the objective after we have access to the row-major
   // matrix

   // first: slack's row, second: removed col,
   // third: slack's coefficient in the row
   // fourth: the rhs of the constraint
   using Tuple = std::tuple<int, int, double, double>;
   std::vector<Tuple> removedSlacks;

   int nslacks = 0;

   // TODO this is buggy on b1c1s1.mps
   /*for (int col = 0; col < ncols; ++col)
   {
      int colsize = rstartT[col + 1] - rstartT[col];

      if (colsize == 1 && !integer[col])
      {
         const int offset = rstartT[col];
         const int row = idxT[offset];
         const double slack_coef = coefsT[offset];

         assert(slack_coef != 0.0);

         double* coefs = coefsT.data() + offset;
         int* indices = idxT.data() + offset;
         int nnz = coefsT.size();

         if (lhs[row] == rhs[row])
         {
            if (objective[col] != 0.0)
               removedSlacks.emplace_back(row, col, slack_coef, rhs[row]);
            ++nslacks;

            // replace by lower bound and change constype
            std::memmove(coefs, coefs + 1,
                         sizeof(double) * (nnz - (offset + 1)));
            std::memmove(indices, indices + 1,
                         sizeof(int) * (nnz - (offset + 1)));
            --nnz;

            coefsT.resize(nnz);
            idxT.resize(nnz);

            for (int i = col + 1; i < ncols + 1; ++i)
               --rstartT[i];
            --rowSize[row];

            assert(rstartT.back() == static_cast<int>(coefsT.size()));

            if (slack_coef > 0.0)
            {
               if (Num::isMinusInf(lb[col]))
                  rhs[row] = Num::infval;
               else
                  rhs[row] -= slack_coef * lb[col];

               if (Num::isInf(ub[col]))
                  lhs[row] = -Num::infval;
               else
                  lhs[row] -= slack_coef * ub[col];
            }
            else
            {
               if (Num::isInf(ub[col]))
                  rhs[row] = Num::infval;
               else
                  rhs[row] -= slack_coef * ub[col];

               if (Num::isMinusInf(lb[col]))
                  lhs[row] = -Num::infval;
               else
                  lhs[row] -= slack_coef * lb[col];
            }

            assert(lhs[row] != rhs[row] || lb[col] == ub[col]);
         }
      }
   }*/

   Message::debug("Removed {} slacks", nslacks);

   // fill the column major matrix
   constMatrixT.coefficients = std::move(coefsT);
   constMatrixT.indices = std::move(idxT);
   constMatrixT.rowStart = std::move(rstartT);
   constMatrixT.ncols = nrows;
   constMatrixT.nrows = ncols;
   stats.ncols = ncols;
   stats.nrows = nrows;

   // get the row-major matrix by transposing
   constMatrix = transpose(constMatrixT, rowSize);

   // adjust the objective for the removed slacks
   // s = (rhs - Ax)/b
   /*for (auto tuple : removedSlacks)
   {
      int row = std::get<0>(tuple);
      int slack = std::get<1>(tuple);
      double slack_coef = std::get<2>(tuple);
      double rhs = std::get<3>(tuple);

      double slack_obj = objective[slack];
      auto [coefs, indices, size] = getRow(row);

      assert(slack_obj != 0.0);
      assert(std::all_of(indices, indices + size,
                         [=](int col) -> bool { return col != slack; }));

      objoffset += slack_obj * rhs / slack_coef;

      for (int i = 0; i < size; ++i)
      {
         const int col = indices[i];
         const double coef = coefs[i];

         objective[col] -= slack_obj * coef / slack_coef;
      }
   }*/

   /*sortRows(constMatrix, [&](int left, int right) {
      return integer[left] && !integer[right];
   });*/

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
            binary.push_back(col);
         else
            generalInt.push_back(col);
      }
      else
         continuous.push_back(col);
   }

   stats.nbin = static_cast<int>(binary.size());
   stats.nint = static_cast<int>(generalInt.size());
   stats.ncont = static_cast<int>(continuous.size());

   stats.avgRowSupport /= nrows;
   stats.avgColSupport /= ncols;
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

   integer = std::move(other.integer);

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
