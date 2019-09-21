#ifndef MAKE_MIP_HPP
#define MAKE_MIP_HPP

#define UNIT_TEST

#include "core/MIP.h"
#include <vector>

inline MIP
make_mip(const std::vector<double>& coefficients, int ncols,
         const std::vector<double>& lhs, const std::vector<double>& rhs,
         const std::vector<double>& lb, const std::vector<double>& ub,
         const dynamic_bitset<>& integer,
         const std::vector<double>& objective)
{
   // 2*x1 + 3*x2 + 4*x3 + 1*x5 <= 3
   // all binary
   MIP mip;

   int nrows = static_cast<int>(coefficients.size()) / ncols;
   mip.lhs = lhs;
   mip.rhs = rhs;
   mip.lb = lb;
   mip.ub = ub;
   mip.objective = objective;
   mip.integer = integer;

   std::vector<int> colsize(ncols, 0);

   mip.constMatrix.nrows = nrows;
   mip.constMatrix.ncols = ncols;
   mip.constMatrix.rowStart.push_back(0);

   for (int row = 0; row < nrows; ++row)
   {
      int nnz = 0;
      for (int col = 0; col < ncols; ++col)
      {
         int index = ncols * row + col;
         if (coefficients[index] != 0.0)
         {
            ++nnz;
            ++colsize[col];
            mip.constMatrix.coefficients.push_back(coefficients[index]);
            mip.constMatrix.indices.push_back(col);
         }
      }
      mip.constMatrix.rowStart.push_back(nnz);
   }

   mip.constMatrixT = MIP::transpose(mip.constMatrix, colsize);

   return mip;
}

inline double
get_coef(const MIP& mip, std::pair<int, int> indices)
{
   int row = indices.first;
   int col = indices.second;

   double coef = 0.0;
   for (int index = mip.constMatrix.rowStart[row];
        index < mip.constMatrix.rowStart[row + 1]; ++index)
   {
      if (mip.constMatrix.indices[index] == col)
      {
         coef = mip.constMatrix.coefficients[index];
         break;
      }
   }

   return coef;
}

inline void
change_coef(MIP& mip, std::pair<int, int> indices, double newcoef)
{
   int row = indices.first;
   int col = indices.second;

   for (int index = mip.constMatrix.rowStart[row];
        index < mip.constMatrix.rowStart[row + 1]; ++index)
   {
      if (mip.constMatrix.indices[index] == col)
      {
         mip.constMatrix.coefficients[index] = newcoef;
         break;
      }
   }

   for (int index = mip.constMatrixT.rowStart[col];
        index < mip.constMatrixT.rowStart[col + 1]; ++index)
   {
      if (mip.constMatrixT.indices[index] == row)
      {
         mip.constMatrixT.coefficients[index] = newcoef;
         break;
      }
   }
}

#endif
