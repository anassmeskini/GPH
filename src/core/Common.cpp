#include "Common.h"
#include "Numerics.h"
#include "SparseMatrix.h"

std::vector<Activity>
computeActivities(const MIP& mip)
{
   int nrows = mip.getNRows();
   std::vector<Activity> activities(nrows);

   const auto& lb = mip.getLB();
   const auto& ub = mip.getUB();

   for (int row = 0; row < nrows; ++row)
   {
      auto [coefs, indices, rowsize] = mip.getRow(row);

      for (int colid = 0; colid < rowsize; ++colid)
      {
         const double coef = coefs[colid];
         const int col = indices[colid];

         if (Num::greater(coef, 0.0))
         {
            if (!Num::isMinusInf(lb[col]))
               activities[row].min += lb[col] * coef;
            else
               activities[row].ninfmin++;

            if (!Num::isInf(ub[col]))
               activities[row].max += ub[col] * coef;
            else
               activities[row].ninfmax++;
         }
         else
         {
            if (!Num::isMinusInf(lb[col]))
               activities[row].max += lb[col] * coef;
            else
               activities[row].ninfmax++;

            if (!Num::isInf(ub[col]))
               activities[row].min += ub[col] * coef;
            else
               activities[row].ninfmin++;
         }
      }
   }

   return activities;
}

std::vector<double>
computeSolActivities(const MIP& mip, const std::vector<double>& sol)
{
   int nrows = mip.getNRows();

   std::vector<double> activities(nrows, 0.0);

   for (int row = 0; row < nrows; ++row)
   {
      auto [coefs, indices, rowsize] = mip.getRow(row);

      for (int colid = 0; colid < rowsize; ++colid)
      {
         const double coef = coefs[colid];
         const int col = indices[colid];

         activities[row] += sol[col] * coef;
      }
   }

   return activities;
}

std::vector<int>
getFractional(const std::vector<double>& sol, const dynamic_bitset<>& integer)
{
   assert(sol.size() == integer.size());
   int ncols = sol.size();
   std::vector<int> fractional;
   fractional.reserve(ncols);

   for (int col = 0; col < ncols; ++col)
   {
      if (integer[col] && !Num::isIntegral(sol[col]))
         fractional.push_back(col);
   }

   return fractional;
}

std::vector<int>
getIndentity(int ncols)
{
   std::vector<int> identity(ncols);

   for (int i = 0; i < ncols; ++i)
      identity[i] = i;
   return identity;
}
