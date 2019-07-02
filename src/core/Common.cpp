#include "Common.h"
#include "Numerics.h"
#include "SparseMatrix.h"

std::vector<Activity>
computeActivities(const MIP<double>& mip)
{
   size_t nrows = mip.getNRows();
   std::vector<Activity> activities(nrows);

   const auto& lb = mip.getLB();
   const auto& ub = mip.getUB();

   for (size_t row = 0; row < nrows; ++row)
   {
      auto [coefs, indices, rowsize] = mip.getRow(row);

      for (size_t colid = 0; colid < rowsize; ++colid)
      {
         const double coef = coefs[colid];
         const size_t col = indices[colid];

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
computeSolActivities(const MIP<double>& mip, const std::vector<double>& sol)
{
   size_t nrows = mip.getNRows();
   std::vector<double> activities(nrows, 0.0);

   for (size_t row = 0; row < nrows; ++row)
   {
      auto [coefs, indices, rowsize] = mip.getRow(row);

      for (size_t colid = 0; colid < rowsize; ++colid)
      {
         const double coef = coefs[colid];
         const size_t col = indices[colid];

         activities[row] += sol[col] * coef;
      }
   }

   return activities;
}
