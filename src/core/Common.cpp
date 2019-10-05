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

         if (coef > 0.0)
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

int
updateSolActivity(std::vector<double>& activities, VectorView colview,
                  const std::vector<double>& lhs,
                  const std::vector<double>& rhs, double deltaVal,
                  std::vector<int>& violatedRows,
                  dynamic_bitset<>& isviolated)
{
   assert(isviolated.size() == lhs.size());
   int nviolated = 0;
   auto [coefs, indices, size] = colview;

   for (int i = 0; i < size; ++i)
   {
      const int row = indices[i];
      const double coef = coefs[i];

      activities[row] += coef * deltaVal;

      if (!Num::isFeasGE(activities[row], lhs[row]) ||
          !Num::isFeasLE(activities[row], rhs[row]))
      {
         if (!isviolated[row])
         {
            violatedRows.push_back(row);
            ++nviolated;
            isviolated[row] = true;
         }
      }
      else if (isviolated[row])
      {
         --nviolated;
         isviolated[row] = false;
      }
   }

   return nviolated;
}

std::vector<int>
getFractional(const std::vector<double>& sol,
              const dynamic_bitset<>& integer)
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
getIdentity(int ncols)
{
   std::vector<int> identity(ncols);

   for (int i = 0; i < ncols; ++i)
      identity[i] = i;
   return identity;
}

bool
hasZeroLockRounding(const std::vector<int>& downLocks,
                    const std::vector<int>& upLocks,
                    const std::vector<int>& fractional)
{
   for (auto col : fractional)
   {
      if (downLocks[col] && upLocks[col])
         return false;
   }
   return true;
}

bool
hasZeroLockRounding(const std::vector<double>& solution,
                    const std::vector<int>& downLocks,
                    const std::vector<int>& upLocks,
                    const dynamic_bitset<>& integer)
{
   for (size_t col = 0; col < solution.size(); ++col)
   {
      if (integer[col] && !Num::isIntegral(solution[col]))
      {
         if (downLocks[col] && upLocks[col])
            return false;
      }
   }
   return true;
}

double
zeroLockRound(std::vector<double>& solution,
              const std::vector<int>& downLocks,
              const std::vector<int>& fractional,
              const std::vector<double>& objective)
{
   double objdiff = 0.0;

   for (int col : fractional)
   {
      double oldval = solution[col];

      if (!downLocks[col])
         solution[col] = Num::floor(solution[col]);
      else
         solution[col] = Num::ceil(solution[col]);

      objdiff += (solution[col] - oldval) * objective[col];
   }

   return objdiff;
}

std::optional<std::pair<std::vector<double>, double>>
minLockRound(const MIP& mip, const std::vector<double>& solution,
             double obj, const std::vector<int>& fractional)
{
   int nrows = mip.getNRows();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();

   const auto& downLocks = mip.getDownLocks();
   const auto& upLocks = mip.getUpLocks();

   const auto& objective = mip.getObj();

   auto localsolution = solution;
   double localobj = obj;

   for (size_t i = 0; i < fractional.size(); ++i)
   {
      int col = fractional[i];

      double oldsolval = localsolution[col];
      if (downLocks[col] < upLocks[col])
         localsolution[col] = Num::floor(solution[col]);
      else
         localsolution[col] = Num::ceil(solution[col]);

      localobj += (localsolution[col] - oldsolval) * objective[col];
   }

   for (int row = 0; row < nrows; ++row)
   {
      auto [coefs, indices, size] = mip.getRow(row);
      double activity = 0.0;

      for (int i = 0; i < size; ++i)
         activity += coefs[i] * localsolution[indices[i]];

      if (!Num::isFeasLE(activity, rhs[row]) ||
          !Num::isFeasGE(activity, lhs[row]))
         return {};

      assert(Num::isFeasLE(activity, rhs[row]) &&
             Num::isFeasGE(activity, lhs[row]));
   }

   return std::make_pair(std::move(localsolution), localobj);
}

void
maxOutSolution(const MIP& mip, std::vector<double>& solution,
               const std::vector<double>& activity)
{
   const auto& objective = mip.getObj();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();
   const auto& integer = mip.getInteger();

   auto localact = activity;

   for (size_t col = 0; col < solution.size(); ++col)
   {
      if (objective[col] == 0.0)
         continue;

      int objsens = 1 - 2 * (objective[col] < 0.0);
      auto [colcoefs, colindices, colsize] = mip.getCol(col);

      double absdelta = objsens * Num::infval;
      for (int id = 0; id < colsize; ++id)
      {
         int row = colindices[id];
         double coef = colcoefs[id];

         int direction = objsens * (1 - 2 * (coef < 0.0));

         if (direction == 1)
         {
            if (!Num::isInf(rhs[row]))
               absdelta = std::min(absdelta, (rhs[row] - localact[row]) /
                                                 std::fabs(coef));
         }
         else
         {
            assert(direction == -1);
            if (!Num::isMinusInf(lhs[row]))
               absdelta = std::min(absdelta, (localact[row] - lhs[row]) /
                                                 std::fabs(coef));
         }
      }

      if (absdelta < 1e-6)
         continue;

      double oldval = solution[col];
      solution[col] += objsens * absdelta;

      if (integer[col])
      {
         if (objsens == 1)
            solution[col] = Num::floor(solution[col]);
         else
            solution[col] = Num::ceil(solution[col]);
      }

      // update the activity
      for (int id = 0; id < colsize; ++id)
      {
         int row = colindices[id];
         double coef = colcoefs[id];

         localact[row] += coef * (solution[col] - oldval);
      }
   }
}

void
roundFeasIntegers(std::vector<double>& sol,
                  const dynamic_bitset<>& integer)
{
   assert(sol.size() == integer.size());

   for (size_t i = 0; i < sol.size(); ++i)
   {
      if (integer[i] && Num::isFeasInt(sol[i]))
         sol[i] = Num::round(sol[i]);
   }
}

