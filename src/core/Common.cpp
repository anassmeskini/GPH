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
getFractional(const std::vector<double>& sol, int ninteger)
{
   int ncols = sol.size();
   std::vector<int> fractional;
   fractional.reserve(ncols);

   for (int col = 0; col < ninteger; ++col)
   {
      if (!Num::isIntegral(sol[col]))
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
                    const std::vector<int>& upLocks, int ninteger)
{
   assert(ninteger <= static_cast<int>(solution.size()));
   for (int col = 0; col < ninteger; ++col)
   {
      if (!Num::isIntegral(solution[col]))
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
roundFeasIntegers(std::vector<double>& sol, int ninteger)
{
   for (int i = 0; i < ninteger; ++i)
   {
      if (Num::isFeasInt(sol[i]))
         sol[i] = Num::round(sol[i]);
   }
}

void
maxOutSolution(const MIP& mip, std::vector<double>& solution, double& cost)
{
   auto st = mip.getStats();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();
   const auto& lb = mip.getLB();
   const auto& ub = mip.getUB();
   const auto& objective = mip.getObj();

   std::vector<double> activities(st.nrows, 0.0);

   for (int row = 0; row < st.nrows; ++row)
   {
      auto [rowcoefs, rowindices, rowsize] = mip.getRow(row);

      for (int i = 0; i < rowsize; ++i)
      {
         int col = rowindices[i];
         double coef = rowcoefs[i];
         activities[row] += coef * solution[col];
      }

      assert(Num::isFeasGE(activities[row], lhs[row]));
      assert(Num::isFeasLE(activities[row], rhs[row]));
   }

   for (int col = 0; col < st.ncols; ++col)
   {
      if (objective[col] == 0)
         continue;

      auto [colcoefs, colindices, colsize] = mip.getCol(col);

      double max_abs_delta = Num::infval;
      for (int i = 0; i < colsize; ++i)
      {
         int row = colindices[i];
         double coef = colcoefs[i];

         if (objective[col] < 0)
         {
            if (coef > 0.0)
               max_abs_delta = std::min(
                   max_abs_delta, (rhs[row] - activities[row]) / coef);
            else
               max_abs_delta = std::min(
                   max_abs_delta, (lhs[row] - activities[row]) / coef);
         }
         else
         {
            if (coef > 0.0)
               max_abs_delta = std::min(
                   max_abs_delta, (activities[row] - lhs[row]) / coef);
            else
               max_abs_delta = std::min(
                   max_abs_delta, (activities[row] - rhs[row]) / coef);
         }

         assert(max_abs_delta >= 0.0);
      }

      if (col < st.nbin + st.nint)
         max_abs_delta = Num::floor(max_abs_delta);

      if (!Num::isFeasGE(max_abs_delta, 0.0))
         continue;

      double delta;

      if (objective[col] < 0.0)
         delta = std::min(max_abs_delta, ub[col] - solution[col]);
      else
         delta = -std::min(max_abs_delta, solution[col] - lb[col]);

      solution[col] += delta;
      cost += objective[col] * delta;

      assert(Num::isFeasGE(solution[col], lb[col]));
      assert(Num::isFeasLE(solution[col], ub[col]));

      for (int i = 0; i < colsize; ++i)
      {
         int row = colindices[i];
         double coef = colcoefs[i];
         activities[row] += coef * delta;
      }
   }
}
