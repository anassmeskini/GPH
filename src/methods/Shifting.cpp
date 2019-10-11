#include "Shifting.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

void
Shifting::search(const MIP& mip, const std::vector<double>& lb,
                 const std::vector<double>& ub,
                 const std::vector<Activity>&, const LPResult& result,
                 const std::vector<double>& solAct,
                 const std::vector<int>& fractional,
                 std::shared_ptr<const LPSolver> lpsolver,
                 SolutionPool& pool)
{
   int nrows = mip.getNRows();
   int ncols = mip.getNCols();

   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();
   const auto& upLocks = mip.getUpLocks();
   const auto& downLocks = mip.getDownLocks();
   auto st = mip.getStats();
   const auto& objective = mip.getObj();

   std::unique_ptr<LPSolver> localsolver;

   int ordering = 0;
   bool feasible = true;
   do
   {
      auto solActivity = solAct;
      auto solution = result.primalSolution;
      auto fracPermutation = fractional;

      switch (ordering)
      {
      case 0:
         std::sort(std::begin(fracPermutation), std::end(fracPermutation),
                   [&](int left, int right) {
                      return std::min(downLocks[left], upLocks[left]) <
                             std::min(downLocks[right], upLocks[right]);
                   });
         break;
      case 1:
         std::sort(std::begin(fracPermutation), std::end(fracPermutation),
                   [&](int left, int right) {
                      return std::max(downLocks[left], upLocks[left]) <
                             std::max(downLocks[right], upLocks[right]);
                   });
         break;
      case 2:
         std::sort(std::begin(fracPermutation), std::end(fracPermutation),
                   [&](int left, int right) {
                      assert(left < ncols);
                      assert(right < ncols);
                      return mip.getColSize(left) < mip.getColSize(right);
                   });
         break;
      case 3:
         std::sort(std::begin(fracPermutation), std::end(fracPermutation),
                   [&](int left, int right) {
                      assert(left < ncols);
                      assert(right < ncols);
                      return mip.getColSize(left) > mip.getColSize(right);
                   });
         break;
      default:
         assert(0);
      }

      for (int i = 0; i < static_cast<int>(fractional.size()); ++i)
      {
         int nviolated = 0;
         std::vector<int> violatedRows;
         dynamic_bitset<> isviolated(nrows, false);
         violatedRows.reserve(nrows);

         int col = fracPermutation[i];

         assert(col < st.nbin + st.nint);
         if (Num::isIntegral(solution[col]))
            continue;

         double oldval = solution[col];
         if (downLocks[col] < upLocks[col])
            solution[col] = Num::floor(solution[col]);
         else
            solution[col] = Num::ceil(solution[col]);

         nviolated += updateSolActivity(solActivity, mip.getCol(col), lhs,
                                        rhs, solution[col] - oldval,
                                        violatedRows, isviolated);

         if (nviolated == 0)
            continue;

         Message::debug_details(
             "Shif: {} rows violated after rouding col {} from {} -> {}",
             nviolated, col, oldval, solution[col]);

         // it's possible to have a cycling change of values of continuous
         // variables so we limit the number of times they can change
         int nchanges = 0;
         for (size_t j = 0;
              j < violatedRows.size() && nchanges <= 50 * ncols; ++j)
         {
            int row = violatedRows[j];
            assert(row < nrows);

            if (!isviolated[row])
               continue;

            assert(!Num::isFeasGE(solActivity[row], lhs[row]) ||
                   !Num::isFeasLE(solActivity[row], rhs[row]));

            Message::debug_details(
                "Shif: trying to correct row {}: {} <= {} <= {}", row,
                lhs[row], solActivity[row], rhs[row]);

            auto [rowcoefs, rowindices, rowsize] = mip.getRow(row);

            violatedRows.clear();
            bool row_corrected = false;

            // first try: only fix integer fractional variables
            for (int k = 0; k < rowsize; ++k)
            {
               int ncol = rowindices[k];
               double ncoef = rowcoefs[k];
               double oldnval = solution[ncol];

               // skip non fractional and cont variables
               if (ncol >= st.nbin + st.nint ||
                   Num::isIntegral(solution[ncol]) || ncol == col)
                  continue;

               if (!Num::isFeasGE(solActivity[row], lhs[row]))
               {
                  if (ncoef > 0.0)
                     solution[ncol] = Num::ceil(solution[ncol]);
                  else
                     solution[ncol] = Num::floor(solution[ncol]);
               }
               else
               {
                  assert(!Num::isFeasLE(solActivity[row], rhs[row]));

                  if (ncoef > 0.0)
                     solution[ncol] = Num::floor(solution[ncol]);
                  else
                     solution[ncol] = Num::ceil(solution[ncol]);
               }

               if (std::fabs(solution[ncol] - oldnval) > 1e-6)
               {
                  Message::debug_details("Shif: changed int col {} (coef "
                                         "{})  value from {} -> {}",
                                         ncol, ncoef, oldnval,
                                         solution[ncol]);

                  nviolated += updateSolActivity(
                      solActivity, mip.getCol(ncol), lhs, rhs,
                      solution[ncol] - oldnval, violatedRows, isviolated);
                  Message::debug_details(
                      "Shif: number of rows violated after col change {}",
                      nviolated);
               }

               auto act = computeSolActivities(mip, solution);
               for (int lrow = 0; lrow < nrows; ++lrow)
                  assert(std::fabs(solActivity[lrow] - act[lrow]) < 1e-6);

               if (Num::isFeasGE(solActivity[row], lhs[row]) &&
                   Num::isFeasLE(solActivity[row], rhs[row]))
               {
                  row_corrected = true;
                  break;
               }
            }

            if (row_corrected)
               continue;

            Message::debug_details("Shif: Trying shifting");

            for (int k = 0; k < rowsize; ++k)
            {
               int ncol = rowindices[k];
               double ncoef = rowcoefs[k];
               double oldnval = solution[ncol];

               if (!Num::isFeasGE(solActivity[row], lhs[row]))
               {
                  if (ncoef > 0.0)
                  {
                     solution[ncol] +=
                         std::min((lhs[row] - solActivity[row]) / ncoef,
                                  ub[ncol] - oldnval);

                     if (ncol < st.nbin + st.nint)
                        solution[ncol] = Num::ceil(solution[ncol]);
                  }
                  else
                  {
                     solution[ncol] +=
                         std::max((lhs[row] - solActivity[row]) / ncoef,
                                  lb[ncol] - oldnval);

                     if (ncol < st.nbin + st.nint)
                        solution[ncol] = Num::floor(solution[ncol]);
                  }
               }
               else
               {
                  assert(!Num::isFeasLE(solActivity[row], rhs[row]));

                  if (ncoef > 0.0)
                  {
                     solution[ncol] +=
                         std::max((rhs[row] - solActivity[row]) / ncoef,
                                  lb[ncol] - oldnval);

                     if (ncol < st.nbin + st.nint)
                        solution[ncol] = Num::floor(solution[ncol]);
                  }
                  else
                  {
                     solution[ncol] +=
                         std::min((rhs[row] - solActivity[row]) / ncoef,
                                  ub[ncol] - oldnval);

                     if (ncol < st.nbin + st.nint)
                        solution[ncol] = Num::ceil(solution[ncol]);
                  }
               }

               if (std::fabs(solution[ncol] - oldnval) > 1e-6)
               {
                  Message::debug_details(
                      "Shif: changed col {} (int? {}, coef "
                      "{})  value from {} -> {}",
                      ncol, ncol < st.nbin + st.nint, ncoef, oldnval,
                      solution[ncol]);

                  nviolated += updateSolActivity(
                      solActivity, mip.getCol(ncol), lhs, rhs,
                      solution[ncol] - oldnval, violatedRows, isviolated);

                  ++nchanges;
                  Message::debug_details(
                      "Shif: number of rows violated after col change {}",
                      nviolated);
               }

#ifndef NDEBUG
               auto act = computeSolActivities(mip, solution);
               for (int lrow = 0; lrow < nrows; ++lrow)
                  assert(std::fabs(solActivity[lrow] - act[lrow]) < 1e-6);
#endif

               if (Num::isFeasGE(solActivity[row], lhs[row]) &&
                   Num::isFeasLE(solActivity[row], rhs[row]))
               {
                  row_corrected = true;
                  break;
               }
            }

            if (!row_corrected)
               break;
         }

         if (nviolated > 0)
         {
            Message::debug(
                "Shif: infeasible, nviolated {} after fixing {} cols",
                nviolated, i + 1);
            feasible = false;
            break;
         }
      }

      if (feasible)
      {
         Message::debug("Shif: feasible");

         for (int row = 0; row < nrows; ++row)
            assert(solActivity[row] >= lhs[row] - 1e-6 &&
                   solActivity[row] <= rhs[row] + 1e-6);

         if (mip.getStats().ncont == 0)
         {
            Message::debug("Shif: 0 cont");

            double cost = 0.0;
            for (int i = 0; i < ncols; ++i)
               cost += objective[i] * solution[i];

            assert(checkFeasibility<double>(mip, solution));
            pool.add(std::move(solution), cost);
         }
         else
         {
            if (!localsolver)
               localsolver = lpsolver->clone();

            for (int col = 0; col < st.nbin + st.nint; ++col)
            {
               assert(Num::isIntegral(solution[col]));
               localsolver->changeBounds(col, solution[col],
                                         solution[col]);
            }

            auto local_result = localsolver->solve(Algorithm::DUAL);
            if (local_result.status == LPResult::OPTIMAL)
            {
               Message::debug("Shif: lp sol feasible");
               assert(checkFeasibility<double>(
                   mip, local_result.primalSolution));
               pool.add(std::move(local_result.primalSolution),
                        local_result.obj);
            }
            else if (local_result.status == LPResult::INFEASIBLE)
               Message::debug("Shif: lp sol infeasible");
            else
               assert(0);
         }
      }

      ++ordering;
   } while (ordering < 4 && feasible);
}
