#include "MinLockRounding.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

void
MinLockRounding::search(const MIP& mip, const std::vector<double>& lb,
                        const std::vector<double>& ub,
                        const std::vector<Activity>&,
                        const LPResult& result,
                        const std::vector<double>& solAct,
                        const std::vector<int>& fractional,
                        std::shared_ptr<const LPSolver> lpsolver,
                        SolutionPool& pool)
{
   int nrows = mip.getNRows();
   int ncols = mip.getNCols();
   int ncont = mip.getStats().ncont;

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
                      assert(left < ncols);
                      assert(right < ncols);
                      return mip.getColSize(left) < mip.getColSize(right);
                   });
         break;
      default:
         assert(0);
      }

      // TODO try std set instead
      std::vector<int> violatedRows;
      dynamic_bitset<> isviolated(nrows, false);
      violatedRows.reserve(nrows / 4);
      for (size_t i = 0; i < fractional.size(); ++i)
      {
         violatedRows.clear();
         isviolated.reset();
         int nviolated = 0;

         int col = fracPermutation[i];
         assert(col < st.nbin + st.nint);

         // if it was set to an integer value in
         // the correction phase
         if (Num::isIntegral(solution[col]))
            continue;

         // round
         double oldval = solution[col];
         if (downLocks[col] < upLocks[col])
            solution[col] = Num::floor(solution[col]);
         else
            solution[col] = Num::ceil(solution[col]);

         nviolated += updateSolActivity(solActivity, mip.getCol(col), lhs,
                                        rhs, solution[col] - oldval,
                                        violatedRows, isviolated);

         auto lpviolated = getNViolated<double, true>;
         assert(nviolated == lpviolated(mip, solution, 1e-9, 1e-6));

         if (nviolated == 0)
            continue;

         Message::debug_details(
             "Round: {} rows violated after rounding col {} from {} -> {}",
             nviolated, col, oldval, solution[col]);

         // it's possible to have a cycling change of values of continuous
         // variables so we limit the number of times they can change
         int ncontchanges = 0;
         for (size_t j = 0;
              j < violatedRows.size() && ncontchanges < 2 * ncont; ++j)
         {
            int row = violatedRows[j];
            assert(row < nrows);

            // if it was corrected after being added
            // to the violated list
            if (!isviolated[row])
               continue;

            assert(!Num::isFeasGE(solActivity[row], lhs[row]) ||
                   !Num::isFeasLE(solActivity[row], rhs[row]));

            Message::debug_details(
                "Round: trying to correct row {}: {} <= {} <= {}", row,
                lhs[row], solActivity[row], rhs[row]);

            auto [rowcoefs, rowindices, rowsize] = mip.getRow(row);

            bool row_corrected = false;

            for (int k = 0; k < rowsize; ++k)
            {
               int ncol = rowindices[k];
               double ncoef = rowcoefs[k];
               double oldnval = solution[ncol];

               // check if fractional
               if (ncol < st.nbin + st.nint &&
                   Num::isIntegral(solution[ncol]))
                  continue;

               if (!Num::isFeasGE(solActivity[row], lhs[row]))
               {
                  if (ncol < st.nbin + st.nint)
                  {
                     if (ncoef > 0.0)
                        solution[ncol] = Num::ceil(solution[ncol]);
                     else
                        solution[ncol] = Num::floor(solution[ncol]);
                  }
                  else
                  {
                     if (ncoef > 0.0)
                        solution[ncol] +=
                            std::min((lhs[row] - solActivity[row]) / ncoef,
                                     ub[ncol] - oldnval);
                     else
                        solution[ncol] +=
                            std::max((lhs[row] - solActivity[row]) / ncoef,
                                     lb[ncol] - oldnval);
                  }
               }
               else
               {
                  assert(!Num::isFeasLE(solActivity[row], rhs[row]));

                  if (ncol < st.nbin + st.nint)
                  {
                     if (ncoef > 0.0)
                        solution[ncol] = Num::floor(solution[ncol]);
                     else
                        solution[ncol] = Num::ceil(solution[ncol]);
                  }
                  else
                  {
                     if (ncoef > 0.0)
                        solution[ncol] +=
                            std::max((rhs[row] - solActivity[row]) / ncoef,
                                     lb[ncol] - oldnval);
                     else
                        solution[ncol] +=
                            std::min((rhs[row] - solActivity[row]) / ncoef,
                                     ub[ncol] - oldnval);
                  }
               }

               if (!Num::isFeasEQ(solution[ncol], oldnval))
               {
                  Message::debug_details(
                      "Round: changed col {} (int?: {}, "
                      "coef {})  value from {} -> {}",
                      ncol, col < st.nbin + st.nint, ncoef, oldnval,
                      solution[ncol]);

                  if (ncol >= st.nbin + st.nint)
                     ++ncontchanges;

                  nviolated += updateSolActivity(
                      solActivity, mip.getCol(ncol), lhs, rhs,
                      solution[ncol] - oldnval, violatedRows, isviolated);
                  Message::debug_details(
                      "Round: number of rows violated after col change {}",
                      nviolated);
               }

               auto act = computeSolActivities(mip, solution);
               for (int lrow = 0; lrow < nrows; ++lrow)
                  assert(std::fabs(solActivity[lrow] - act[lrow]) < 1e-6);

               // if the row was corrected
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
                "Round: infeasible, nviolated {} after fixing {} cols",
                nviolated, i + 1);
            feasible = false;
            break;
         }
      }

      if (feasible)
      {
         Message::debug("Round: feasible");

         if (ncont == 0)
         {
            Message::debug("Round: 0 cont");
            assert(checkFeasibility<double>(mip, solution));

            double cost = 0.0;
            for (int i = 0; i < ncols; ++i)
               cost += objective[i] * solution[i];

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
               Message::debug("Round: lp sol feasible");
               assert(checkFeasibility<double>(
                   mip, local_result.primalSolution, 1e-6, 1e-6));
               pool.add(std::move(local_result.primalSolution),
                        local_result.obj);
            }
            else if (local_result.status == LPResult::INFEASIBLE)
               Message::debug("Round: lp sol infeasible");
            else
               assert(0);
         }
      }
   } while (++ordering < 2 && feasible);
}
