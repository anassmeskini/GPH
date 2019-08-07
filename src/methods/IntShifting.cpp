#include "IntShifting.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

void IntShifting::search(const MIP &mip, const std::vector<double> &lb,
                         const std::vector<double> &ub,
                         const std::vector<Activity> &, const LPResult &result,
                         const std::vector<double> &solAct,
                         const std::vector<int> &fractional,
                         std::shared_ptr<const LPSolver> lpsolver,
                         SolutionPool &pool)
{
   // if there is no continuous variables
   // this heuristic es equal to shifting
   auto stats = mip.getStatistics();
   if (stats.ncont == 0)
      return;

   int nrows = mip.getNRows();
   int ncols = mip.getNCols();

   const auto &lhs = mip.getLHS();
   const auto &rhs = mip.getRHS();
   const auto &upLocks = mip.getUpLocks();
   const auto &downLocks = mip.getDownLocks();
   const auto &integer = mip.getInteger();

   std::unique_ptr<LPSolver> localsolver;

   auto locallhs = lhs;
   auto localrhs = rhs;

   // we replace the continuous variables with
   // their minimal contributions
   for (int row = 0; row < nrows; ++row)
   {
      auto [rowcoefs, rowindices, rowsize] = mip.getRow(row);

      for (int id = 0; id < rowsize; ++id)
      {
         int col = rowindices[id];
         int coef = rowcoefs[id];

         if (integer[col])
            continue;

         if (coef > 0.0)
         {
            if (!Num::isMinusInf(locallhs[row]))
            {
               if (!Num::isInf(ub[col]))
                  locallhs[row] = locallhs[row] - coef * ub[col];
               else
                  locallhs[row] = -Num::infval;
            }

            if (!Num::isInf(localrhs[row]))
            {
               if (!Num::isMinusInf(lb[col]))
                  localrhs[row] = localrhs[row] - coef * lb[col];
               else
                  localrhs[row] = Num::infval;
            }
         }
         else
         {
            if (!Num::isMinusInf(locallhs[row]))
            {
               if (!Num::isMinusInf(lb[col]))
                  locallhs[row] = locallhs[row] - coef * lb[col];
               else
                  locallhs[row] = -Num::infval;
            }

            if (!Num::isInf(localrhs[row]))
            {
               if (!Num::isInf(ub[col]))
                  localrhs[row] = localrhs[row] - coef * ub[col];
               else
                  localrhs[row] = Num::infval;
            }
         }
      }
   }

   // compute the solution activities that ignores continuous activities
   auto intSolAct = solAct;
   const auto &contVars = mip.getContVars();
   for (auto col : contVars)
   {
      auto [colcoefs, colindices, colsize] = mip.getCol(col);
      double solval = result.primalSolution[col];

      for (int id = 0; id < colsize; ++id)
      {
         int row = colindices[id];
         double coef = colcoefs[id];

         intSolAct[row] -= coef * solval;
      }
   }

   int ordering = 0;
   bool feasible = true;
   do
   {
      auto solActivity = intSolAct;
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

         assert(integer[col]);
         if (Num::isIntegral(solution[col]))
            continue;

         double oldval = solution[col];
         if (downLocks[col] < upLocks[col])
            solution[col] = Num::floor(solution[col]);
         else
            solution[col] = Num::ceil(solution[col]);

         nviolated += updateSolActivity(solActivity, mip.getCol(col), locallhs,
                                        localrhs, solution[col] - oldval,
                                        violatedRows, isviolated);

         if (nviolated == 0)
            continue;

         Message::debug_details(
             "ShifInt: {} rows violated after rouding col {} from {} -> {}",
             nviolated, col, oldval, solution[col]);

         // it's possible to have a cycling change of values of continuous
         // variables so we limit the number of times they can change
         int nchanges = 0;
         for (size_t j = 0; j < violatedRows.size() && nchanges <= 50 * ncols;
              ++j)
         {
            int row = violatedRows[j];
            assert(row < nrows);

            if (!isviolated[row])
               continue;

            assert(!Num::isFeasGE(solActivity[row], locallhs[row]) ||
                   !Num::isFeasLE(solActivity[row], localrhs[row]));

            Message::debug_details(
                "ShifInt: trying to correct row {}: {} <= {} <= {}", row,
                locallhs[row], solActivity[row], localrhs[row]);

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
               if (!integer[ncol] ||
                   (integer[ncol] && Num::isIntegral(solution[ncol])))
                  continue;

               if (!Num::isFeasGE(solActivity[row], locallhs[row]))
               {
                  if (ncoef > 0.0)
                     solution[ncol] = Num::ceil(solution[ncol]);
                  else
                     solution[ncol] = Num::floor(solution[ncol]);
               }
               else
               {
                  assert(!Num::isFeasLE(solActivity[row], localrhs[row]));

                  if (ncoef > 0.0)
                     solution[ncol] = Num::floor(solution[ncol]);
                  else
                     solution[ncol] = Num::ceil(solution[ncol]);
               }

               if (std::fabs(solution[ncol] - oldnval) > 1e-6)
               {
                  Message::debug_details("ShifInt: changed int col {} (coef "
                                         "{})  value from {} -> {}",
                                         ncol, ncoef, oldnval, solution[ncol]);

                  nviolated += updateSolActivity(
                      solActivity, mip.getCol(ncol), locallhs, localrhs,
                      solution[ncol] - oldnval, violatedRows, isviolated);
                  Message::debug_details(
                      "ShifInt: number of rows violated after col change {}",
                      nviolated);
               }

               auto act = computeSolActivities(mip, solution);
               for (int lrow = 0; lrow < nrows; ++lrow)
                  assert(std::fabs(solActivity[lrow] - act[lrow]) < 1e-6);

               if (Num::isFeasGE(solActivity[row], locallhs[row]) &&
                   Num::isFeasLE(solActivity[row], localrhs[row]))
               {
                  row_corrected = true;
                  break;
               }
            }

            if (row_corrected)
               continue;

            Message::debug_details("ShifInt: Trying shifting");

            for (int k = 0; k < rowsize; ++k)
            {
               int ncol = rowindices[k];
               double ncoef = rowcoefs[k];
               double oldnval = solution[ncol];

               if (!integer[ncol] || ncol == col)
                  continue;

               if (!Num::isFeasGE(solActivity[row], locallhs[row]))
               {
                  if (ncoef > 0.0)
                  {
                     solution[ncol] +=
                         std::min((locallhs[row] - solActivity[row]) / ncoef,
                                  ub[ncol] - oldnval);

                     solution[ncol] = Num::ceil(solution[ncol]);
                  }
                  else
                  {
                     solution[ncol] +=
                         std::max((locallhs[row] - solActivity[row]) / ncoef,
                                  lb[ncol] - oldnval);

                     solution[ncol] = Num::floor(solution[ncol]);
                  }
               }
               else
               {
                  assert(!Num::isFeasLE(solActivity[row], localrhs[row]));

                  if (ncoef > 0.0)
                  {
                     solution[ncol] +=
                         std::max((localrhs[row] - solActivity[row]) / ncoef,
                                  lb[ncol] - oldnval);

                     solution[ncol] = Num::floor(solution[ncol]);
                  }
                  else
                  {
                     solution[ncol] +=
                         std::min((localrhs[row] - solActivity[row]) / ncoef,
                                  ub[ncol] - oldnval);

                     solution[ncol] = Num::ceil(solution[ncol]);
                  }
               }

               if (std::fabs(solution[ncol] - oldnval) > 1e-6)
               {
                  Message::debug_details(
                      "ShifInt: changed col {} (int? {}, coef "
                      "{})  value from {} -> {}",
                      ncol, integer[ncol], ncoef, oldnval, solution[ncol]);

                  nviolated += updateSolActivity(
                      solActivity, mip.getCol(ncol), locallhs, localrhs,
                      solution[ncol] - oldnval, violatedRows, isviolated);

                  ++nchanges;
                  Message::debug_details(
                      "ShifInt: number of rows violated after col change {}",
                      nviolated);
               }

               auto act = computeSolActivities(mip, solution);
               for (int lrow = 0; lrow < nrows; ++lrow)
                  assert(std::fabs(solActivity[lrow] - act[lrow]) < 1e-6);

               if (Num::isFeasGE(solActivity[row], locallhs[row]) &&
                   Num::isFeasLE(solActivity[row], localrhs[row]))
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
                "ShifInt: infeasible, nviolated {} after fixing {} cols",
                nviolated, i + 1);
            feasible = false;
            break;
         }
      }

      if (feasible)
      {
         Message::debug("ShifInt: relaxed problem feasible");

         if (!localsolver)
            localsolver = lpsolver->clone();

         for (int col = 0; col < ncols; ++col)
         {
            if (integer[col])
            {
               assert(Num::isIntegral(solution[col]));
               localsolver->changeBounds(col, solution[col], solution[col]);
            }
         }

         auto local_result = localsolver->solve();
         if (local_result.status == LPResult::OPTIMAL)
         {
            Message::debug("ShifInt: lp sol feasible");

            pool.add(std::move(local_result.primalSolution), local_result.obj);
         }
         else if (local_result.status == LPResult::INFEASIBLE)
            Message::debug("ShifInt: lp sol infeasible");
         else
            assert(0);
      }

      ++ordering;
   } while (ordering < 4 && feasible);
}
