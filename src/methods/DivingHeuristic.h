#ifndef DIVING_HEUR_HPP
#define DIVING_HEUR_HPP

#include "core/Heuristic.h"
#include "core/MIP.h"
#include "core/Propagation.h"
#include "io/Message.h"

static std::string
operator+(std::string_view str1, std::string_view str2)
{
   return std::string(str1) + std::string(str2);
}

// TODO better stopping criteria

template <typename SELECTION>
class DivingHeuristic : public HeuristicMethod
{
 public:
   DivingHeuristic() : HeuristicMethod(SELECTION::name + "Diving") {}

   ~DivingHeuristic() override = default;

   void search(const MIP& mip, const std::vector<double>& lb,
               const std::vector<double>& ub,
               const std::vector<Activity>& activities,
               const LPResult& result, const std::vector<double>&,
               const std::vector<int>& fractional,
               std::shared_ptr<const LPSolver> solver,
               SolutionPool& pool) override
   {
      int ncols = mip.getNCols();
      const auto& integer = mip.getInteger();
      const auto& objective = mip.getObj();
      const auto& downLocks = mip.getDownLocks();
      const auto& upLocks = mip.getUpLocks();

      auto local_activities = activities;
      auto locallb = lb;
      auto localub = ub;
      auto localfrac = fractional;
      auto localsol = result.primalSolution;
      double localobj = -1.0;
      std::unique_ptr<LPSolver> localsolver = solver->clone();

      auto heur_name = SELECTION::name + "Diving";

      bool feasible = true;
      int iter = 0;

      std::vector<int> buffer;

      bool hasZeroLockFractionals = false;
      int prev_simplex_iter = -1;
      bool limit_reached = false;
      do
      {
         ++iter;
         prev_simplex_iter = result.niter;

         auto [varToFix, direction, nFrac] =
             SELECTION::select(mip, locallb, localub, localsol);

         Message::debug_details("{}: iter {}, nFrac {}", heur_name, iter,
                                nFrac);

         if (varToFix < 0)
         {
            if (nFrac > 0)
               hasZeroLockFractionals = true;
            else
               assert(checkFeasibility<double>(mip, localsol, 1e-9, 1e-6));

            break;
         }

         double oldlb = locallb[varToFix];
         double oldub = localub[varToFix];

         if (Num::isMinusInf(locallb[varToFix]) ||
             Num::isInf(localub[varToFix]))
         {
            locallb[varToFix] = Num::floor(localsol[varToFix]);
            localub[varToFix] = Num::ceil(localsol[varToFix]);
         }

         if (direction == 1)
            locallb[varToFix] = localub[varToFix];
         else
         {
            assert(direction == -1);
            localub[varToFix] = locallb[varToFix];
         }

         bool status = propagate_get_changed_cols(
             mip, locallb, localub, local_activities, varToFix, oldlb,
             oldub, buffer);

         Message::debug_details(
             "{}: fixed var {} -> {}, propagation: {} changed {}",
             heur_name, varToFix, locallb[varToFix], status,
             buffer.size());

         if (!status)
         {
            feasible = false;
            break;
         }

         assert(locallb[varToFix] == localub[varToFix]);

         for (auto col : buffer)
            localsolver->changeBounds(col, locallb[col], localub[col]);
         buffer.clear();

         auto local_result = localsolver->solve(Algorithm::DUAL);

         if (local_result.status != LPResult::OPTIMAL)
            feasible = false;
         else
         {
            localsol = std::move(local_result.primalSolution);
            localobj = local_result.obj;

            roundFeasIntegers(localsol, integer);

            auto checklpFeas = checkFeasibility<double, true>;
            assert(checklpFeas(mip, localsol, 1e-9, 1e-6));
            assert(Num::isFeasEQ(localsol[varToFix], locallb[varToFix]));
         }

         if (iter > ncols * iter_per_col_max ||
             result.niter > simplex_iter_growth_max * prev_simplex_iter)
            limit_reached = true;

      } while (feasible && !limit_reached);

      if (feasible && !limit_reached && hasZeroLockFractionals)
      {
         localobj = 0.0;

         for (int col = 0; col < ncols; ++col)
         {
            localobj += objective[col] * localsol[col];
            if (!integer[col] || Num::isIntegral(localsol[col]))
               continue;

            double oldval = localsol[col];
            if (downLocks[col] == 0)
               localsol[col] = Num::floor(localsol[col]);
            else
            {
               assert(upLocks[col] == 0);
               localsol[col] = Num::ceil(localsol[col]);
            }

            localobj += objective[col] * (localsol[col] - oldval);
         }

         Message::debug("{}: found solution val {}", heur_name, localobj);
         pool.add(std::move(localsol), localobj);
      }
      else if (feasible && !limit_reached)
      {
         Message::debug("{}: found solution val {}", heur_name, localobj);
         assert(checkFeasibility<double>(mip, localsol, 1e-6, 1e-6));
         pool.add(std::move(localsol), localobj);
      }
      else
      {
         Message::debug("{}: infeasible after {} iterations", heur_name,
                        iter);
      }
   }

 private:
   constexpr static double iter_per_col_max = 0.25;
   constexpr static double simplex_iter_growth_max = 1.05;
};

#endif
