#ifndef DIVING_HEUR_HPP
#define DIVING_HEUR_HPP

#include "core/Heuristic.h"
#include "core/MIP.h"
#include "io/Message.h"

std::string
operator+(std::string_view str1, const char* str2)
{
   return std::string(str1) + str2;
}

std::string
operator+(std::string_view str1, const std::string& str2)
{
   return std::string(str1) + str2;
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
               std::shared_ptr<const LPSolver> solver, SolutionPool& pool)
   {
      int ncols = mip.getNCols();
      const auto& integer = mip.getInteger();
      const auto& objective = mip.getObj();
      const auto& downLocks = mip.getDownLocks();
      const auto& upLocks = mip.getUpLocks();

      // try solution = lower bound
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

      bool hasZeroLockFractionals = false;
      do
      {
         ++iter;
         auto [varToFix, direction, nFrac] =
             SELECTION::select(mip, localsol);

         Message::debug_details("{}: iter {}, nFrac {}", heur_name, iter,
                                nFrac);

         // TODO
         if (varToFix < 0)
         {
            if (nFrac > 0)
               hasZeroLockFractionals = true;
            break;
         }

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

         Message::debug_details("{}: fixed var {} -> {}", heur_name,
                                varToFix, locallb[varToFix]);

         --nFrac;

         // TODO propagate

         assert(locallb[varToFix] == localub[varToFix]);
         localsolver->changeBounds(varToFix, locallb[varToFix],
                                   localub[varToFix]);
         auto res = localsolver->solve(Algorithm::DUAL);

         if (res.status != LPResult::OPTIMAL)
            feasible = false;
         else
         {
            localsol = std::move(res.primalSolution);
            localobj = res.obj;

            assert(Num::isFeasEQ(localsol[varToFix], locallb[varToFix]));
         }

      } while (feasible && iter < ncols / 4);

      if (feasible && hasZeroLockFractionals)
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
      }

      if (feasible)
      {
         // pool.add(std::move(localsol), localobj);
         Message::debug("{}: found solution val {}", heur_name, localobj);
         pool.add(std::move(localsol), localobj);
      }
      else
      {
         Message::debug("{}: infeasible after {} iterations", heur_name,
                        iter);
      }
   }
};

#endif
