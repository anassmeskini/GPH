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

template <typename SELECTION>
class DivingHeuristic : public FeasibilityHeuristic
{
 public:
   DivingHeuristic() : FeasibilityHeuristic(SELECTION::name + "Diving") {}

   ~DivingHeuristic() override = default;

   void
   setParam(const std::string& param,
            const std::variant<std::string, int, double>& value) override
   {
      if (param == "propagate")
         propagate = static_cast<bool>(std::get<int>(value));
      else if (param == "backtrack")
         backtrack = static_cast<bool>(std::get<int>(value));
      else if (param == "iter_per_col_max")
         iter_per_col_max = std::get<double>(value);
   }

   void search(const MIP& mip, const std::vector<double>& lb,
               const std::vector<double>& ub,
               const std::vector<Activity>& activities,
               const LPResult& result, const std::vector<double>&,
               const std::vector<int>& fractional,
               std::shared_ptr<const LPSolver> solver, TimeLimit tlimit,
               SolutionPool& pool) override
   {
      auto heur_name = SELECTION::name + "Diving";
      int ncols = mip.getNCols();
      auto st = mip.getStats();
      const auto& objective = mip.getObj();
      const auto& downLocks = mip.getDownLocks();

      // avoid warnings
#ifndef NDEBUG
      const auto& upLocks = mip.getUpLocks();
#endif

      auto local_activities = activities;
      auto locallb = lb;
      auto localub = ub;
      auto localfrac = fractional;
      auto localsol = result.primalSolution;
      double localobj = -1.0;
      std::unique_ptr localsolver = solver->clone();

      std::vector<int> changedCols;
      std::vector<double> oldlbs;
      std::vector<double> oldubs;
      std::vector<Activity> old_activities;

      bool hasZeroLockFractionals = false;
      bool limit_reached = false;
      bool feasible = true;
      int iter = 0;
      do
      {
         ++iter;

         auto [varToFix, direction, nFrac] =
             SELECTION::select(mip, locallb, localub, localsol);

         Message::debug_details("{}: iter {}, nFrac {}", heur_name, iter,
                                nFrac);

         // no more variables to fix
         if (varToFix < 0)
         {
            if (nFrac > 0)
            {
               hasZeroLockFractionals = true;
#ifndef NDEBUG
               bool checklpFeas =
                   checkFeasibility<double, true>(mip, localsol);
               assert(checklpFeas);
#endif
            }
            else
               assert(checkFeasibility<double>(mip, localsol));

            break;
         }

         // remember the bounds and activities for backtracking
         if (propagate && backtrack)
         {
            oldlbs = locallb;
            oldubs = localub;
            old_activities = local_activities;
         }

         double fixedVarOldlb = locallb[varToFix];
         double fixedVarOldub = localub[varToFix];

         if (Num::isMinusInf(locallb[varToFix]) ||
             Num::isInf(localub[varToFix]))
         {
            locallb[varToFix] = Num::floor(localsol[varToFix]);
            localub[varToFix] = Num::ceil(localsol[varToFix]);
         }

         // fix the variable
         if (direction == 1)
            locallb[varToFix] = localub[varToFix];
         else
         {
            assert(direction == -1);
            localub[varToFix] = locallb[varToFix];
         }

         // propagate the change
         if (propagate)
         {
            bool status = propagate_get_changed_cols(
                mip, locallb, localub, local_activities, varToFix,
                fixedVarOldlb, fixedVarOldub, changedCols);

            Message::debug_details(
                "{}: fixed var {} -> {}, propagation: {} changed {}",
                heur_name, varToFix, locallb[varToFix], status,
                changedCols.size());

            if (!status && backtrack)
            {
               Message::debug_details("{}: backtraking", heur_name);

               // backtrack
               for (int col : changedCols)
               {
                  locallb[col] = oldlbs[col];
                  localub[col] = oldubs[col];
               }
               changedCols.clear();

               assert(locallb[varToFix] == oldlbs[varToFix]);
               assert(localub[varToFix] == oldubs[varToFix]);

               // correct the activities
               local_activities = old_activities;

               // try the opposite direction
               if (direction == 1)
                  localub[varToFix] = locallb[varToFix];
               else
                  locallb[varToFix] = localub[varToFix];

               status = propagate_get_changed_cols(
                   mip, locallb, localub, local_activities, varToFix,
                   fixedVarOldlb, fixedVarOldub, changedCols);

               if (!status)
               {
                  feasible = false;
                  Message::debug_details(
                      "{}: infeasible after backtracking + propagation",
                      heur_name);
                  break;
               }
            }
            else if (!status)
            {
               assert(0);
               feasible = false;
               break;
            }
         }

         assert(locallb[varToFix] == localub[varToFix]);

         // apply propagation changes and solve
         for (auto col : changedCols)
            localsolver->changeBounds(col, locallb[col], localub[col]);
         changedCols.clear();

         auto local_result = localsolver->solve(Algorithm::DUAL);

         if (local_result.status != LPResult::OPTIMAL)
         {
            // TODO backtrack ?
            feasible = false;
            Message::debug_details("{}: LP infeasible", heur_name);
         }
         else
         {
            localsol = std::move(local_result.primalSolution);
            localobj = local_result.obj;

            roundFeasIntegers(localsol, st.nbin + st.nint);
#ifndef NDEBUG
            bool checklpFeas =
                checkFeasibility<double, true>(mip, localsol);
            assert(checklpFeas);
            assert(Num::isFeasEQ(localsol[varToFix], locallb[varToFix]));
#endif
         }

         if (iter > ncols * iter_per_col_max ||
             tlimit.reached(Timer::now()))
            limit_reached = true;

      } while (feasible && !limit_reached);

      // add the solution if its feasible
      if (feasible && !limit_reached && hasZeroLockFractionals)
      {
         localobj = 0.0;

         for (int col = 0; col < ncols; ++col)
         {
            localobj += objective[col] * localsol[col];

            if (col < st.nbin + st.nint && !Num::isIntegral(localsol[col]))
            {
               assert(!downLocks[col] || !upLocks[col]);

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

         Message::debug("{}: found solution val {}", heur_name, localobj);
         assert(checkFeasibility<double>(mip, localsol));
         pool.add(std::move(localsol), localobj);
      }
      else if (feasible && !limit_reached)
      {
         Message::debug("{}: found solution val {}", heur_name, localobj);
         assert(checkFeasibility<double>(mip, localsol));
         pool.add(std::move(localsol), localobj);
      }
      else if (limit_reached)
         Message::debug("{}: limit reached after {} iterations", heur_name,
                        iter);
      else
         Message::debug("{}: infeasible after {} iterations", heur_name,
                        iter);
   }

 private:
   double iter_per_col_max = 0.3;
   bool propagate = true;
   bool backtrack = true;
};

#endif
