#include "MinFracRounding.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

void
MinFracRounding::search(const MIP& mip, const std::vector<double>& lb,
                        const std::vector<double>& ub,
                        const std::vector<Activity>& activities,
                        const LPResult& result, const std::vector<double>&,
                        const std::vector<int>& fractional,
                        std::shared_ptr<const LPSolver> solver,
                        SolutionPool& pool)
{
   int ncols = mip.getNCols();

   auto st = mip.getStats();
   const auto& objective = mip.getObj();

   std::unique_ptr<LPSolver> localsolver;

   auto solution = result.primalSolution;

   auto process_sol = [&](std::vector<double>& sol) {
      if (mip.getStats().ncont == 0)
      {
         Message::debug("FracRound: feasible, no continuous variables");
         double cost = 0.0;

         for (int col = 0; col < ncols; ++col)
            cost += sol[col] * objective[col];
         pool.add(std::move(sol), cost);
      }
      else
      {
         Message::debug("FracRound: feasible, solving lp");
         if (!localsolver)
            localsolver = solver->clone();

         for (int col = 0; col < st.nbin + st.nint; ++col)
         {
            localsolver->changeBounds(col, sol[col], sol[col]);
         }

         auto res = localsolver->solve(Algorithm::DUAL);

         // ??
         if (res.status == LPResult::OPTIMAL)
            pool.add(std::move(res.primalSolution), res.obj);
         else
            Message::debug("FracRound: lp infeasible");
      }
   };

   for (size_t i = 0; i < fractional.size(); ++i)
   {
      int col = fractional[i];

      double floor = Num::floor(solution[col]);
      double frac = solution[col] - floor;
      assert(col < st.nbin + st.nint && frac > 0.0);

      solution[col] =
          floor + (frac > 0.5) + (frac == 0.5 && objective[col] < 0.0);
   }

   if (getNViolated<double>(mip, solution) == 0)
   {
      process_sol(solution);
      return;
   }

   Message::debug("FracRound: infeasible, trying propagation");

   auto locallb = lb;
   auto localub = ub;
   auto local_activities = activities;
   bool feasible = true;

   for (int col = 0; col < st.nbin + st.nint && feasible; ++col)
   {
      if (!Num::isFeasEQ(locallb[col], localub[col]))
      {
         double oldlb = locallb[col];
         double oldub = localub[col];
         double floor = Num::floor(result.primalSolution[col]);
         double frac = solution[col] - floor;

         locallb[col] =
             floor + (frac > 0.5) + (frac == 0.5 && objective[col] < 0.0);
         localub[col] = locallb[col];

         if (!propagate(mip, locallb, localub, local_activities, col,
                        oldlb, oldub))
            feasible = false;
      }
   }

   if (feasible)
      process_sol(locallb);
   else
      Message::debug("FracRound: no solution found");
}
