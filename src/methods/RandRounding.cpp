#include "RandRounding.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

#include <random>

void
RandRounding::search(const MIP& mip, const std::vector<double>& lb,
                     const std::vector<double>& ub,
                     const std::vector<Activity>& activities,
                     const LPResult& result, const std::vector<double>&,
                     const std::vector<int>& fractional,
                     std::shared_ptr<const LPSolver> solver,
                     TimeLimit tlimit, SolutionPool& pool)
{
   auto st = mip.getStats();
   auto objective = mip.getObj();

   std::unique_ptr<LPSolver> localsolver;

   static std::default_random_engine gen;
   std::uniform_real_distribution<double> dist(0.0, 1.0);

   auto process_sol = [&](std::vector<double>& sol) {
      if (mip.getStats().ncont == 0)
      {
         double cost = 0.0;

         for (int col = 0; col < st.ncols; ++col)
            cost += sol[col] * objective[col];
         pool.add(std::move(sol), cost);

         Message::debug(
             "RandRound: feasible, no continuous variables, cost {}",
             cost);
      }
      else
      {
         Message::debug("RandRound: feasible, solving lp");
         if (localsolver)
            localsolver.release();

         localsolver = solver->clone();

         for (int col = 0; col < st.nbin + st.nint; ++col)
            localsolver->changeBounds(col, sol[col], sol[col]);

         auto res = localsolver->solve(Algorithm::DUAL);

         // ??
         if (res.status == LPResult::OPTIMAL)
         {
            pool.add(std::move(res.primalSolution), res.obj);
            Message::debug("RandRound: lp feasible, cost {}", res.obj);
         }
         else
            Message::debug("RandRound: lp infeasible");
      }
   };

   std::vector<double> locallb_partial = lb;
   std::vector<double> localub_partial = ub;
   std::vector<Activity> local_activities_partial = activities;

   bool feasible = true;
   for (int col = 0; col < st.nbin + st.nint && feasible; ++col)
   {
      double solval = result.primalSolution[col];
      if (Num::isIntegral(solval))
      {
         double oldlb = locallb_partial[col];
         double oldub = localub_partial[col];

         assert(Num::isFeasGE(solval, lb[col]));
         assert(Num::isFeasLE(solval, ub[col]));

         locallb_partial[col] = solval;
         localub_partial[col] = solval;

         feasible = propagate(mip, locallb_partial, localub_partial,
                              local_activities_partial, col, oldlb, oldub);
      }

      if (tlimit.reached(Timer::now()))
         return;
   }

   if (!feasible)
   {
      Message::debug_details(
          "RandRound: infeasible after propagating integral lp values");
      return;
   }

   std::vector<double> locallb(st.ncols);
   std::vector<double> localub(st.ncols);
   std::vector<Activity> local_activities(st.nrows);

   int iter = 0;
   do
   {
      ++iter;
      Message::debug("RandRound: iter {}", iter);

      locallb = locallb_partial;
      localub = localub_partial;
      local_activities = local_activities_partial;

      for (size_t i = 0; i < fractional.size() && feasible; ++i)
      {
         int col = fractional[i];
         if (!Num::isFeasEQ(locallb[col], localub[col]))
         {
            double rand = dist(gen);
            double floor = Num::floor(result.primalSolution[col]);
            double frac = result.primalSolution[col] - floor;

            double oldlb = locallb[col];
            double oldub = localub[col];

            locallb[col] = floor + static_cast<double>(rand <= frac);
            localub[col] = locallb[col];

            feasible = propagate(mip, locallb, localub, local_activities,
                                 col, oldlb, oldub);

            if (tlimit.reached(Timer::now()))
               return;
         }
      }

      if (feasible)
         process_sol(locallb);
      else
         Message::debug("RandRound: No solution found");
   } while (feasible && iter < 5);
}
