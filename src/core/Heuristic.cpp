#include "Heuristic.h"
#include "MySolver.h"
#include "Numerics.h"
#include "Timer.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

#include <cassert>
#include <mutex>

Search::Search(std::initializer_list<HeuristicMethod*> list)
{
   for (auto handle : list)
      heuristics.emplace_back(handle);

   heuristics_solutions.resize(heuristics.size() + 1);
}

void
Search::run(const MIP& mip)
{
   auto st = mip.getStatistics();
   auto lpSolver = std::make_shared<MySolver>(mip);

   // variables to be captured by the lambda
   LPResult result;

   auto t0 = Timer::now();
   result = lpSolver->solve();
   auto t1 = Timer::now();

   // TODO
   if (result.status != LPResult::OPTIMAL &&
       result.status != LPResult::INFEASIBLE)
      assert(0);
   //
   // TODO should compute both
   auto lpSolAct = computeSolActivities(mip, result.primalSolution);
   auto fractional =
       getFractional(result.primalSolution, mip.getInteger());
   auto activities = computeActivities(mip);

   double percfrac = 100.0 * static_cast<double>(fractional.size()) /
                     (st.nbin + st.nint);
   Message::print(
       "solving lp took {:0.2f} -> {} | obj: {:0.4e} | fractional: "
       "{} ({:0.1f}%)",
       Timer::seconds(t1, t0), to_str(result.status), result.obj,
       fractional.size(), percfrac);

   if (auto optSol = minLockRound(mip, result.primalSolution, result.obj,
                                  fractional))
   {
      auto& sol = optSol.value();
      Message::print("Root lp can be rounded, obj {}", sol.second);

      heuristics_solutions.back().add(std::move(sol.first), sol.second);
   }

   auto run = [&](tbb::blocked_range<size_t>& range) -> void {
      for (size_t i = range.begin(); i != range.end(); ++i)
         heuristics[i]->execute(mip, mip.getLB(), mip.getUB(), activities,
                                result, lpSolAct, fractional, lpSolver,
                                heuristics_solutions[i]);
   };

   tbb::parallel_for(tbb::blocked_range<size_t>{0, heuristics.size()},
                     std::move(run));
   auto tend = Timer::now();

   double min_cost = std::numeric_limits<double>::max();
   std::pair<int, int> min_cost_sol{-1, -1};
   int nsol = 0;
   for (size_t i = 0; i < heuristics_solutions.size(); ++i)
   {
      nsol += heuristics_solutions[i].size();

      for (size_t j = 0; j < heuristics_solutions[i].size(); ++j)
      {
         if (heuristics_solutions[i][j].second < min_cost)
         {
            min_cost = heuristics_solutions[i][j].second;
            min_cost_sol = {i, j};
         }
      }
   }

   if (min_cost_sol != std::make_pair(-1, -1))
   {
      assert(min_cost_sol.first != -1 && min_cost_sol.second != -1);

      double gap = 100.0 * std::fabs(min_cost - result.obj) /
                   (std::fabs(min_cost) + 1e-6);

      if (gap < 999.9)
         Message::print(
             "Found {} solutions | gap {:0.2f}% after {:0.2} sec.", nsol,
             gap, Timer::seconds(tend, t0));
      else
         Message::print("Found {} solutions | gap --- after {:0.2} sec.",
                        nsol, Timer::seconds(tend, t0));
   }
   else
      Message::print("No solution found");

   Message::print("Summary:");
   Message::print("\tRoot: found {} solutions",
                  heuristics_solutions.back().size());
   for (size_t i = 0; i < heuristics.size(); ++i)
      Message::print("\t{}: {:0.1f} sec. found {} solutions",
                     heuristics[i]->getName(), heuristics[i]->getRunTime(),
                     heuristics_solutions[i].size());
}
