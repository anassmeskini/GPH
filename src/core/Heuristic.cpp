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

   runtime.resize(heuristics.size(), 0.0);
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

      SOLFormat::write("rounded_lp.sol", sol.first, mip.getVarNames());
      pool.add(std::move(sol.first), sol.second);
   }

   auto run = [&](tbb::blocked_range<size_t>& range) -> void {
      for (size_t i = range.begin(); i != range.end(); ++i)
         heuristics[i]->search(mip, mip.getLB(), mip.getUB(), activities,
                               result, lpSolAct, fractional, lpSolver,
                               pool);
   };

   tbb::parallel_for(tbb::blocked_range<size_t>{0, heuristics.size()},
                     std::move(run));
   auto tend = Timer::now();

   // TODO
   if (pool.pool.size())
   {
      double mincost = pool.pool[0].objective;
      for (size_t i = 1; i < pool.pool.size(); ++i)
         mincost = std::min(mincost, pool.pool[i].objective);

      double gap = 100.0 * std::fabs(mincost - result.obj) /
                   (std::fabs(mincost) + 1e-6);

      if (gap < 999.9)
         Message::print(
             "Found {} solutions | gap {:0.2f}% after {:0.2} sec.",
             pool.pool.size(), gap, Timer::seconds(tend, t0));
      else
         Message::print("Found {} solutions | gap --- after {:0.2} sec.",
                        pool.pool.size(), Timer::seconds(tend, t0));
   }
   else
      Message::print("No solution found");
}
