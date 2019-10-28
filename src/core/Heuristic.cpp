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

std::tuple<int, int, double, int>
Search::getSolSummary() const
{
   int min_cost_heur = -1;
   int min_cost_sol = -1;
   double min_cost = std::numeric_limits<double>::max();
   int nsols = 0;

   for (size_t i = 0; i < heuristics_solutions.size(); ++i)
   {
      nsols += heuristics_solutions[i].size();

      for (size_t j = 0; j < heuristics_solutions[i].size(); ++j)
      {
         if (heuristics_solutions[i][j].second < min_cost)
         {
            min_cost = heuristics_solutions[i][j].second;
            min_cost_heur = i;
            min_cost_sol = j;
         }
      }
   }

   return {min_cost_heur, min_cost_sol, min_cost, nsols};
}

bool
Search::checkSolFeas(const MIP& mip) const
{
   for (size_t i = 0; i < heuristics_solutions.size(); ++i)
   {
      for (size_t j = 0; j < heuristics_solutions[i].size(); ++j)
      {
         if (!checkFeasibility(mip, heuristics_solutions[i][j].first, 1e-6,
                               1e-9))
         {
            Message::debug("{} solution number {} is INFEASIBLE",
                           heuristics[i]->getName(), j);
            return false;
         }
      }
   }
   return true;
}

std::optional<std::vector<double>>
Search::run(const MIP& mip, int seconds)
{
   TimeLimit tlimit(Timer::now(), seconds);
   auto st = mip.getStats();
   auto lpSolver = std::make_shared<MySolver>(mip);

   // variables to be captured by the lambda
   LPResult result;

   auto t0 = Timer::now();
   result = lpSolver->solve(Algorithm::DUAL);
   auto t1 = Timer::now();

   // TODO
   if (result.status != LPResult::OPTIMAL &&
       result.status != LPResult::INFEASIBLE)
      assert(0);

   auto lpFeas = checkFeasibility<double, true>;
   assert(lpFeas(mip, result.primalSolution, 1e-9, 1e-6));

   roundFeasIntegers(result.primalSolution, st.nbin + st.nint);

   //
   // TODO should compute both
   auto lpSolAct = computeSolActivities(mip, result.primalSolution);
   auto fractional =
       getFractional(result.primalSolution, st.nbin + st.nint);
   auto activities = computeActivities(mip);

   double percfrac = 100.0 * static_cast<double>(fractional.size()) /
                     (st.nbin + st.nint);
   Message::print("solving LP took {:0.2f} sec. npivots {} -> {} | obj: "
                  "{:0.4e} | fractional: "
                  "{} ({:0.1f}%)",
                  Timer::seconds(t1, t0), result.niter,
                  to_str(result.status), result.obj, fractional.size(),
                  percfrac);

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
                                tlimit, heuristics_solutions[i]);
   };

   tbb::parallel_for(tbb::blocked_range<size_t>{0, heuristics.size()},
                     std::move(run));
   auto tend = Timer::now();

   assert(checkSolFeas(mip));

   auto [min_cost_heur, min_cost_sol, min_cost, nsols] = getSolSummary();

   if (nsols > 0)
   {
      assert(min_cost_sol != -1 && min_cost_heur != -1);

      double gap = 100.0 * std::fabs(min_cost - result.obj) /
                   (std::fabs(result.obj) + 1e-6);

      if (gap < 10000.0)
         Message::print(
             "Found {} solutions | gap {:0.2f}% after {:0.2} sec.", nsols,
             gap, Timer::seconds(tend, t0));
      else
         Message::print("Found {} solutions | gap --- after {:0.2} sec.",
                        nsols, Timer::seconds(tend, t0));
   }
   else
      Message::print("No solution found after {} sec.",
                     Timer::seconds(tend, t0));

   Message::print("Summary:");
   Message::print("\tRoot: found {} solutions",
                  heuristics_solutions.back().size());
   for (size_t i = 0; i < heuristics.size(); ++i)
   {
      std::string best;
      if (i == static_cast<size_t>(min_cost_heur))
         best = "\t[best]";
      Message::print("\t{}: {:0.1f} sec. found {} solutions {}",
                     heuristics[i]->getName(), heuristics[i]->getRunTime(),
                     heuristics_solutions[i].size(), best);
   }

   if (nsols > 0)
      return heuristics_solutions[min_cost_heur][min_cost_sol].first;

   return {};
}
