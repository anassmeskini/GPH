#include "Heuristic.h"
#include "MySolver.h"
#include "Numerics.h"
#include "Timer.h"
#include "io/Message.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

#include <cassert>
#include <iostream>

Search::Search(std::initializer_list<HeuristicMethod*> list)
{
   for (auto handle : list)
      heuristics.emplace_back(handle);
}

void
Search::run(const MIP& mip)
{
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

   Message::print("solving lp took {:0.2f} -> {} | obj: {}",
                  Timer::seconds(t1, t0),
                  to_str(result.status),
                  result.obj);

   // TODO should compute both
   auto lpSolAct = computeSolActivities(mip, result.primalSolution);
   auto fractional = getFractional(result.primalSolution, mip.getInteger());
   auto activities = computeActivities(mip);

   tbb::parallel_for(tbb::blocked_range<size_t>{ 0, heuristics.size() },
                     [&](tbb::blocked_range<size_t>& range) {
                        for (size_t i = range.begin(); i != range.end(); ++i)
                           heuristics[i]->search(mip,
                                                 mip.getLB(),
                                                 mip.getUB(),
                                                 activities,
                                                 result,
                                                 lpSolAct,
                                                 fractional,
                                                 lpSolver);
                     });
}

Search::~Search() {}
