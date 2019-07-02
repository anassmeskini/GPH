#include "Heuristic.h"
#include "LPFactory.h"
#include "Numerics.h"
#include "Timer.h"
#include "io/Message.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

#include <cassert>
#include <iostream>

Heuristics::Heuristics(std::initializer_list<HeuristicMethod*> list)
{
   for (auto handle : list)
      heuristics.emplace_back(handle);
}

void
Heuristics::run(MIP<double>&& mip)
{
   // int ncols = mip.getNCols();
   LPFactory lpFactory(mip);

   // variables to be captured by the lambda
   LPResult result;

   auto t0 = Timer::now();
   result = lpFactory.solve();
   auto t1 = Timer::now();

   Message::print("solving lp took {:0.2f} -> {} | obj: {}",
                  Timer::seconds(t1, t0),
                  to_str(result.status),
                  result.obj);

   std::vector<double> lpSolAct =
     computeSolActivities(mip, result.primalSolution);
   std::vector<size_t> fractional;

   if (result.status != LPResult::OPTIMAL &&
       result.status != LPResult::INFEASIBLE)
   {
      Message::print("unkown status");
      assert(0);
   }

   tbb::parallel_for(tbb::blocked_range<size_t>{ 0, heuristics.size() },
                     [&](tbb::blocked_range<size_t>& range) {
                        for (size_t i = range.begin(); i != range.end(); ++i)
                           heuristics[i]->search(mip,
                                                 mip.getLB(),
                                                 mip.getUB(),
                                                 result,
                                                 lpSolAct,
                                                 fractional,
                                                 lpFactory);
                     });
}

Heuristics::~Heuristics()
{
   // TODO
}
