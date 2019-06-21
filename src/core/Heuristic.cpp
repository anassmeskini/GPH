#include "Heuristic.h"
#include "LPFactory.h"
#include "Numerics.h"
#include "io/Message.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

#include <cassert>
#include <iostream>

Heuristics::Heuristics(const std::initializer_list<HeuristicMethod*>& list)
{
   for (auto* handle : list)
      heuristics.emplace_back(handle);
}

void
Heuristics::run(MIP<double>&& mip)
{
   LPFactory lpFactory(mip);

   auto solver = lpFactory.getOriginal();
   auto result = solver->solve();

   Message::print("LP solver return status: {}\n", to_str(result.status));

   // TODO throw exception ??
   if (result.status != LPResult::OPTIMAL)
      assert(0);

   else if (result.status == LPResult::OPTIMAL)
      Message::print("obj: {}\n", result.obj);

   ProblemView problem(std::move(mip), std::move(result.primalSolution));

   auto run = [this, &problem](tbb::blocked_range<size_t>& range) {
      for (size_t i = range.begin(); i != range.end(); ++i)
         heuristics[i]->search(problem);
   };

   tbb::parallel_for(tbb::blocked_range<size_t>{ 0, heuristics.size() },
                     std::move(run));
}
