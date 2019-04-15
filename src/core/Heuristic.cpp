#include "Heuristic.h"
#include "AvaiLPSolver.h"
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

std::optional<std::vector<double>>
Heuristics::run(MIP<double>&& mip)
{
   LPFactory lpFactory(mip);

   auto solver = lpFactory.getOriginal();
   auto result = solver->solve();

   ProblemView problem(std::move(mip), std::move(result.primalSolution));

   Message::print("LP solver return status: {}\n", to_str(result.status));

   // TODO throw exception ??
   if (result.status != LPResult::OPTIMAL)
      return {};

   else if (result.status == LPResult::OPTIMAL)
      Message::print("obj: {}\n", result.obj);

   auto run = [this, &problem](tbb::blocked_range<size_t>& range) {
      for (size_t i = range.begin(); i != range.end(); ++i)
         heuristics[i]->search(problem);
   };

   tbb::parallel_for(tbb::blocked_range<size_t>{ 0, heuristics.size() },
                     std::move(run));

   return {};
}
