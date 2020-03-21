#include "TrivialRounding.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

#include <random>

void
TrivialRounding::search(const MIP& mip, const std::vector<double>&,
                        const std::vector<double>&,
                        const std::vector<Activity>&,
                        const LPResult& result, const std::vector<double>&,
                        const std::vector<int>& fractional,
                        std::shared_ptr<const LPSolver>, TimeLimit,
                        SolutionPool& pool)
{

   std::optional optSol =
       minLockRound(mip, result.primalSol, result.obj, fractional);

   if (optSol)
   {
      auto& sol = optSol.value();
      pool.add(std::move(sol.first), sol.second);
   }
}
