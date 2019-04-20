#include "LPFactory.h"
#include "AvaiLPSolver.h"

#include <mutex>

LPFactory::LPFactory(const MIP<double>& mip)
  : original(new AvaiLPSolver(mip))
{
}

std::unique_ptr<LPSolver<double>>
LPFactory::get() const
{
   std::lock_guard<tbb::mutex> lock(copyLock);

   return original->clone();
}

std::shared_ptr<LPSolver<double>>
LPFactory::getOriginal()
{
   return { original };
}

void
LPFactory::solve()
{
   original->solve();
}
