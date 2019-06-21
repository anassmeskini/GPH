#include "LPFactory.h"
#include "MySolver.h"

#include <mutex>

LPFactory::LPFactory(const MIP<double>& mip)
  : original(new MySolver(mip))
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
