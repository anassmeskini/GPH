#include "LPFactory.h"
#include "AvaiLPSolver.h"

LPFactory::LPFactory(const MIP<double>& mip)
  : original(new AvaiLPSolver(mip))
{}

std::unique_ptr<LPSolver<double>>
LPFactory::get() const
{
   copyLock.lock();

   auto ptr = original->clone();

   copyLock.unlock();

   return ptr;
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
