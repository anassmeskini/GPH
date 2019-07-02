#include "LPFactory.h"
#include "MySolver.h"

#include <mutex>

LPFactory::LPFactory(const MIP<double>& mip)
  : solver(new MySolver(mip))
{
}

LPFactory::LPFactory(const LPFactory& other) : solver(std::move(other.solver->clone()))
{
}

LPFactory::LPFactory(LPFactory&& source) : solver(source.solver)
{
}

std::unique_ptr<LPSolver<double>>
LPFactory::get() const
{
   std::lock_guard<tbb::mutex> lock(copyLock);

   return solver->clone();
}

LPResult
LPFactory::solve()
{
   std::lock_guard<tbb::mutex> lock(copyLock);

   return solver->solve();
}

void
LPFactory::branch(size_t col, double val, Direction direction)
{
   solver->branch(col, val, direction);
}


