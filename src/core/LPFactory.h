#ifndef LPFACTORY_HPP
#define LPFACTORY_HPP

#include "LPSolver.h"
#include "MIP.h"
#include "tbb/mutex.h"

class LPFactory
{
   public:
   LPFactory(const MIP<double>&);

   LPFactory(const LPFactory&);

   LPFactory(LPFactory&&);

   std::unique_ptr<LPSolver<double>> get() const;

   LPResult solve();

   void branch(size_t col, double val, Direction direction);

   private:
   std::shared_ptr<LPSolver<double>> solver;
   mutable tbb::mutex copyLock;
};

#endif
