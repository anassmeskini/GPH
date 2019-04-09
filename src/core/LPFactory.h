#ifndef LPFACTORY_HPP
#define LPFACTORY_HPP

#include "LPSolver.h"
#include "MIP.h"
#include "tbb/mutex.h"

class LPFactory
{
 public:
   LPFactory(const MIP<double>&);

   std::unique_ptr<LPSolver<double>> get() const;

   std::shared_ptr<LPSolver<double>> getOriginal();

   void solve();

 private:
   std::shared_ptr<LPSolver<double>> original;
   mutable tbb::mutex copyLock;
};

#endif
