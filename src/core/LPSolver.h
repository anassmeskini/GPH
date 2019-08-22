#ifndef LPSOLVER_HPP
#define LPSOLVER_HPP

#include "MIP.h"
#include "fmt/format.h"
#include <memory>
#include <mutex>
#include <tbb/mutex.h>
#include <vector>

struct LPResult
{
   enum Status
   {
      INFEASIBLE,
      UNBOUNDED,
      OPTIMAL,
      OTHER
   } status;

   std::vector<double> primalSolution;
   std::vector<double> dualSolution;
   double obj;
};

std::string to_str(LPResult::Status);

enum class Algorithm
{
   PRIMAL,
   DUAL,
};

class LPSolver
{
 public:
   LPSolver() = default;

   virtual ~LPSolver() = default;

   virtual LPResult solve(Algorithm) = 0;

   std::unique_ptr<LPSolver> clone() const
   {
      std::unique_lock guard(copyLock);
      return makeCopy();
   }

   virtual void changeBounds(int column, double lb, double ub) = 0;

   virtual void changeBounds(const std::vector<double>&,
                             const std::vector<double>&) = 0;

   virtual void changeObjective(int column, double coef) = 0;

 private:
   virtual std::unique_ptr<LPSolver> makeCopy() const = 0;

   mutable tbb::mutex copyLock;
};

#endif
