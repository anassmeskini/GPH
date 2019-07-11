#ifndef _LPSOLVER_HPP
#define _LPSOLVER_HPP

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

enum class LPAlgorithm
{
   PRIMAL_SPX,
   DUAL_SPX,
   BARRIER,
   AUTO
};

enum class Direction
{
   UP,
   Down
};

class LPSolver
{
   public:
   virtual ~LPSolver() = default;

   virtual LPResult solve() = 0;

   std::unique_ptr<LPSolver> clone() const
   {
      std::lock_guard guard(copyLock);
      return makeCopy();
   }

   virtual void branch(int column, double val, Direction direction) = 0;

   virtual std::unique_ptr<LPSolver> makeCopy() const = 0;

   private:
   mutable tbb::mutex copyLock;
};

#endif
