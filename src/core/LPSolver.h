#ifndef _LPSOLVER_HPP
#define _LPSOLVER_HPP

#include "MIP.h"
#include "fmt/format.h"
#include <memory>
#include <vector>

struct LPResult
{
   enum Status
   {
      INFEASIBLE,
      UNBOUNDED,
      OPTIMAL,
      ILL_CONDITIONED,
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

template<typename REAL>
class LPSolver
{
   public:
   virtual ~LPSolver() = default;

   virtual LPResult solve() = 0;

   virtual LPResult solve(LPAlgorithm);

   virtual std::unique_ptr<LPSolver<REAL>> clone() const = 0;
};

template<typename REAL>
LPResult LPSolver<REAL>::solve(LPAlgorithm)
{
   fmt::print("Info: specified LP algorithm was ignored");

   return solve();
}

#endif
