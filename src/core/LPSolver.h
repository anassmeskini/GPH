#ifndef _LPSOLVER_HPP
#define _LPSOLVER_HPP

#include "MIP.h"
#include <vector>

// TODO templetize this class
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

template<typename REAL>
class LPSolver
{
 public:
   LPSolver(const MIP<REAL>&);

   virtual ~LPSolver(){};

   virtual LPResult solve() = 0;

 protected:
   const MIP<REAL>& mip;
};

template<typename REAL>
LPSolver<REAL>::LPSolver(const MIP<REAL>& _mip)
  : mip(_mip)
{}

#endif
