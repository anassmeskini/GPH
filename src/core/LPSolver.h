#ifndef _LPSOLVER_HPP
#define _LPSOLVER_HPP

#include "MIP.h"
#include <memory>
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

template<typename REAL>
class LPSolver
{
public:
  virtual ~LPSolver(){};

  virtual LPResult solve() = 0;

  virtual std::unique_ptr<LPSolver<REAL>> clone() const = 0;
};

#endif
