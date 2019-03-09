#ifndef _DUMMY_SOLVER_HPP
#define _DUMMY_SOLVER_HPP

#include "LPSolver.h"
#include "MIP.h"

class DSolver : public LPSolver<double> {
  public:
   DSolver(const MIP<double>& mip) : LPSolver(mip) {}

   LPResult solve() { return {}; }
};

#endif
