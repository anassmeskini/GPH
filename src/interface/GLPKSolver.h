#ifndef GLPK_SOLVER_HPP
#define GLPK_SOLVER_HPP

#ifdef GLPK_FOUND

#include <glpk.h>
#include "core/LPSolver.h"
#include "core/MIP.h"

class GLPKSolver : public LPSolver<double> {
  public:
   GLPKSolver(const MIP<double>&, bool noOutput = true);

   ~GLPKSolver() override;

   virtual LPResult solve() override;

  private:
   glp_prob* problem;
};

#endif

#endif
