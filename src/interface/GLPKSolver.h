#ifndef GLPK_SOLVER_HPP
#define GLPK_SOLVER_HPP

#ifdef GLPK_FOUND

#include "core/LPSolver.h"
#include "core/MIP.h"
#include <glpk.h>

class GLPKSolver : public LPSolver
{
 public:
   GLPKSolver(const MIP&);

   GLPKSolver(const GLPKSolver&);

   ~GLPKSolver() override;

   LPResult solve() override;

   std::unique_ptr<LPSolver> makeCopy() const override;

   void changeBounds(int column, double lb, double ub) override;
   :q!

:q!
:visual
visual

   void changeBounds(const std::vector<double>&,
                     const std::vector<double>&) override;

 private:
   glp_prob* problem;

   int ncols;
   int nrows;
};

#endif

#endif
