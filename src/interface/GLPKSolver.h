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

   virtual ~GLPKSolver() override;

   virtual LPResult solve() override;

   virtual std::unique_ptr<LPSolver> makeCopy() const override;

   virtual void branch(int column, double val, Direction direction) override;

   virtual void changeBounds(int column, double lb, double ub) override;

   virtual void changeBounds(const std::vector<double>&,
                             const std::vector<double>&) override;

   private:
   glp_prob* problem;

   int ncols;
   int nrows;
};

#endif

#endif
