#ifndef GLPK_SOLVER_HPP
#define GLPK_SOLVER_HPP

#ifdef GLPK_FOUND

#include "core/LPSolver.h"
#include "core/MIP.h"
#include <glpk.h>

class GLPKSolver : public LPSolver<double>
{
   public:
   GLPKSolver(const MIP<double>&);

   GLPKSolver(const GLPKSolver&);

   virtual ~GLPKSolver() override;

   virtual LPResult solve() override;

   virtual LPResult solve(LPAlgorithm) override;

   virtual std::unique_ptr<LPSolver<double>> makeCopy() const override;

   virtual void branch(int column, double val, Direction direction) override;

   private:
   glp_prob* problem;

   size_t ncols;
   size_t nrows;
};

#endif

#endif
