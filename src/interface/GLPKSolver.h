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

   virtual std::unique_ptr<LPSolver<double>> clone() const override;

 private:
   glp_prob* problem;

   size_t ncols;
   size_t nrows;
};

#endif

#endif
