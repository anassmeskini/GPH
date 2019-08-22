#ifndef SOPLEX_HPP
#define SOPLEX_HPP

#ifdef SOPLEX_FOUND
#include "core/LPSolver.h"
#include "soplex.h"

class SPXSolver : public LPSolver
{
 public:
   SPXSolver(const MIP&);

   SPXSolver(const SPXSolver&);

   ~SPXSolver() override = default;

   LPResult solve(Algorithm) override;

   std::unique_ptr<LPSolver> makeCopy() const override;

   void changeBounds(int column, double lb, double ub) override;

   void changeBounds(const std::vector<double>&,
                     const std::vector<double>&) override;

   void changeObjective(int, double) override;

 private:
   soplex::SoPlex mysoplex;
   int ncols;
   int nrows;
};

#endif

#endif
