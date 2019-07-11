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

   virtual ~SPXSolver() override;

   virtual LPResult solve() override;

   virtual std::unique_ptr<LPSolver> makeCopy() const override;

   virtual void branch(int column, double val, Direction direction) override;

   private:
   soplex::SoPlex mysoplex;
   int ncols;
   int nrows;
};

#endif

#endif
