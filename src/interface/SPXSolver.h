#ifndef SOPLEX_HPP
#define SOPLEX_HPP

#ifdef SOPLEX_FOUND
#include "core/LPSolver.h"
#include "soplex.h"

class SPXSolver : public LPSolver<double> {
  public:
   SPXSolver(const MIP<double>&);

   virtual ~SPXSolver();

   virtual LPResult solve() override;

  private:
   soplex::SoPlex mysoplex;
};

#endif

#endif
