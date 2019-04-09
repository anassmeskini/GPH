#ifndef SOPLEX_HPP
#define SOPLEX_HPP

#ifdef SOPLEX_FOUND
#include "core/LPSolver.h"
#include "soplex.h"

class SPXSolver : public LPSolver<double>
{
public:
  SPXSolver(const MIP<double>&);

  SPXSolver(const SPXSolver&);

  virtual ~SPXSolver() override;

  virtual LPResult solve() override;

  virtual std::unique_ptr<LPSolver<double>> clone() const override;

private:
  soplex::SoPlex mysoplex;
  size_t ncols;
  size_t nrows;
};

#endif

#endif
