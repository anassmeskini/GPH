#ifndef CPLEX_SOLVER_HPP
#define CPLEX_SOLVER_HPP

#ifdef CONCERT_CPLEX_FOUND

#include "core/LPSolver.h"
#include "core/MIP.h"
#include <ilcplex/ilocplex.h>

class CPXSolver : public LPSolver<double>
{
public:
  CPXSolver(const MIP<double>&);

  CPXSolver(const CPXSolver&);

  ~CPXSolver() override;

  virtual LPResult solve() override;

  virtual std::unique_ptr<LPSolver<double>> clone() const override;

  // using IntParam = std::pair<IloCplex::IntParam, int>;

  // using BoolParam = std::pair<IloCplex::BoolParam, bool>;

private:
  IloEnv env;
  IloModel model;
  IloNumVarArray variables;
  IloRangeArray constraints;
  IloCplex cplex;

  size_t ncols;
  size_t nrows;
  bool deleteEnv;
};

#endif

#endif
