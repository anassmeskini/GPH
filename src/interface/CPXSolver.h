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

   virtual LPResult solve(LPAlgorithm) override;

   virtual std::unique_ptr<LPSolver<double>> clone() const override;

   private:
   IloEnv env;
   IloModel model;
   IloNumVarArray variables;
   IloRangeArray constraints;
   IloCplex cplex;
   bool deleteEnv;

   size_t ncols;
   size_t nrows;
};

#endif

#endif
