#ifndef CPLEX_SOLVER_HPP
#define CPLEX_SOLVER_HPP

#ifdef CONCERT_CPLEX_FOUND

#include "core/LPSolver.h"
#include "core/MIP.h"
#include <ilcplex/ilocplex.h>

class CPXSolver : public LPSolver
{
   public:
   CPXSolver(const MIP&);

   CPXSolver(const CPXSolver&);

   ~CPXSolver() override;

   virtual LPResult solve() override;

   virtual std::unique_ptr<LPSolver> makeCopy() const override;

   virtual void branch(int column, double val, Direction direction) override;

   virtual void changeBounds(int column, double lb, double ub) override;

   virtual void changeBounds(const std::vector<double>&,
                             const std::vector<double>&) override;

   private:
   IloEnv env;
   IloModel model;
   IloNumVarArray variables;
   IloRangeArray constraints;
   IloCplex cplex;
   bool deleteEnv;

   int ncols;
   int nrows;
};

#endif

#endif
