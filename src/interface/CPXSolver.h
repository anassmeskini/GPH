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

   LPResult solve(Algorithm) override;

   std::unique_ptr<LPSolver> makeCopy() const override;

   void changeBounds(int column, double lb, double ub) override;

   void changeBounds(const std::vector<double>&,
                     const std::vector<double>&) override;

   void changeObjective(int, double) override;

 private:
   IloEnv env;
   IloModel model;
   IloObjective objective;
   IloNumVarArray variables;
   IloRangeArray constraints;
   IloCplex cplex;

   // mapping of the order in which the IloRange are obtained
   // when iterating with IloIterator<IloRange> to the row number in the
   // problem
   // used to copy a model and keep the correct order in the IloRangeArray
   HashMap<int, int> exidTorowid;

   int ncols;
   int nrows;
};

#endif

#endif
