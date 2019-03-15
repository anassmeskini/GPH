#ifndef CPLEX_SOLVER_HPP
#define CPLEX_SOLVER_HPP

#ifdef CONCERT_CPLEX_FOUND

#include <ilcplex/ilocplex.h>
#include "core/LPSolver.h"
#include "core/MIP.h"

class CPXSolver : public LPSolver<double> {
  public:
   CPXSolver(const MIP<double>&, bool noOutput = true);

   // CPXSolver(MIP<double>&&);

   ~CPXSolver() override;

   virtual LPResult solve() override;

   using IntParam = std::pair<IloCplex::IntParam, int>;

   using BoolParam = std::pair<IloCplex::BoolParam, bool>;

   // void setIntParams(const std::initializer_list<IntParam>&);

   // void setBoolParams(const std::initializer_list<BoolParam>&);

  private:
   IloEnv env;
   IloCplex cplex;
   IloNumVarArray variables;
   IloRangeArray constraints;
};

#endif

#endif
