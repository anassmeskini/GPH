#include "CPXSolver.h"
#include <iostream>

#ifdef CONCERT_CPLEX_FOUND

CPXSolver::CPXSolver(const MIP<double>& linearProgram, bool noOutput)
    : mip(linearProgram),
      LPSolver(linearProgram),
      variables(env),
      constraints(env) {
   IloModel model(env);
   IloNumExpr objExpr(env);

   const auto& lb = linearProgram.getLB();
   const auto& ub = linearProgram.getUB();
   const auto& obj = linearProgram.getObj();
   const auto& lhs = linearProgram.getLHS();
   const auto& rhs = linearProgram.getRHS();

   // add variables to the model and build the objective expression
   for (size_t var = 0; var < linearProgram.getNCols(); ++var) {
      variables.add(IloNumVar(env, lb[var], ub[var]));
      objExpr += variables[var] * obj[var];
   }

   model.add(IloMinimize(env, objExpr));

   for (size_t row = 0; row < linearProgram.getNRows(); ++row) {
      auto rowView = linearProgram.getRow(row);
      const double* coefs = rowView.coefs;
      const size_t* indices = rowView.indices;

      IloNumExpr rowExpr(env);
      for (size_t id = 0; id < rowView.size; ++id)
         rowExpr += rowView.coefs[id] * variables[indices[id]];

      constraints.add(lhs[row] <= rowExpr <= rhs[row]);
   }

   model.add(constraints);
   cplex = model;

   if (noOutput) cplex.setOut(env.getNullStream());
}

LPResult CPXSolver::solve() {
   cplex.solve();

   LPResult result;
   auto status = cplex.getStatus();

   if (status == IloAlgorithm::Optimal) {
      result.status = LPResult::OPTIMAL;

      IloNumArray vals(env);
      // get primal vals
      size_t ncols = mip.getNCols();
      cplex.getValues(vals, variables);
      for (size_t i = 0; i < ncols; ++i)
         result.primalSolution.push_back(vals[i]);

      // get dual vals
      size_t nrows = mip.getNRows();
      cplex.getDuals(vals, constraints);
      for (size_t i = 0; i < ncols; ++i) result.dualSolution.push_back(vals[i]);

      result.obj = cplex.getObjValue();
   } else if (IloAlgorithm::Unbounded)
      result.status = LPResult::UNBOUNDED;
   else if (IloAlgorithm::Infeasible)
      result.status = LPResult::INFEASIBLE;
   else
      result.status = LPResult::OTHER;

   return result;
}

/*
void CPXSolver::setIntParams(const std::initializer_list<IntParam>& params) {
   for (auto param : params) cplex.setParam(param.first, param.second);
}

void CPXSolver::setBoolParams(const std::initializer_list<BoolParam>& params) {
   for (auto param : params) cplex.setParam(param.first, param.second);
}*/

CPXSolver::~CPXSolver() { env.end(); }

#endif
