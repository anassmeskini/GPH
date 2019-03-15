#include "interface/CPXSolver.h"
#include <iostream>

#ifdef CONCERT_CPLEX_FOUND

CPXSolver::CPXSolver(const MIP<double>& linearProgram, bool noOutput)
    : LPSolver(linearProgram), variables(env), constraints(env) {
   IloModel model(env);
   IloNumExpr objExpr(env);

   const auto& lb = linearProgram.getLB();
   const auto& ub = linearProgram.getUB();
   const auto& obj = linearProgram.getObj();
   const auto& lhs = linearProgram.getLHS();
   const auto& rhs = linearProgram.getRHS();

   const auto& varNames = linearProgram.getVarNames();
   const auto& consNames = linearProgram.getConsNames();

   // add variables to the model and build the objective expression
   for (size_t var = 0; var < linearProgram.getNCols(); ++var) {
      assert(var < varNames.size());
      variables.add(IloNumVar(env, lb[var], ub[var], varNames[var].c_str()));
      assert(var < obj.size());
      objExpr += variables[var] * obj[var];
   }

   model.add(IloMinimize(env, objExpr));

   for (size_t row = 0; row < linearProgram.getNRows(); ++row) {
      auto rowView = linearProgram.getRow(row);
      const size_t* indices = rowView.indices;

      IloNumExpr rowExpr(env);
      for (size_t id = 0; id < rowView.size; ++id)
         rowExpr += rowView.coefs[id] * variables[indices[id]];

      IloConstraint cons(lhs[row] <= rowExpr <= rhs[row]);
      cons.setName(consNames[row].c_str());
      constraints.add(lhs[row] <= rowExpr <= rhs[row]);
   }

   model.add(constraints);
   cplex = model;

   if (noOutput) cplex.setOut(env.getNullStream());
   cplex.exportModel("exp.mps");
}

LPResult CPXSolver::solve() {
   bool success = cplex.solve();

   LPResult result;
   auto cpxstatus = cplex.getStatus();

   if (cpxstatus == IloAlgorithm::Optimal) {
      assert(success);
      result.status = LPResult::OPTIMAL;

      IloNumArray vals(env);
      // get primal vals
      size_t ncols = this->mip.getNCols();
      cplex.getValues(vals, variables);
      for (size_t i = 0; i < ncols; ++i)
         result.primalSolution.push_back(vals[i]);

      // get dual vals
      size_t nrows = this->mip.getNRows();
      cplex.getDuals(vals, constraints);
      for (size_t i = 0; i < nrows; ++i) result.dualSolution.push_back(vals[i]);

      result.obj = cplex.getObjValue();
   } else if (cpxstatus == IloAlgorithm::Unbounded) {
      result.status = LPResult::UNBOUNDED;
   } else if (cpxstatus == IloAlgorithm::Infeasible) {
      result.status = LPResult::INFEASIBLE;
   } else {
      result.status = LPResult::OTHER;
   }

   return result;
}

/*
void CPXSolver::setIntParams(const std::initializer_list<IntParam>& params) {
   for (auto param : params) cplex.setParam(param.first, param.second);
}

void CPXSolver::setBoolParams(const std::initializer_list<BoolParam>& params)
{ for (auto param : params) cplex.setParam(param.first, param.second);
}*/

CPXSolver::~CPXSolver() { env.end(); }

#endif
