#include "CPXSolver.h"
#include <iostream>

#ifdef CONCERT_CPLEX_FOUND

CPXSolver::CPXSolver(const MIP<double>& mip)
  : model(env)
  , variables(env)
  , constraints(env)
  , ncols(mip.getNCols())
  , nrows(mip.getNRows())
  , deleteEnv(true)
{
  IloNumExpr objExpr(env);

  const auto& lb = mip.getLB();
  const auto& ub = mip.getUB();
  const auto& obj = mip.getObj();
  const auto& lhs = mip.getLHS();
  const auto& rhs = mip.getRHS();

  const auto& varNames = mip.getVarNames();
  const auto& consNames = mip.getConsNames();

  // add variables to the model and build the objective expression
  for (size_t var = 0; var < mip.getNCols(); ++var) {
    assert(var < varNames.size());
    variables.add(IloNumVar(env, lb[var], ub[var], varNames[var].c_str()));
    assert(var < obj.size());
    objExpr += variables[var] * obj[var];
  }

  model.add(IloMinimize(env, objExpr));

  for (size_t row = 0; row < mip.getNRows(); ++row) {
    auto rowView = mip.getRow(row);
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
}

CPXSolver::CPXSolver(const CPXSolver& cpxsolver)
  : env(cpxsolver.model.getEnv())
  , model(env)
  , variables(env)
  , constraints(env)
  , cplex(env)
  , ncols(cpxsolver.ncols)
  , nrows(cpxsolver.nrows)
  , deleteEnv(false)
{
  model = cpxsolver.model;
  variables = cpxsolver.variables;
  constraints = cpxsolver.constraints;
  cplex.extract(model);
}

LPResult
CPXSolver::solve()
{
  cplex.setOut(env.getNullStream());

  bool success = cplex.solve();

  LPResult result;
  auto cpxstatus = cplex.getStatus();

  if (cpxstatus == IloAlgorithm::Optimal) {
    assert(success);
    result.status = LPResult::OPTIMAL;

    IloNumArray vals(env);
    // get primal vals
    cplex.getValues(vals, variables);
    for (size_t i = 0; i < ncols; ++i)
      result.primalSolution.push_back(vals[i]);

    // get dual vals
    cplex.getDuals(vals, constraints);
    for (size_t i = 0; i < nrows; ++i)
      result.dualSolution.push_back(vals[i]);

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

std::unique_ptr<LPSolver<double>>
CPXSolver::clone() const
{
  return std::unique_ptr<LPSolver<double>>(new CPXSolver(*this));
}

/*
void CPXSolver::setIntParams(const std::initializer_list<IntParam>& params) {
   for (auto param : params) cplex.setParam(param.first, param.second);
}

void CPXSolver::setBoolParams(const std::initializer_list<BoolParam>& params)
{ for (auto param : params) cplex.setParam(param.first, param.second);
}*/

CPXSolver::~CPXSolver()
{
  if (deleteEnv)
    env.end();
}

#endif
