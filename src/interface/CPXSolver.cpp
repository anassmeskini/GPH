#include "CPXSolver.h"
#include "core/Common.h"
#include <iostream>

#ifdef CONCERT_CPLEX_FOUND

CPXSolver::CPXSolver(const MIP& mip)
  : model(env)
  , variables(env)
  , constraints(env)
  , deleteEnv(true)
  , ncols(mip.getNCols())
  , nrows(mip.getNRows())
{
   IloNumExpr objExpr(env);

   const auto& lb = mip.getLB();
   const auto& ub = mip.getUB();
   const auto& obj = mip.getObj();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();

   const auto& varNames = mip.getVarNames();

   // add variables to the model and build the objective expression
   for (int var = 0; var < mip.getNCols(); ++var)
   {
      assert(static_cast<size_t>(var) < varNames.size());
      assert(static_cast<size_t>(var) < obj.size());

      variables.add(IloNumVar(env, lb[var], ub[var], varNames[var].c_str()));
      objExpr += variables[var] * obj[var];
   }

   model.add(IloMinimize(env, objExpr));

   for (int row = 0; row < mip.getNRows(); ++row)
   {
      auto rowView = mip.getRow(row);
      const int* indices = rowView.indices;

      IloNumExpr rowExpr(env);
      for (int id = 0; id < rowView.size; ++id)
         rowExpr += rowView.coefs[id] * variables[indices[id]];

      // IloConstraint cons(lhs[row] <= rowExpr <= rhs[row]);
      constraints.add(lhs[row] <= rowExpr <= rhs[row]);
   }

   model.add(constraints);
   try
   {
      cplex = model;
   }
   catch (IloAlgorithm::CannotExtractException& ex)
   {
      // TODO
      assert(0);
   }
}

CPXSolver::CPXSolver(const CPXSolver& cpxsolver)
  : env(cpxsolver.model.getEnv())
  , model(env)
  , variables(env)
  , constraints(env)
  , cplex(env)
  , deleteEnv(false)
  , ncols(cpxsolver.ncols)
  , nrows(cpxsolver.nrows)
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

   cplex.solve();

   LPResult result;
   auto cpxstatus = cplex.getStatus();

   if (cpxstatus == IloAlgorithm::Optimal)
   {
      result.status = LPResult::OPTIMAL;

      IloNumArray vals(env);
      // get primal vals
      cplex.getValues(vals, variables);
      for (int i = 0; i < ncols; ++i)
         result.primalSolution.push_back(vals[i]);

      // get dual vals
      cplex.getDuals(vals, constraints);
      for (int i = 0; i < nrows; ++i)
         result.dualSolution.push_back(vals[i]);

      result.obj = cplex.getObjValue();
   }
   else if (cpxstatus == IloAlgorithm::Unbounded)
   {
      result.status = LPResult::UNBOUNDED;
   }
   else if (cpxstatus == IloAlgorithm::Infeasible)
   {
      result.status = LPResult::INFEASIBLE;
   }
   else
   {
      result.status = LPResult::OTHER;
   }

   return result;
}

std::unique_ptr<LPSolver>
CPXSolver::makeCopy() const
{
   return std::make_unique<CPXSolver>(*this);
}

void
CPXSolver::branch(int column, double val, Direction direction)
{
   IloNumExpr rowExpr(env);
   rowExpr += 1.0 * variables[column];

   constexpr double inf = std::numeric_limits<double>::infinity();
   if (direction == Direction::UP)
      model.add(IloConstraint(val <= rowExpr <= inf));
   else
      model.add(IloConstraint(-inf <= rowExpr <= val));
}

void
CPXSolver::changeBounds(int column, double lb, double ub)
{
   variables[column].setBounds(lb, ub);
}

void
CPXSolver::changeBounds(const std::vector<double>& lb,
                        const std::vector<double>& ub)
{
   IloNumArray ilolb(env, ncols);
   IloNumArray iloub(env, ncols);

   for (int i = 0; i < ncols; ++i)
   {
      ilolb[i] = lb[i];
      iloub[i] = ub[i];
   }

   variables.setBounds(ilolb, iloub);
}

CPXSolver::~CPXSolver()
{
   if (deleteEnv)
      env.end();
}

#endif
