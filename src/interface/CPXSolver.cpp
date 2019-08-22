#include "CPXSolver.h"
#include "core/Common.h"
#include <iostream>

#include "io/Message.h"

#ifdef CONCERT_CPLEX_FOUND

CPXSolver::CPXSolver(const MIP& mip)
    : model(env), variables(env), constraints(env), ncols(mip.getNCols()),
      nrows(mip.getNRows())
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

      variables.add(
          IloNumVar(env, lb[var], ub[var], varNames[var].c_str()));
      objExpr += variables[var] * obj[var];
   }

   objective = IloMinimize(env, objExpr);
   model.add(objective);

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
   cplex.setOut(env.getNullStream());
}

CPXSolver::CPXSolver(const CPXSolver& cpxsolver)
    : env(), model(env), variables(env), constraints(env), cplex(env),
      ncols(cpxsolver.ncols), nrows(cpxsolver.nrows)
{
   model = cpxsolver.model.getClone(env);
   for (IloIterator<IloNumVar> it(env); it.ok(); ++it)
      variables.add(*it);

   for (IloIterator<IloRange> it(env); it.ok(); ++it)
      constraints.add(*it);

   cplex.extract(model);
   cplex.setOut(env.getNullStream());
}

LPResult
CPXSolver::solve(Algorithm alg)
{
   switch (alg)
   {
   case Algorithm::PRIMAL:
      cplex.setParam(IloCplex::Param::RootAlgorithm, CPX_ALG_PRIMAL);
      break;
   case Algorithm::DUAL:
      cplex.setParam(IloCplex::Param::RootAlgorithm, CPX_ALG_DUAL);
      break;
   default:
      assert(0);
   }

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
CPXSolver::changeBounds(int column, double lb, double ub)
{
   variables[column].setBounds(lb, ub);
   assert(variables[column].getLb() == lb);
   assert(variables[column].getUb() == ub);
}

void
CPXSolver::changeObjective(int column, double coef)
{
   objective.setLinearCoef(variables[column], coef);
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

CPXSolver::~CPXSolver() { env.end(); }

#endif
