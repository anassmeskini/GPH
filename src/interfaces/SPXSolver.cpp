#include "SPXSolver.h"
#include "core/Common.h"

#ifdef SOPLEX_FOUND

SPXSolver::SPXSolver(const MIP& mip)
    : ncols(mip.getNCols()), nrows(mip.getNRows())
{
   using namespace soplex;

   const auto& lb = mip.getLB();
   const auto& ub = mip.getUB();
   const auto& obj = mip.getObj();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();

   mysoplex.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
   DSVector dummycol(0);

   const double myinf = std::numeric_limits<double>::infinity();

   for (int var = 0; var < mip.getNCols(); ++var)
   {
      double soLB = lb[var];
      double soUB = ub[var];
      if (soLB == -myinf)
         soLB = -infinity;
      if (soUB == myinf)
         soUB = infinity;

      mysoplex.addColReal(LPCol(obj[var], dummycol, soUB, soLB));
   }

   for (int row = 0; row < mip.getNRows(); ++row)
   {
      auto rowview = mip.getRow(row);

      DSVector vector(rowview.size);
      for (int id = 0; id < rowview.size; ++id)
         vector.add(rowview.indices[id], rowview.coefs[id]);

      double soLHS = lhs[row];
      double soRHS = rhs[row];
      if (soLHS == -myinf)
         soLHS = -infinity;
      if (soRHS == myinf)
         soRHS = infinity;

      mysoplex.addRowReal(LPRow(soLHS, vector, soRHS));
   }

   mysoplex.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
}

SPXSolver::SPXSolver(const SPXSolver& spxsolver)
    : mysoplex(spxsolver.mysoplex), ncols(spxsolver.ncols),
      nrows(spxsolver.nrows)
{
}

LPResult
SPXSolver::solve(Algorithm alg)
{
   using namespace soplex;

   switch (alg)
   {
   case Algorithm::PRIMAL:
      mysoplex.setIntParam(SoPlex::IntParam::ALGORITHM,
                           SoPlex::ALGORITHM_PRIMAL);
      break;
   case Algorithm::DUAL:
      mysoplex.setIntParam(SoPlex::IntParam::ALGORITHM,
                           SoPlex::ALGORITHM_DUAL);
      break;
   default:
      assert(0);
   }

   LPResult result;
   SPxSolver::Status stat;
   DVector prim(ncols);
   DVector dual(nrows);

   stat = mysoplex.optimize();

   if (stat == SPxSolver::OPTIMAL)
   {
      result.status = LPResult::OPTIMAL;
      mysoplex.getPrimalReal(prim);
      mysoplex.getDualReal(dual);

      // TODO use memcpy
      for (int i = 0; i < ncols; ++i)
         result.primalSol.push_back(prim[i]);

      for (int i = 0; i < nrows; ++i)
         result.dualSol.push_back(dual[i]);

      result.obj = mysoplex.objValueReal();
      result.niter = mysoplex.numIterations();
   }
   else if (stat == SPxSolver::INFEASIBLE)
      result.status = LPResult::INFEASIBLE;
   else
      result.status = LPResult::OTHER;

   return result;
}

std::unique_ptr<LPSolver>
SPXSolver::makeCopy() const
{
   return std::make_unique<SPXSolver>(*this);
}

void
SPXSolver::changeBounds(int column, double lb, double ub)
{
   mysoplex.changeBoundsReal(column, lb, ub);
}

void
SPXSolver::changeBounds(const std::vector<double>& lb,
                        const std::vector<double>& ub)
{
   using namespace soplex;

   DVector soplb(ncols);
   DVector sopub(ncols);

   for (int i = 0; i < ncols; ++i)
   {
      soplb[i] = lb[i];
      sopub[i] = ub[i];
   }

   mysoplex.changeBoundsReal(soplb, sopub);
}

void
SPXSolver::changeObjective(int column, double coef)
{
   mysoplex.changeObjReal(column, coef);
}

#endif
