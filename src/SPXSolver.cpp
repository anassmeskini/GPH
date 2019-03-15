#include "SPXSolver.h"

#ifdef SOPLEX_FOUND

SPXSolver::SPXSolver(const MIP<double>& mip) : LPSolver(mip) {
   using namespace soplex;

   const auto& lb = mip.getLB();
   const auto& ub = mip.getUB();
   const auto& obj = mip.getObj();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();

   mysoplex.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
   DSVector dummycol(0);

   for (size_t var = 0; var < mip.getNCols(); ++var)
      mysoplex.addColReal(LPCol(obj[var], dummycol, lb[var], ub[var]));

   for (size_t row = 0; row < mip.getNRows(); ++row) {
      auto rowview = mip.getRow(row);

      DSVector vector(rowview.size);
      for (size_t id = 0; id < rowview.size; ++id)
         vector.add(rowview.indices[id], rowview.coefs[id]);

      mysoplex.addRowReal(LPRow(lhs[row], vector, rhs[row]));
   }
}

LPResult SPXSolver::solve() {
   using namespace soplex;

   LPResult result;
   SPxSolver::Status stat;
   DVector prim(mip.getNCols());
   DVector dual(mip.getNRows());

   stat = mysoplex.optimize();

   if (stat == SPxSolver::OPTIMAL) {
      result.status = LPResult::OPTIMAL;
      mysoplex.getPrimalReal(prim);
      mysoplex.getDualReal(dual);

      // TODO use memcpy
      for (size_t i = 0; i < mip.getNRows(); ++i)
         result.primalSolution.push_back(prim[i]);

      for (size_t i = 0; i < mip.getNCols(); ++i)
         result.primalSolution.push_back(dual[i]);
   } else
      result.status = LPResult::OTHER;

   return result;
}

SPXSolver::~SPXSolver() {}

#endif
