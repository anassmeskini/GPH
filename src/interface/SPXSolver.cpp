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

   const double myinf = std::numeric_limits<double>::infinity();

   for (size_t var = 0; var < mip.getNCols(); ++var) {
      double soLB = lb[var];
      double soUB = ub[var];
      if (soLB == -myinf) soLB = -infinity;
      if (soUB == myinf) soUB = infinity;

      mysoplex.addColReal(LPCol(obj[var], dummycol, soUB, soLB));
   }

   for (size_t row = 0; row < mip.getNRows(); ++row) {
      auto rowview = mip.getRow(row);

      DSVector vector(rowview.size);
      for (size_t id = 0; id < rowview.size; ++id)
         vector.add(rowview.indices[id], rowview.coefs[id]);

      double soLHS = lhs[row];
      double soRHS = rhs[row];
      if (soLHS == -myinf) soLHS = -infinity;
      if (soRHS == myinf) soRHS = infinity;

      mysoplex.addRowReal(LPRow(soLHS, vector, soRHS));
   }

   mysoplex.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
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
      for (size_t i = 0; i < mip.getNCols(); ++i)
         result.primalSolution.push_back(prim[i]);

      for (size_t i = 0; i < mip.getNRows(); ++i)
         result.dualSolution.push_back(dual[i]);

      result.obj = mysoplex.objValueReal();
   } else
      result.status = LPResult::OTHER;

   return result;
}

SPXSolver::~SPXSolver() {}

#endif
