#include "MinLockRounding.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"

std::optional<Solution>
MinLockRounding::search(const MIP& mip,
                        const std::vector<double>& lb,
                        const std::vector<double>& ub,
                        const std::vector<Activity>&,
                        const LPResult& result,
                        const std::vector<double>& solAct,
                        const std::vector<int>&,
                        std::shared_ptr<LPSolver>)
{
   int nrows = mip.getNRows();
   int ncols = mip.getNCols();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();
   const auto& upLocks = mip.getUpLocks();
   const auto& downLocks = mip.getDownLocks();
   const auto& objective = mip.getObj();
   const auto& integer = mip.getInteger();

   double cost = 0.0;
   std::vector<double> solActivity = solAct;
   auto lbcopy = lb;
   auto ubcopy = ub;

   auto colPermutation = getIndentity(ncols);
   auto solution = result.primalSolution;

   std::sort(std::begin(colPermutation),
             std::end(colPermutation),
             [&](int lhs, int rhs) {
                return (integer[lhs] && !integer[rhs]) ||
                       std::min(downLocks[lhs], upLocks[lhs]) <
                         std::min(downLocks[rhs], upLocks[rhs]);
             });

   for (int i = 0; i < ncols; ++i)
   {
      int col = colPermutation[i];
      if (integer[col])
      {
         if (downLocks[col] < upLocks[col])
            solution[col] = Num::floor(solution[col]);
         else
            solution[col] = Num::ceil(solution[col]);
      }
   }

   // row and activity
   std::vector<std::pair<int, double>> violatedRows;

   for (int row = 0; row < nrows; ++row)
   {
      auto [coefs, indices, size] = mip.getRow(row);
      double activity = 0.0;

      for (int i = 0; i < size; ++i)
         activity += coefs[i] * solution[indices[i]];

      if (Num::greater(activity, rhs[row]) || Num::less(activity, lhs[row]))
         violatedRows.emplace_back(row, activity);
   }

   if (violatedRows.empty())
   {
      for (size_t i = 0; i < objective.size(); ++i)
         cost += objective[i] * lbcopy[i];

      Message::debug("lb sol feasible, obj: {}", cost);
   }
   else
      Message::debug("violated rows: {}", violatedRows.size());

   return {};
}
