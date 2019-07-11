#include "BoundSolution.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"

// TODO remove slacks
// TODO reorder variables

std::optional<Solution>
BoundSolution::search(const MIP& mip,
                      const std::vector<double>& lb,
                      const std::vector<double>& ub,
                      const std::vector<Activity>& activities,
                      const LPResult&,
                      const std::vector<double>&,
                      const std::vector<int>&,
                      std::shared_ptr<LPSolver>)
{
   int nrows = mip.getNRows();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();
   const auto& objective = mip.getObj();

   double cost = 0.0;

   // copies
   auto activities_copy = activities;
   std::vector<double> solActivity(nrows, 0.0);

   auto lbcopy = lb;
   auto ubcopy = ub;

   bool lbsolfeasible = true;
   // try solution = lb
   for (int row = 0; row < nrows && lbsolfeasible; ++row)
   {
      auto [rowcoefs, rowindices, rowsize] = mip.getRow(row);

      for (int colid = 0; colid < rowsize; ++colid)
      {
         int col = rowindices[colid];
         double coef = rowcoefs[colid];

         if (lbcopy[col] != ub[col])
         {
            bool feasible =
              updateActivities<ChangedBound::UPPER>(mip.getCol(col),
                                                    ubcopy[col],
                                                    lbcopy[col],
                                                    activities_copy,
                                                    lhs,
                                                    rhs);

            ubcopy[col] = lbcopy[col];

            if (!feasible)
            {
               Message::debug("unfeasible", row, nrows);
            }

            if (!propagate(mip, lbcopy, ubcopy, activities_copy, col))
            {
               Message::debug("propagation failed!", row, nrows);
               lbsolfeasible = false;
               break;
            }
         }

         solActivity[row] += lbcopy[col] * coef;
      }

      if (Num::greater(solActivity[row], rhs[row]) ||
          Num::less(solActivity[row], lhs[row]))
      {
         Message::debug("LB sol infeasible, row {} ,nrows {}", row, nrows);
         lbsolfeasible = false;
         break;
      }
   }

   if (lbsolfeasible)
   {
      for (size_t i = 0; i < objective.size(); ++i)
         cost += objective[i] * lbcopy[i];

      Message::print("lb sol feasible, obj: {}", cost);
   }

   // ub solution
   activities_copy = activities;

   solActivity = std::vector<double>(nrows, 0.0);

   lbcopy = lb;
   ubcopy = ub;

   bool ubsolfeasible = true;
   for (int row = 0; row < nrows && ubsolfeasible; ++row)
   {
      auto [rowcoefs, rowindices, rowsize] = mip.getRow(row);

      for (int colid = 0; colid < rowsize; ++colid)
      {
         int col = rowindices[colid];
         double coef = rowcoefs[colid];

         if (lbcopy[col] != ub[col])
         {
            bool feasible =
              updateActivities<ChangedBound::LOWER>(mip.getCol(col),
                                                    lbcopy[col],
                                                    ubcopy[col],
                                                    activities_copy,
                                                    lhs,
                                                    rhs);

            lbcopy[col] = ubcopy[col];

            if (!feasible)
            {
               Message::debug("unfeasible", row, nrows);
            }

            if (!propagate(mip, lbcopy, ubcopy, activities_copy, col))
            {
               Message::debug("propagation failed!", row, nrows);
               ubsolfeasible = false;
               break;
            }
         }

         solActivity[row] += ubcopy[col] * coef;
      }

      if (Num::greater(solActivity[row], rhs[row]) ||
          Num::less(solActivity[row], lhs[row]))
      {
         Message::debug("UB sol infeasible, row {} ,nrows {}", row, nrows);
         ubsolfeasible = false;
         break;
      }
   }

   if (ubsolfeasible)
   {
      cost = 0.0;
      for (size_t i = 0; i < objective.size(); ++i)
         cost += objective[i] * ubcopy[i];

      Message::print("ub sol feasible, obj: {}", cost);
   }

   return {};
}
