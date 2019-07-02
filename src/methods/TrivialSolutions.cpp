#include "TrivialSolutions.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"

// TODO remove slacks
// TODO reorder variables

std::optional<std::vector<double>>
TrivialSolutions::search(const MIP<double>& problem,
                         const std::vector<double>& lb,
                         const std::vector<double>& ub,
                         const LPResult&,
                         const std::vector<double>&,
                         const std::vector<size_t>&,
                         const LPFactory&)
{
   size_t nrows = problem.getNRows();
   const auto& lhs = problem.getLHS();
   const auto& rhs = problem.getRHS();

   // copies
   std::vector<Activity> activities = computeActivities(problem);

   std::vector<double> solActivity(nrows, 0.0);

   auto lbcopy = lb;
   auto ubcopy = ub;

   // try solution = lb
   for (size_t row = 0; row < nrows; ++row)
   {
      auto rowview = problem.getRow(row);
      const double* rowcoefs = rowview.coefs;
      const size_t* rowindices = rowview.indices;
      size_t rowsize = rowview.size;

      for (size_t colid = 0; colid < rowsize; ++colid)
      {
         size_t col = rowindices[colid];
         double coef = rowcoefs[colid];

         if (lbcopy[col] != ub[col])
         {
            bool feasible =
              updateActivities<ChangedBound::UPPER>(problem.getCol(col),
                                                    ubcopy[col],
                                                    lbcopy[col],
                                                    activities,
                                                    lhs,
                                                    rhs);

            ubcopy[col] = lbcopy[col];

            if (!feasible)
            {
               Message::debug("unfeasible", row, nrows);
            }

            if (!propagate(problem, lbcopy, ubcopy, activities, col))
            {
               Message::debug("propagation failed!", row, nrows);
               return {};
            }
         }

         solActivity[row] += lbcopy[col] * coef;
      }

      if (Num::greater(solActivity[row], rhs[row]) ||
          Num::less(solActivity[row], lhs[row]))
      {
         Message::debug("LB sol infeasible, row {} ,nrows {}", row, nrows);
         break;
      }
   }

   auto objective = problem.getObj();

   double obj = 0.0;
   for (size_t i = 0; i < objective.size(); ++i)
      obj += objective[i] * lbcopy[i];

   Message::print("done, obj: {}", obj);

   return {};
}
