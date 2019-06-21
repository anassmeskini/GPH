#include "TrivialSolutions.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"

#include <math.h>

// TODO remove slacks

std::optional<std::vector<double>>
TrivialSolutions::search(const ProblemView& problem)
{
   size_t ncols = problem.getNCols();
   size_t nrows = problem.getNRows();

   const auto& lhs = problem.getLHS();
   const auto& rhs = problem.getRHS();

   const auto& downLocks = problem.getDownLocks();
   const auto& upLocks = problem.getUpLocks();

   std::vector<Activity> activities = problem.getActivities();
   std::vector<double> lb = problem.getLB();
   std::vector<double> ub = problem.getUB();

   std::vector<double> solActivity(nrows, 0.0);

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

         if (lb[col] != ub[col])
         {
            bool feasible = updateActivities<ChangedBound::UPPER>(
              problem.getCol(col), ub[col], lb[col], activities, lhs, rhs);

            ub[col] = lb[col];

            if (!feasible)
            {
               Message::debug("unfeasible", row, nrows);
            }

            if (!propagate(problem, lb, ub, activities, col))
            {
               Message::debug("propagation failed!", row, nrows);
               return {};
            }
         }

         solActivity[row] += lb[col] * coef;
      }

      if (Num::greater(solActivity[row], rhs[row]) ||
          Num::less(solActivity[row], lhs[row]))
      {
         Message::debug("LB sol infeasible, row {} ,nrows {}", row, nrows);
         break;
      }
   }

   Message::print("done");

   return {};
}
