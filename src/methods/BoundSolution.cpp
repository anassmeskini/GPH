#include "BoundSolution.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"

// TODO fix binary first

void
BoundSolution::search(const MIP& mip, const std::vector<double>& lb,
                      const std::vector<double>& ub,
                      const std::vector<Activity>& activities,
                      const LPResult&, const std::vector<double>&,
                      const std::vector<int>&,
                      std::shared_ptr<const LPSolver> solver,
                      SolutionPool& pool)
{
   int ncols = mip.getNCols();
   const auto& integer = mip.getInteger();
   const auto& objective = mip.getObj();

   std::unique_ptr<LPSolver> localsolver;

   // try solution = lower bound
   auto local_activities = activities;
   auto locallb = lb;
   auto localub = ub;
   bool feasible = true;

   std::vector<int> inflbIntVars;
   for (int col = 0; col < ncols; ++col)
   {
      if (integer[col] && locallb[col] != localub[col])
      {
         if (Num::isMinusInf(locallb[col]))
         {
            inflbIntVars.push_back(col);
            continue;
         }

         double oldub = localub[col];
         localub[col] = locallb[col];

         if (!propagate(mip, locallb, localub, local_activities, col,
                        oldub, locallb[col]))
         {
            feasible = false;
            break;
         }
      }
   }

   for (auto col : inflbIntVars)
   {
      double oldlb = locallb[col];
      double oldub = localub[col];
      // check if the bound have been changed by propagation
      if (!Num::isMinusInf(oldlb))
      {
         localub[col] = locallb[col];
      }
      else
      {
         if (Num::isInf(localub[col]))
         {
            locallb[col] = 0.0;
            localub[col] = 0.0;
         }
         else
            locallb[col] = localub[col];

         if (!propagate(mip, locallb, localub, local_activities, col,
                        oldlb, oldub))
         {
            feasible = false;
            break;
         }
      }
   }

   if (feasible)
   {
      if (mip.getStatistics().ncont == 0)
      {
         Message::debug("Bnd: lb sol feasible");

         // TODO compute the objective while fixing
         double obj = 0.0;
         for (int i = 0; i < ncols; ++i)
            obj += objective[i] * locallb[i];

         pool.add(std::move(locallb), obj);
      }
      else
      {
         Message::debug("Bnd: solving local lp");

         if (!localsolver)
            localsolver = solver->clone();
         localsolver->changeBounds(locallb, localub);

         auto localresult = localsolver->solve(Algorithm::DUAL);
         if (localresult.status == LPResult::OPTIMAL)
         {
            Message::debug("Bnd: lb: lp feasible");
            pool.add(std::move(localresult.primalSolution),
                     localresult.obj);
         }
         else if (localresult.status == LPResult::INFEASIBLE)
            Message::debug("Bnd: lb: lp infeasible");
      }
   }

   // try solution = upper bound
   local_activities = activities;
   locallb = lb;
   localub = ub;
   feasible = true;

   std::vector<int> infubIntVars;
   for (int col = 0; col < ncols; ++col)
   {
      if (integer[col] && locallb[col] != localub[col])
      {
         if (Num::isInf(localub[col]))
         {
            infubIntVars.push_back(col);
            continue;
         }

         double oldlb = locallb[col];
         locallb[col] = localub[col];

         if (!propagate(mip, locallb, localub, local_activities, col,
                        oldlb, localub[col]))
         {
            feasible = false;
            break;
         }
      }
   }

   for (auto col : infubIntVars)
   {
      double oldlb = locallb[col];
      double oldub = localub[col];
      // check if the bound have been changed by propagation
      if (!Num::isInf(oldub))
      {
         locallb[col] = localub[col];
      }
      else
      {
         if (Num::isMinusInf(locallb[col]))
         {
            locallb[col] = 0.0;
            localub[col] = 0.0;
         }
         else
            localub[col] = locallb[col];

         if (!propagate(mip, locallb, localub, local_activities, col,
                        oldlb, oldub))
         {
            feasible = false;
            break;
         }
      }
   }

   if (feasible)
   {
      if (mip.getStatistics().ncont == 0)
      {
         Message::debug("Bnd: ub sol feasible");
         double obj = 0.0;

         for (int i = 0; i < ncols; ++i)
            obj += objective[i] * locallb[i];

         pool.add(std::move(locallb), obj);
      }
      else
      {
         Message::debug("Bnd: ub: solving local lp");

         if (!localsolver)
            localsolver = solver->clone();
         localsolver->changeBounds(locallb, localub);

         auto localresult = localsolver->solve(Algorithm::DUAL);
         if (localresult.status == LPResult::OPTIMAL)
         {
            Message::debug("Bnd: ub: lp feasible, obj {:0.2e}",
                           localresult.obj);
            pool.add(std::move(localresult.primalSolution),
                     localresult.obj);
         }
         else if (localresult.status == LPResult::INFEASIBLE)
            Message::debug("Bnd: ub: lp infeasible");
      }
   }
}
