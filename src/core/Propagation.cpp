#include "Propagation.h"
#include "Common.h"
#include "Numerics.h"
#include "io/Message.h"

template<>
bool
updateActivities<ChangedBound::LOWER>(VectorView colview,
                                      double oldlb,
                                      double newlb,
                                      std::vector<Activity>& activities,
                                      const std::vector<double>& lhs,
                                      const std::vector<double>& rhs)
{
   auto [colcoefs, colindices, colsize] = colview;
   bool lbfinite = !Num::isInf(oldlb);
   assert(!Num::isInf(newlb));

   for (int i = 0; i < colsize; ++i)
   {
      int row = colindices[i];
      const double coef = colcoefs[i];

      if (Num::greater(coef, 0.0))
      {
         if (lbfinite)
            activities[row].min += coef * (newlb - oldlb);
         else
         {
            activities[row].min += coef * newlb;
            --activities[row].ninfmin;
         }
      }
      else
      {
         if (lbfinite)
            activities[row].max += coef * (newlb - oldlb);
         else
         {
            activities[row].max += coef * newlb;
            --activities[row].ninfmax;
         }
      }

      if ((Num::greater(activities[row].min, rhs[row]) &&
           !activities[row].ninfmin) ||
          (Num::less(activities[row].max, lhs[row]) &&
           !activities[row].ninfmax))
      {
         /*Message::debug("Act update infeasiblity: row {} activities {} {},
            sides {} {}", row, activities[row].min, activities[row].max,
                        lhs[row],
                        rhs[row]);*/
         return false;
      }
   }

   return true;
}

template<>
bool
updateActivities<ChangedBound::UPPER>(VectorView colview,
                                      double oldub,
                                      double newub,
                                      std::vector<Activity>& activities,
                                      const std::vector<double>& lhs,
                                      const std::vector<double>& rhs)
{
   const double* colcoefs = colview.coefs;
   const int* colindices = colview.indices;
   const int colsize = colview.size;

   assert(!Num::isInf(newub));
   bool ubfinite = !Num::isInf(oldub);

   for (int i = 0; i < colsize; ++i)
   {
      int row = colindices[i];
      const double coef = colcoefs[i];

      if (Num::greater(coef, 0.0))
      {
         if (ubfinite)
            activities[row].max += coef * (newub - oldub);
         else
         {
            activities[row].max += coef * newub;
            --activities[row].ninfmax;
         }
      }
      else
      {
         if (ubfinite)
            activities[row].min += coef * (newub - oldub);
         else
         {
            activities[row].min += coef * newub;
            --activities[row].ninfmin;
         }
      }

      // TODO check if activity.min > activity.max
      if ((Num::greater(activities[row].min, rhs[row]) &&
           !activities[row].ninfmin) ||
          (Num::less(activities[row].max, lhs[row]) &&
           !activities[row].ninfmax))
      {
         /*Message::debug(
           "Act update infeasility: row {} activities {} {}, sides {} {}",
           row,
           activities[row].min,
           activities[row].max,
           lhs[row],
           rhs[row]);*/
         return false;
      }
   }

   return true;
}

bool
updateActivities(VectorView colview,
                 double oldlb,
                 double newlb,
                 double oldub,
                 double newub,
                 std::vector<Activity>& activities,
                 const std::vector<double>& lhs,
                 const std::vector<double>& rhs)
{
   const double* colcoefs = colview.coefs;
   const int* colindices = colview.indices;
   const int colsize = colview.size;

   assert(!Num::isMinusInf(newlb) && !Num::isInf(newub));

   bool lbfinite = !Num::isInf(oldlb);
   bool ubfinite = !Num::isInf(oldub);

   for (int i = 0; i < colsize; ++i)
   {
      int row = colindices[i];
      const double coef = colcoefs[i];

      if (Num::greater(coef, 0.0))
      {
         if (lbfinite)
            activities[row].min += coef * (newlb - oldlb);
         else
         {
            activities[row].min += coef * newlb;
            --activities[row].ninfmin;
         }

         if (ubfinite)
            activities[row].max += coef * (newub - oldub);
         else
         {
            activities[row].max += coef * newub;
            --activities[row].ninfmax;
         }
      }
      else
      {
         if (ubfinite)
            activities[row].min += coef * (newub - oldub);
         else
         {
            activities[row].min += coef * newub;
            --activities[row].ninfmin;
         }

         if (lbfinite)
            activities[row].max += coef * (newlb - oldlb);
         else
         {
            activities[row].max += coef * newlb;
            --activities[row].ninfmax;
         }
      }

      if ((Num::greater(activities[row].min, rhs[row]) &&
           !activities[row].ninfmin) ||
          (Num::less(activities[row].max, lhs[row]) &&
           !activities[row].ninfmax))
      {
         /*Message::debug(
           "Act update infeasibility: row {} activities {} {}, sides {} {}",
           row,
           activities[row].min,
           activities[row].max,
           lhs[row],
           rhs[row]);*/
         return false;
      }
   }

   return true;
}

static bool
propagateRow(const MIP& problem,
             int row,
             std::vector<Activity>& activities,
             std::vector<double>& lb,
             std::vector<double>& ub,
             std::vector<int>& changedCols)
{
   auto [rowcoefs, rowindices, rowsize] = problem.getRow(row);

   auto& activity = activities[row];

   const auto& integer = problem.getInteger();

   const auto& lhs = problem.getLHS();
   const auto& rhs = problem.getRHS();

   double impliedlb;
   double impliedub;

   for (int i = 0; i < rowsize; ++i)
   {
      int col = rowindices[i];
      double coef = rowcoefs[i];

      double maxPartialActivity = activity.max;
      double minPartialActivity = activity.min;

      bool implbfinite = false;
      bool impubfinite = false;

      const bool ubInf = Num::isInf(ub[col]);
      const bool lbInf = Num::isMinusInf(lb[col]);

      if (Num::greater(coef, 0.0))
      {
         if (!ubInf)
            maxPartialActivity -= coef * ub[col];

         if (!lbInf)
            minPartialActivity -= coef * lb[col];

         if (integer[col])
         {
            impliedlb = Num::ceil((lhs[row] - maxPartialActivity) / coef);
            impliedub = Num::floor((rhs[row] - minPartialActivity) / coef);
         }
         else
         {
            impliedlb = (lhs[row] - maxPartialActivity) / coef;
            impliedub = (rhs[row] - minPartialActivity) / coef;
         }

         implbfinite =
           activity.ninfmax == 0 || (activity.ninfmax == 1 && ubInf);
         impubfinite =
           activity.ninfmin == 0 || (activity.ninfmin == 1 && lbInf);
      }
      else
      {
         assert(coef != 0.0);

         if (!lbInf)
            maxPartialActivity -= coef * lb[col];
         if (!ubInf)
            minPartialActivity -= coef * ub[col];

         if (integer[col])
         {
            impliedlb = Num::ceil((rhs[row] - minPartialActivity) / coef);
            impliedub = Num::floor((lhs[row] - maxPartialActivity) / coef);
         }
         else
         {
            impliedlb = (rhs[row] - minPartialActivity) / coef;
            impliedub = (lhs[row] - maxPartialActivity) / coef;
         }

         implbfinite =
           activity.ninfmin == 0 || (activity.ninfmin == 1 && lbInf);
         impubfinite =
           activity.ninfmax == 0 || (activity.ninfmax == 1 && ubInf);
      }

      auto colview = problem.getCol(col);

      if ((Num::greater(impliedlb, lb[col]) && implbfinite) &&
          (Num::less(impliedub, ub[col]) && impubfinite))
      {
         /*Message::debug(
           "Propagation changed Bds of col {} from [{}, {}] -> {{}, {}]",
           col,
           lb[col],
           ub[col],
           impliedlb,
           impliedub);*/
         // update right and left
         if (!updateActivities(colview,
                               lb[col],
                               impliedlb,
                               ub[col],
                               impliedub,
                               activities,
                               lhs,
                               rhs))
            return false;

         lb[col] = impliedlb;
         ub[col] = impliedub;
         changedCols.push_back(col);
      }
      else if (Num::greater(impliedlb, lb[col]) && implbfinite)
      {
         /*Message::debug("Propagation changed LB of col {} from {} -> {}",
                        col,
                        lb[col],
                        impliedlb);*/
         // update left
         if (!updateActivities<ChangedBound::LOWER>(
               colview, lb[col], impliedlb, activities, lhs, rhs))
            return false;

         lb[col] = impliedlb;
         changedCols.push_back(col);
      }
      else if (Num::less(impliedub, ub[col]) && impubfinite)
      {
         /*Message::debug("Propagation changed UB of col {}: {} --> {}",
                        col,
                        ub[col],
                        impliedub);*/
         // update right
         if (!updateActivities<ChangedBound::UPPER>(
               colview, ub[col], impliedub, activities, lhs, rhs))
            return false;

         ub[col] = impliedub;
         changedCols.push_back(col);
      }
   }

   return true;
}

// assumes that the lb and ub reftlect the fixing, and the activities are
// up-to-date
int
propagate(const MIP& problem,
          std::vector<double>& lb,
          std::vector<double>& ub,
          std::vector<Activity>& activities,
          int changedcol)
{
   assert(lb[changedcol] == ub[changedcol]);

   int ncols = problem.getNCols();

   std::vector<int> changedCols;
   // is this a bug?
   // changedCols.reserve(ncols);

   changedCols.push_back(changedcol);

   for (size_t i = 0; i < changedCols.size(); ++i)
   {
      const int col = changedCols[i];

      auto [colcoefs, colindices, colsize] = problem.getCol(col);

      for (int j = 0; j < colsize; ++j)
      {
         const int row = colindices[j];

         // propagate row
         if (!propagateRow(problem, row, activities, lb, ub, changedCols))
            return 0;
      }

      if (changedCols.size() >= static_cast<size_t>(ncols))
      {
         int newsize = static_cast<int>(changedCols.size()) - (i + 1);
         assert(newsize >= 0);

         std::memmove(changedCols.data(), changedCols.data() + i + 1, newsize);
         changedCols.resize(newsize);
         i = 0;
      }
   }

   return 1;
}
