#include "Propagation.h"
#include "Common.h"
#include "Numerics.h"
#include "io/Message.h"

#include <memory>
template<>
bool
updateActivities<ChangedBound::LOWER>(VectorView<double> colview,
                                      double oldlb,
                                      double newlb,
                                      std::vector<Activity>& activities,
                                      const std::vector<double>& lhs,
                                      const std::vector<double>& rhs)
{
   const double* colcoefs = colview.coefs;
   const size_t* colindices = colview.indices;
   const size_t colsize = colview.size;

   assert(!Num::isInf(newlb));
   bool lbfinite = !Num::isInf(oldlb);

   for (size_t i = 0; i < colsize; ++i)
   {
      size_t row = colindices[i];
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

      if ((!Num::less(activities[row].min, rhs[row]) &&
           !activities[row].ninfmin) ||
          (!Num::greater(activities[row].max, lhs[row]) &&
           !activities[row].ninfmax))
      {
         Message::debug("Act update: row {} activities {} {}, sides {} {}",
                        row,
                        activities[row].min,
                        activities[row].max,
                        lhs[row],
                        rhs[row]);
         return false;
      }
   }

   return true;
}

template<>
bool
updateActivities<ChangedBound::UPPER>(VectorView<double> colview,
                                      double oldub,
                                      double newub,
                                      std::vector<Activity>& activities,
                                      const std::vector<double>& lhs,
                                      const std::vector<double>& rhs)
{
   const double* colcoefs = colview.coefs;
   const size_t* colindices = colview.indices;
   const size_t colsize = colview.size;

   assert(!Num::isInf(newub));
   bool ubfinite = !Num::isInf(oldub);

   for (size_t i = 0; i < colsize; ++i)
   {
      size_t row = colindices[i];
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

      if ((!Num::less(activities[row].min, rhs[row]) &&
           !activities[row].ninfmin) ||
          (!Num::greater(activities[row].max, lhs[row]) &&
           !activities[row].ninfmax))
      {
         Message::debug("Act update: row {} activities {} {}, sides {} {}",
                        row,
                        activities[row].min,
                        activities[row].max,
                        lhs[row],
                        rhs[row]);
         return false;
      }
   }

   return true;
}

bool
updateActivities(VectorView<double> colview,
                 double oldlb,
                 double newlb,
                 double oldub,
                 double newub,
                 std::vector<Activity>& activities,
                 const std::vector<double>& lhs,
                 const std::vector<double>& rhs)
{
   const double* colcoefs = colview.coefs;
   const size_t* colindices = colview.indices;
   const size_t colsize = colview.size;

   assert(!Num::isMinusInf(newlb) && !Num::isInf(newub));

   bool lbfinite = !Num::isInf(oldlb);
   bool ubfinite = !Num::isInf(oldub);

   for (size_t i = 0; i < colsize; ++i)
   {
      size_t row = colindices[i];
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

      if ((!Num::less(activities[row].min, rhs[row]) &&
           !activities[row].ninfmin) ||
          (!Num::greater(activities[row].max, lhs[row]) &&
           !activities[row].ninfmax))
      {
         Message::debug("Act update: row {} activities {} {}, sides {} {}",
                        row,
                        activities[row].min,
                        activities[row].max,
                        lhs[row],
                        rhs[row]);
         return false;
      }
   }

   return true;
}

static bool
propagateRow(const ProblemView& problem,
             size_t row,
             std::vector<Activity>& activities,
             std::vector<double>& lb,
             std::vector<double>& ub)
{
   auto rowview = problem.getRow(row);
   const double* rowcoefs = rowview.coefs;
   const size_t* rowindices = rowview.indices;
   const size_t rowsize = rowview.size;

   auto& activity = activities[row];

   auto& integer = problem.integer();

   auto& lhs = problem.getLHS();
   auto& rhs = problem.getRHS();

   double impliedlb;
   double impliedub;

   for (size_t i = 0; i < rowsize; ++i)
   {
      size_t col = rowindices[i];
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
         Message::debug(
           "Propagation changed Bs of col {} from [{}, {}] -> {{}, {}]",
           col,
           lb[col],
           ub[col],
           impliedlb,
           impliedub);
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
      }
      else if (Num::greater(impliedlb, lb[col]) && implbfinite)
      {
         Message::debug("Propagation changed LB of col {} from {} -> {}",
                        col,
                        lb[col],
                        impliedlb);
         // update left
         if (!updateActivities<ChangedBound::LOWER>(
               colview, lb[col], impliedlb, activities, lhs, rhs))
            return false;

         lb[col] = impliedlb;
      }
      else if (Num::less(impliedub, ub[col]) && impubfinite)
      {
         Message::debug("Propagation changed UB of col {}: {} --> {}",
                        col,
                        ub[col],
                        impliedub);
         // update right
         if (!updateActivities<ChangedBound::UPPER>(
               colview, ub[col], impliedub, activities, lhs, rhs))
            return false;

         ub[col] = impliedub;
      }
   }

   return true;
}

// assumes that the lb and ub reftlect the fixing, and the activities are
// up-to-date
int
propagate(const ProblemView& problem,
          std::vector<double>& lb,
          std::vector<double>& ub,
          std::vector<Activity>& activities,
          size_t changedcol)
{
   assert(lb[changedcol] == ub[changedcol]);

   size_t nrows = problem.getNRows();
   size_t ncols = problem.getNCols();

   std::vector<size_t> changedCols;
   bitset colChanged(ncols, false);
   changedCols.reserve(nrows);

   changedCols.push_back(changedcol);
   colChanged[changedcol] = true;

   for (size_t i = 0; i < changedCols.size(); ++i)
   {
      const size_t col = changedCols[i];

      assert(colChanged[col]);

      auto colview = problem.getCol(col);
      const size_t* colindices = colview.indices;
      const size_t colsize = colview.size;

      for (size_t j = 0; j < colsize; ++j)
      {
         const size_t row = colindices[j];

         // propagate row
         // will update the candidate list through the callback
         if (!propagateRow(problem, row, activities, lb, ub))
            return false;
      }

      colChanged[col] = false;

      if (changedCols.size() >= ncols)
      {
         int newsize = static_cast<int>(changedCols.size()) - (i + 1);
         assert(newsize >= 0);

         std::memmove(changedCols.data(), changedCols.data() + i + 1, newsize);
         changedCols.resize(newsize);
         i = 0;
      }
   }

   return 0;
}
