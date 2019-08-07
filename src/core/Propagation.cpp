#include "Propagation.h"
#include "Common.h"
#include "Numerics.h"
#include "io/Message.h"

template <>
void
updateActivities<ChangedBound::LOWER>(
    VectorView colview, double oldlb, double newlb,
    std::vector<Activity>& activities) noexcept
{
   auto [colcoefs, colindices, colsize] = colview;
   bool lbfinite = !Num::isInf(oldlb);
   assert(!Num::isInf(newlb));

   for (int i = 0; i < colsize; ++i)
   {
      int row = colindices[i];
      const double coef = colcoefs[i];

      if (coef > 0.0)
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
   }
}

template <>
void
updateActivities<ChangedBound::UPPER>(
    VectorView colview, double oldub, double newub,
    std::vector<Activity>& activities) noexcept
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

      if (coef > 0.0)
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
   }
}

void
updateActivities(VectorView colview, double oldlb, double newlb,
                 double oldub, double newub,
                 std::vector<Activity>& activities) noexcept
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

      if (coef > 0.0)
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
   }
}

static bool
propagateRow(const MIP& problem, int row,
             std::vector<Activity>& activities, std::vector<double>& lb,
             std::vector<double>& ub,
             std::vector<int>& changedCols) noexcept
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

      if (coef > 0.0)
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
             (!Num::isMinusInf(lhs[row])) &&
             (activity.ninfmax == 0 || (activity.ninfmax == 1 && ubInf));
         impubfinite =
             (!Num::isInf(rhs[row])) &&
             (activity.ninfmin == 0 || (activity.ninfmin == 1 && lbInf));
      }
      else
      {
         assert(coef != 0.0);

         if (!ubInf)
            minPartialActivity -= coef * ub[col];
         if (!lbInf)
            maxPartialActivity -= coef * lb[col];

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
             (!Num::isInf(rhs[row])) &&
             (activity.ninfmin == 0 || (activity.ninfmin == 1 && ubInf));
         impubfinite =
             (!Num::isMinusInf(lhs[row])) &&
             (activity.ninfmax == 0 || (activity.ninfmax == 1 && lbInf));
      }

      auto colview = problem.getCol(col);

      if (((impliedlb > lb[col]) && implbfinite) &&
          ((impliedub < ub[col]) && impubfinite))
      {
         // update right and left
         updateActivities(colview, lb[col], impliedlb, ub[col], impliedub,
                          activities);

         lb[col] = impliedlb;
         ub[col] = impliedub;
         changedCols.push_back(col);
      }
      else if ((impliedlb > lb[col]) && implbfinite)
      {
         // update left
         updateActivities<ChangedBound::LOWER>(colview, lb[col], impliedlb,
                                               activities);

         lb[col] = impliedlb;
         changedCols.push_back(col);
      }
      else if ((impliedub < ub[col]) && impubfinite)
      {
         // update right
         updateActivities<ChangedBound::UPPER>(colview, ub[col], impliedub,
                                               activities);

         ub[col] = impliedub;
         changedCols.push_back(col);
      }
   }

   return true;
}

int
propagate(const MIP& mip, std::vector<double>& lb, std::vector<double>& ub,
          std::vector<Activity>& activities, int changedcol, double oldlb,
          double oldub)
{
   int ncols = mip.getNCols();

   std::vector<int> changedCols;
   changedCols.reserve(ncols);

   changedCols.push_back(changedcol);

   updateActivities(mip.getCol(changedcol), oldlb, lb[changedcol], oldub,
                    lb[changedcol], activities);

   for (size_t i = 0; i < changedCols.size(); ++i)
   {
      const int col = changedCols[i];

      auto [colcoefs, colindices, colsize] = mip.getCol(col);

      for (int j = 0; j < colsize; ++j)
      {
         const int row = colindices[j];

         // propagate row
         if (!propagateRow(mip, row, activities, lb, ub, changedCols))
            return 0;
      }
   }
   return 1;
}
