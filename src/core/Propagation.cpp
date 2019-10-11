#include "Propagation.h"
#include "Common.h"
#include "Numerics.h"
#include "io/Message.h"

template <>
bool
updateActivities<ChangedBound::LOWER>(
    VectorView colview, double oldlb, double newlb,
    std::vector<Activity>& activities, const std::vector<double>& lhs,
    const std::vector<double>& rhs) noexcept
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

      if ((activities[row].ninfmin == 0 &&
           !Num::isFeasLE(activities[row].min, rhs[row])) ||
          (activities[row].ninfmax == 0 &&
           !Num::isFeasGE(activities[row].max, lhs[row])))
         return false;
   }

   return true;
}

template <>
bool
updateActivities<ChangedBound::UPPER>(
    VectorView colview, double oldub, double newub,
    std::vector<Activity>& activities, const std::vector<double>& lhs,
    const std::vector<double>& rhs) noexcept
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

      if ((activities[row].ninfmin == 0 &&
           !Num::isFeasLE(activities[row].min, rhs[row])) ||
          (activities[row].ninfmax == 0 &&
           !Num::isFeasGE(activities[row].max, lhs[row])))
         return false;
   }
   return true;
}

bool
updateActivities(VectorView colview, double oldlb, double newlb,
                 double oldub, double newub,
                 std::vector<Activity>& activities,
                 const std::vector<double>& lhs,
                 const std::vector<double>& rhs) noexcept
{
   auto [colcoefs, colindices, colsize] = colview;

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

         if ((activities[row].ninfmin == 0 &&
              !Num::isFeasLE(activities[row].min, rhs[row])) ||
             (activities[row].ninfmax == 0 &&
              !Num::isFeasGE(activities[row].max, lhs[row])))
            return false;
      }
   }
   return true;
}

static bool
propagateRow(const MIP& problem, int row,
             std::vector<Activity>& activities, std::vector<double>& lb,
             std::vector<double>& ub,
             std::vector<int>& changedCols) noexcept
{
   auto [rowcoefs, rowindices, rowsize] = problem.getRow(row);

   auto& activity = activities[row];

   auto st = problem.getStats();

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

         if (col < st.nbin + st.nint)
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

         if (col < st.nbin + st.nint)
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

      if ((impliedlb > lb[col] + 1e-6 && implbfinite) &&
          (impliedub < ub[col] - 1e-6 && impubfinite))
      {
         auto colview = problem.getCol(col);
         // update right and left
         if (!updateActivities(colview, lb[col], impliedlb, ub[col],
                               impliedub, activities, lhs, rhs))
            return false;

         lb[col] = impliedlb;
         ub[col] = impliedub;
         changedCols.push_back(col);
      }
      else if (impliedlb > lb[col] + 1e-6 && implbfinite)
      {
         auto colview = problem.getCol(col);
         // update right and left
         if (!updateActivities<ChangedBound::LOWER>(
                 colview, lb[col], impliedlb, activities, lhs, rhs))
            return false;

         lb[col] = impliedlb;
         changedCols.push_back(col);
      }
      else if (impliedub < ub[col] - 1e-6 && impubfinite)
      {
         auto colview = problem.getCol(col);
         // update right and left
         if (!updateActivities<ChangedBound::UPPER>(
                 colview, ub[col], impliedub, activities, lhs, rhs))
            return false;

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

   if (!updateActivities(mip.getCol(changedcol), oldlb, lb[changedcol],
                         oldub, ub[changedcol], activities, mip.getLHS(),
                         mip.getRHS()))
      return false;

   for (size_t i = 0; i < changedCols.size(); ++i)
   {
      const int col = changedCols[i];

      auto [colcoefs, rowindices, colsize] = mip.getCol(col);

      for (int j = 0; j < colsize; ++j)
      {
         const int row = rowindices[j];

         // propagate row
         if (!propagateRow(mip, row, activities, lb, ub, changedCols))
            return 0;
      }
   }

   return 1;
}

bool
propagate_get_changed_cols(const MIP& mip, std::vector<double>& lb,
                           std::vector<double>& ub,
                           std::vector<Activity>& activities,
                           int changedcol, double oldlb, double oldub,
                           std::vector<int>& changedCols)
{
   assert(changedCols.empty());
   int ncols = mip.getNCols();

   changedCols.reserve(ncols);

   changedCols.push_back(changedcol);

   if (!updateActivities(mip.getCol(changedcol), oldlb, lb[changedcol],
                         oldub, ub[changedcol], activities, mip.getLHS(),
                         mip.getRHS()))
      return false;

   for (size_t i = 0; i < changedCols.size(); ++i)
   {
      const int col = changedCols[i];

      auto [colcoefs, colindices, colsize] = mip.getCol(col);

      for (int j = 0; j < colsize; ++j)
      {
         const int row = colindices[j];

         // propagate row
         if (!propagateRow(mip, row, activities, lb, ub, changedCols))
            return false;
      }
   }
   return true;
}
