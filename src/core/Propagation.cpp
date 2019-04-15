#include "Propagation.h"
#include "Common.h"
#include "Numerics.h"
#include "io/Message.h"

#include <memory>

// TODO remove the need for a call back, update the activities inside

// the callback is called whenever a bound thightening is possible on the row
// It should be used to update the activities of the rows in the support of
// the modified columns, the callback returns whether it found a conflict
// i.e. max{Act_i} < lhs_i or min{Act_i} > rhs_i
template<typename CALLBLACK>
static bool
propagateRow(VectorView<double> rowview,
             double lhs,
             double rhs,
             Activity activity,
             const std::vector<double>& lb,
             const std::vector<double>& ub,
             const bitset& integer,
             CALLBLACK&& callback)
{
   int nbndchg = 0;

   const double* rowcoefs = rowview.coefs;
   const size_t* rowindices = rowview.indices;
   const size_t rowsize = rowview.size;

   for (size_t i = 0; i < rowsize; ++i)
   {
      size_t col = rowindices[i];
      double coef = rowcoefs[i];

      double impliedlb;
      double impliedub;
      double maxPartialActivity = activity.max;
      double minPartialActivity = activity.min;

      if (Num::greater(coef, 0.0))
      {
         maxPartialActivity -= coef * ub[col];
         minPartialActivity -= coef * lb[col];

         if (integer[col])
         {
            impliedlb = Num::ceil((lhs - maxPartialActivity) / coef);
            impliedub = Num::floor((rhs - minPartialActivity) / coef);
         }
         else
         {
            impliedlb = (lhs - maxPartialActivity) / coef;
            impliedub = (rhs - minPartialActivity) / coef;
         }
      }
      else
      {
         assert(coef != 0.0);
         maxPartialActivity -= coef * lb[col];
         minPartialActivity -= coef * ub[col];

         if (integer[col])
         {
            impliedlb = Num::ceil((rhs - minPartialActivity) / coef);
            impliedub = Num::floor((lhs - maxPartialActivity) / coef);
         }
         else
         {
            impliedlb = (rhs - minPartialActivity) / coef;
            impliedub = (lhs - maxPartialActivity) / coef;
         }

         // TODO
      }

      if (!Num::less(activity.min, rhs) || !Num::greater(activity.max, lhs))
      {
         Message::debug("Propagation detected infesiblity 1");
         return false;
      }

      bool updateActivities = false;

      if (Num::greater(impliedlb, lb[col]))
      {

         ++nbndchg;
         updateActivities = true;
         if (Num::greater(coef, 0.0))
            activity.min += coef * (impliedlb - lb[col]);
         else
            activity.max += coef * (impliedlb - lb[col]);
      }

      if (Num::greater(impliedub, ub[col]))
      {
         ++nbndchg;
         updateActivities = true;
         if (Num::greater(coef, 0.0))
            activity.max += coef * (impliedub - ub[col]);
         else
            activity.min += coef * (impliedub - ub[col]);
      }

      if (updateActivities)
      {
         if (callback(col, impliedlb, impliedub))
         {
            Message::debug("Propagation detected infesiblity 2");
            return false;
         }
      }
   }

   Message::debug("Propagation changed {} bounds", nbndchg);

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

   const auto& lhs = problem.getLHS();
   const auto& rhs = problem.getRHS();
   const auto& integer = problem.integer();
   size_t nrows = problem.getNRows();
   size_t ncols = problem.getNCols();

   std::vector<size_t> changedCols;
   bitset colChanged(ncols, false);
   changedCols.reserve(nrows);

   changedCols.push_back(changedcol);
   colChanged[changedcol] = true;

   const auto propagation_callback =
     [&](size_t col, double newlb, double newub) -> bool {
      assert(newlb != lb[col] || newub != ub[col]);

      // add it to the candidate list of propagation if not already there
      if (!colChanged[col])
      {
         changedCols.push_back(col);
         colChanged[col] = true;
      }

      auto colview = problem.getCol(col);
      auto colindices = colview.indices;
      auto colcoefs = colview.coefs;
      size_t colsize = colview.size;

      Message::debug("updating activities {}", colsize);

      for (size_t i = 0; i < colsize; ++i)
      {
         const size_t row = colindices[i];
         const double coef = colcoefs[i];

         if (Num::greater(coef, 0.0))
         {
            activities[row].min += coef * (newlb - lb[col]);
            activities[row].max += coef * (newub - ub[col]);
         }
         else
         {
            activities[row].min += coef * (newub - ub[col]);
            activities[row].max += coef * (newlb - lb[col]);
         }

         // stop and return
         if (!Num::less(activities[row].min, rhs[row]) ||
             !Num::greater(activities[row].max, lhs[row]))
            return true;
      }

      return false;
   };

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
         if (!propagateRow(problem.getRow(row),
                           lhs[row],
                           rhs[row],
                           activities[row],
                           lb,
                           ub,
                           integer,
                           propagation_callback))
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

enum class ChangedBound : uint8_t
{
   LOWER,
   UPPER,
};

template<ChangedBound chgbd>
void
updateActivities(VectorView<double> colview,
                 double oldb,
                 double newb,
                 std::vector<Activity>& activities);

template<>
void
updateActivities<ChangedBound::LOWER>(VectorView<double> colview,
                                      double oldb,
                                      double newb,
                                      std::vector<Activity>& activities)
{
   const double* colcoefs = colview.coefs;
   const size_t* colindices = colview.indices;
   const size_t colsize = colview.size;

   for (size_t i = 0; i < colsize; ++i)
   {
      size_t row = colindices[i];
      const double coef = colcoefs[i];

      if (Num::greater(coef, 0.0))
         activities[row].min += coef * (newb - oldb);
      else
         activities[row].max += coef * (newb - oldb);
   }
}

template<>
void
updateActivities<ChangedBound::UPPER>(VectorView<double> colview,
                                      double oldb,
                                      double newb,
                                      std::vector<Activity>& activities)
{
   const double* colcoefs = colview.coefs;
   const size_t* colindices = colview.indices;
   const size_t colsize = colview.size;

   for (size_t i = 0; i < colsize; ++i)
   {
      size_t row = colindices[i];
      const double coef = colcoefs[i];

      if (Num::greater(coef, 0.0))
         activities[row].max += coef * (newb - oldb);
      else
         activities[row].min += coef * (newb - oldb);
   }
}

// TODO make a specialization for changes only in one bound
void
updateActivities(VectorView<double> colview,
                 double oldlb,
                 double newlb,
                 double oldub,
                 double newub,
                 std::vector<Activity>& activities)
{
   const double* colcoefs = colview.coefs;
   const size_t* colindices = colview.indices;
   const size_t colsize = colview.size;

   for (size_t i = 0; i < colsize; ++i)
   {
      size_t row = colindices[i];
      const double coef = colcoefs[i];

      if (Num::greater(coef, 0.0))
      {
         activities[row].min += coef * (newlb - oldlb);
         activities[row].max += coef * (newub - oldub);
      }
      else
      {
         activities[row].min += coef * (newub - oldub);
         activities[row].max += coef * (newlb - oldlb);
      }
   }
}
