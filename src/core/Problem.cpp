#include "Problem.h"
#include "Numerics.h"

ProblemView::ProblemView(MIP<double>&& _mip, std::vector<double>&& _lpSol)
  : mip(std::move(_mip))
  , ncols(mip.getNCols())
  , nrows(mip.getNRows())
  , activities(nrows)
  , lpSol(_lpSol)
  , lpSolActivity(nrows)
  , valuesType(ncols)
  , downLocks(ncols, 0)
  , upLocks(ncols, 0)
{
   const auto& integer = mip.getInteger();
   const auto& lb = mip.getLB();
   const auto& ub = mip.getUB();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();

   for (size_t i = 0; i < ncols; ++i)
   {
      if (integer[i])
      {
         if (Num::isIntegral(lpSol[i]))
            valuesType[i] = ValueType::INTEGER_VALUE;
         else
         {
            fractional.push_back(i);

            if (lb[i] == 0.0 && ub[i] == 1.0)
               valuesType[i] = ValueType::BINARY_FRACTIONAL;
            else
               valuesType[i] = ValueType::INTEGER_FRACTIONAL;
         }
      }
      else
         valuesType[i] = ValueType::CONTINUOUS_VARIABLE;
   }

   // compute activites and locks
   for (size_t row = 0; row < nrows; ++row)
   {
      double solActivity = 0.0;
      double minActivity = 0.0;
      double maxActivity = 0.0;

      // TODO use fold expressions
      auto rowview = mip.getRow(row);
      size_t rowsize = rowview.size;
      const double* coefs = rowview.coefs;
      const size_t* indices = rowview.indices;

      bool lhsfinite = !Num::isMinusInf(lhs[row]);
      bool rhsfinite = !Num::isInf(rhs[row]);

      for (size_t id = 0; id < rowsize; ++id)
      {
         size_t col = indices[id];
         double coef = coefs[id];
         assert(coef != 0.0);

         solActivity += coef * lpSol[col];

         if (Num::greater(coef, 0.0))
         {
            minActivity += coef * lb[col];
            maxActivity += coef * ub[col];

            if (lhsfinite)
               ++downLocks[col];
            if (rhsfinite)
               ++upLocks[col];
         }
         else
         {
            minActivity += coef * ub[col];
            maxActivity += coef * lb[col];

            if (lhsfinite)
               ++upLocks[col];
            if (rhsfinite)
               ++downLocks[col];
         }
      }

      lpSolActivity[row] = solActivity;
      activities[row].min = minActivity;
      activities[row].max = maxActivity;
   }
}

size_t
ProblemView::getNCols() const
{
   return ncols;
}

size_t
ProblemView::getNRows() const
{
   return nrows;
}

VectorView<double>
ProblemView::getRow(size_t row) const
{
   return mip.getRow(row);
}

VectorView<double>
ProblemView::getCol(size_t col) const
{
   return mip.getCol(col);
}

const std::vector<size_t>&
ProblemView::getUpLocks() const
{
   return upLocks;
}

const std::vector<size_t>&
ProblemView::getDownLocks() const
{
   return downLocks;
}

const bitset&
ProblemView::integer() const
{
   return mip.getInteger();
}

const std::vector<Activity>&
ProblemView::getActivities() const
{
   return activities;
}

const std::vector<double>&
ProblemView::getLPSol() const
{
   return lpSol;
}

ValueType
ProblemView::getValueType(size_t col) const
{
   return valuesType[col];
}

const std::vector<double>&
ProblemView::getSolActivity() const
{
   return lpSolActivity;
}

const std::vector<double>&
ProblemView::getLB() const
{
   return mip.getLB();
}
const std::vector<double>&
ProblemView::getUB() const
{
   return mip.getUB();
}

const std::vector<double>&
ProblemView::getLHS() const
{
   return mip.getLHS();
}

const std::vector<double>&
ProblemView::getRHS() const
{
   return mip.getRHS();
}
