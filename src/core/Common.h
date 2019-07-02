#ifndef _COMMON_HPP_
#define _COMMON_HPP_

#include "Bitset.h"
#include "MIP.h"
#include "Numerics.h"
#include "SparseMatrix.h"

#include <cassert>
#include <chrono>
#include <memory>

template<typename REAL>
struct VectorView
{
   VectorView(const REAL* _array, const size_t* _indices, size_t _size)
     : coefs(_array)
     , indices(_indices)
     , size(_size)
   {
   }

   VectorView(const VectorView<REAL>&) = default;

   const REAL* coefs;
   // use int?
   const size_t* indices;
   size_t size;
};

struct Activity
{
   double min = 0.0;
   double max = 0.0;

   int ninfmin = 0;
   int ninfmax = 0;
};

std::vector<Activity>
computeActivities(const MIP<double>& mip);

std::vector<double>
computeSolActivities(const MIP<double>& mip, const std::vector<double>& sol);

template<typename REAL>
bool
checkFeasibility(const MIP<REAL>& mip,
                 const std::vector<REAL>& sol,
                 REAL boundtol,
                 REAL constol)
{
   const auto& ub = mip.getUB();
   const auto& lb = mip.getLB();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();
   const auto& integer = mip.getInteger();

   assert(sol.size() == mip.getNCols());

   for (size_t colid = 0; colid < mip.getNCols(); ++colid)
   {
      if (sol[colid] > ub[colid] + boundtol ||
          sol[colid] < lb[colid] - boundtol)
         return false;

      if (integer[colid] && !Num::isIntegral(sol[colid]))
         return false;
   }

   for (size_t rowid = 0; rowid < mip.getNRows(); ++rowid)
   {
      auto& row = mip.getRow(rowid);
      REAL activity{ 0.0 };

      for (size_t id = 0; id < row.size; ++id)
         activity += row.coefs[id] * sol[row.indices[id]];

      if (activity > rhs[rowid] + constol || activity < lhs[rowid] - constol)
         return false;
   }

   return true;
}

#endif
