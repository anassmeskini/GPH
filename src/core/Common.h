#ifndef _COMMON_HPP_
#define _COMMON_HPP_

#include "MIP.h"
#include "Numerics.h"
#include "SparseMatrix.h"

#include <cassert>
#include <chrono>
#include <memory>

struct Activity
{
   double min = 0.0;
   double max = 0.0;

   int ninfmin = 0;
   int ninfmax = 0;
};

struct LPSolInfo
{
   std::vector<double> rowActivities;
   std::vector<int> fractional;
};

std::vector<Activity>
computeActivities(const MIP& mip);

std::vector<double>
computeSolActivities(const MIP& mip, const std::vector<double>& sol);

int
updateSolActivity(std::vector<double>&, VectorView,
                  const std::vector<double>&, const std::vector<double>&,
                  double, std::vector<int>&, dynamic_bitset<>&);

std::vector<int>
getFractional(const std::vector<double>&, const dynamic_bitset<>&);

std::vector<int>
getIndentity(int ncols);

template <typename REAL, bool LP = false>
bool
checkFeasibility(const MIP& mip, const std::vector<double>& sol,
                 REAL boundtol = 1e-9, REAL constol = 1e-6)
{
   const auto& ub = mip.getUB();
   const auto& lb = mip.getLB();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();
   const auto& integer = mip.getInteger();

   assert(sol.size() == static_cast<size_t>(mip.getNCols()));

   for (int colid = 0; colid < mip.getNCols(); ++colid)
   {
      if (sol[colid] > ub[colid] + boundtol ||
          sol[colid] < lb[colid] - boundtol)
         return false;

      if constexpr (!LP)
      {
         if (integer[colid] && !Num::isFeasInt(sol[colid]))
            return false;
      }
   }

   for (int rowid = 0; rowid < mip.getNRows(); ++rowid)
   {
      auto [coefs, indices, size] = mip.getRow(rowid);
      REAL activity{0.0};

      for (int id = 0; id < size; ++id)
         activity += coefs[id] * sol[indices[id]];

      if (activity > rhs[rowid] + constol ||
          activity < lhs[rowid] - constol)
         return false;
   }

   return true;
}

template <typename REAL, bool LP = false>
int
getNViolated(const MIP& mip, const std::vector<double>& sol,
             REAL boundtol = 1e-9, REAL constol = 1e-6)
{
   const auto& ub = mip.getUB();
   const auto& lb = mip.getLB();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();
   const auto& integer = mip.getInteger();

   assert(sol.size() == static_cast<size_t>(mip.getNCols()));

   int ninfeas = 0;
   for (int colid = 0; colid < mip.getNCols(); ++colid)
   {
      if (sol[colid] > ub[colid] + boundtol ||
          sol[colid] < lb[colid] - boundtol)
         ++ninfeas;

      if constexpr (!LP)
      {
         if (integer[colid] && !Num::isFeasInt(sol[colid]))
            ++ninfeas;
      }
   }

   for (int rowid = 0; rowid < mip.getNRows(); ++rowid)
   {
      auto [coefs, indices, size] = mip.getRow(rowid);
      REAL activity{0.0};

      for (int id = 0; id < size; ++id)
         activity += coefs[id] * sol[indices[id]];

      if (activity > rhs[rowid] + constol ||
          activity < lhs[rowid] - constol)
         ++ninfeas;
   }

   return ninfeas;
}

template <typename COMP>
void
sortRows(SparseMatrix& mat, COMP&& comp)
{
   const std::vector<int> identity = getIndentity(mat.ncols);
   std::vector<int> permutation;
   std::vector<int> indcopy;
   std::vector<double> coefcopy;

   permutation.reserve(mat.ncols);
   indcopy.reserve(mat.ncols);
   coefcopy.reserve(mat.ncols);

   for (int row = 0; row < mat.nrows; ++row)
   {
      const int rowStart = mat.rowStart[row];
      const int size = mat.rowStart[row + 1] - mat.rowStart[row];
      int* indices = mat.indices.data() + rowStart;
      double* coefs = mat.coefficients.data() + rowStart;

      assert(size < mat.ncols);

      permutation.resize(size);
      indcopy.resize(size);
      coefcopy.resize(size);

      std::memcpy(permutation.data(), identity.data(), sizeof(int) * size);

      std::sort(std::begin(permutation), std::end(permutation),
                [&](int left, int right) {
                   return comp(indices[left], indices[right]);
                });

      for (size_t i = 0; i < permutation.size(); ++i)
      {
         indcopy[i] = indices[permutation[i]];
         coefcopy[i] = coefs[permutation[i]];
      }

      std::memcpy(indices, indcopy.data(), sizeof(int) * size);
      std::memcpy(coefs, coefcopy.data(), sizeof(double) * size);
   }
}

bool
hasZeroLockRounding(const std::vector<int>&,  // down locks
                    const std::vector<int>&,  // up locks
                    const std::vector<int>&); // fractional columns

bool
hasZeroLockRounding(const std::vector<double>&, // lp solution
                    const std::vector<int>&,    // down locks
                    const std::vector<int>&,    // up locks
                    const dynamic_bitset<>&);   // is integer

// assumes the solution has a zero lock rounding
// returns the diff in the cost
double
zeroLockRound(std::vector<double>&,        // lp solution
              const std::vector<int>&,     // down locks
              const std::vector<int>&,     // fractional
              const std::vector<double>&); // objective

std::optional<std::pair<std::vector<double>, double>>
minLockRound(const MIP& mip,             // mip
             const std::vector<double>&, // solution
             double,                     // objective
             const std::vector<int>&);   // fractional columns

void
maxOutSolution(const MIP& mip, std::vector<double>& solution,
               const std::vector<double>& activity);
void
roundFeasIntegers(std::vector<double>& sol,
                  const dynamic_bitset<>& integer);

template <typename T, typename PREDICATE>
bool all_of(const std::vector<T>& v1, const std::vector<T>& v2, PREDICATE&& pred)
{
   assert(v1.size() == v2.size());

   for(size_t i = 0; i < v1.size(); ++i)
   {
      if(!pred(v1[i], v2[i]))
         return false;
   }
   return true;
}
#endif
