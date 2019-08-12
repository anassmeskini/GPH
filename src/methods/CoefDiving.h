#ifndef COEF_DIVING_HPP
#define COEF_DIVING_HPP

#include "DivingHeuristic.h"
#include "core/Numerics.h"

#include <string>
#include <vector>

struct CoefDivingSelection
{
   static std::tuple<int, int, int>
   select(const MIP& mip, const std::vector<double>& solution)
   {
      const int ncols = mip.getNCols();
      const auto& integer = mip.getInteger();
      const auto& downLocks = mip.getDownLocks();
      const auto& upLocks = mip.getUpLocks();

      int varToFix = -1;
      int direction = 0;
      int minNLocks = std::numeric_limits<int>::max();
      int nFrac = 0;

      for (int col = 0; col < ncols; ++col)
      {
         if (integer[col] && !Num::isIntegral(solution[col]))
         {
            ++nFrac;
            if (std::min(downLocks[col], upLocks[col]) == 0)
               continue;

            if (downLocks[col] < minNLocks)
            {
               varToFix = col;
               minNLocks = downLocks[col];
               direction = -1;
            }

            if (upLocks[col] < minNLocks)
            {
               varToFix = col;
               minNLocks = upLocks[col];
               direction = 1;
            }
         }
      }

      return {varToFix, direction, nFrac};
   }

   static constexpr std::string_view name{"Coef"};
};

using CoefDiving = DivingHeuristic<CoefDivingSelection>;

#endif
