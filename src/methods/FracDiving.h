#ifndef FRAC_DIVING_HPP
#define FRAC_DIVING_HPP

#include "DivingHeuristic.h"
#include "core/Numerics.h"

#include <string>
#include <vector>

struct FracDivingSelection
{
   static std::tuple<int, int, int>
   select(const MIP& mip, const std::vector<double>& lb,
          const std::vector<double>& ub,
          const std::vector<double>& solution)
   {
      const int ncols = mip.getNCols();
      const auto& integer = mip.getInteger();
      const auto& downLocks = mip.getDownLocks();
      const auto& upLocks = mip.getUpLocks();

      int varToFix = -1;
      int direction = 0;
      double minFrac = std::numeric_limits<double>::max();
      int nFrac = 0;

      for (int col = 0; col < ncols; ++col)
      {
         if (Num::isFeasEQ(lb[col], ub[col]))
            continue;

         if (integer[col] && !Num::isIntegral(solution[col]))
         {
            ++nFrac;
            if (std::min(downLocks[col], upLocks[col]) == 0)
               continue;

            const double fractionality =
                solution[col] - Num::floor(solution[col]);

            if (fractionality < minFrac)
            {
               varToFix = col;
               minFrac = fractionality;
               direction = -1;
            }
            else if (1.0 - fractionality < minFrac)
            {
               varToFix = col;
               minFrac = 1.0 - fractionality;
               direction = 1;
            }
         }
      }

      return {varToFix, direction, nFrac};
   }

   static constexpr std::string_view name{"Frac"};
};

using FracDiving = DivingHeuristic<FracDivingSelection>;

#endif
