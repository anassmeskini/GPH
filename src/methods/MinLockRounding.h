#ifndef MINLOCK2_HPP
#define MINLOCK2_HPP
#include "core/Heuristic.h"

#include <vector>

class MinLockRounding : public HeuristicMethod
{
 public:
   MinLockRounding() : HeuristicMethod("Rounding") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               SolutionPool&) override;

   ~MinLockRounding() override = default;
};

#endif
