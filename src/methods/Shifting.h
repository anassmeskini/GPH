#ifndef SHIFTING_HPP
#define SHIFTING_HPP
#include "core/Heuristic.h"

#include <vector>

class Shifting : public FeasibilityHeuristic
{
 public:
   Shifting() : FeasibilityHeuristic("Shifting") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               TimeLimit, SolutionPool&) override;

   ~Shifting() override = default;
};

#endif
