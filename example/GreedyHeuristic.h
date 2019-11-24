#include "core/Heuristic.h"

class GreedyHeuristic final : public FeasibilityHeuristic
{
 public:
   GreedyHeuristic() : FeasibilityHeuristic("SetCoverGreedy") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               TimeLimit, SolutionPool&) override;

 private:
   static constexpr int itermax = 1000;
   static constexpr int cand_list_max = 10;
   static constexpr double priority = 0.5;
   static constexpr double improvement = 0.2;
};
