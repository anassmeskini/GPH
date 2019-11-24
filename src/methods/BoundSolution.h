#ifndef HEURISTIC_HPP
#define HEURISTIC_HPP
#include "core/Heuristic.h"

#include <vector>

class BoundSolution final : public FeasibilityHeuristic
{
 public:
   BoundSolution() : FeasibilityHeuristic("BoundSolution") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               TimeLimit, SolutionPool&) override;

   ~BoundSolution() override = default;

 private:
   bool tryUBSolution(const MIP&, std::vector<double>& lb,
                      std::vector<double>& ub,
                      const std::vector<Activity>& activities,
                      TimeLimit) const;

   bool tryLBSolution(const MIP&, std::vector<double>& lb,
                      std::vector<double>& ub,
                      const std::vector<Activity>& activities,
                      TimeLimit) const;

   bool tryOptimisticSolution(const MIP&, std::vector<double>& lb,
                              std::vector<double>& ub,
                              const std::vector<Activity>& activities,
                              TimeLimit) const;
};

#endif
