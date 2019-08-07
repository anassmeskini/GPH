#ifndef HEURISTIC_HPP
#define HEURISTIC_HPP
#include "core/Heuristic.h"

#include <vector>

class BoundSolution : public HeuristicMethod
{
   public:
   virtual void search(const MIP&,
                       const std::vector<double>&,
                       const std::vector<double>&,
                       const std::vector<Activity>&,
                       const LPResult&,
                       const std::vector<double>&,
                       const std::vector<int>&,
                       std::shared_ptr<const LPSolver>,
                       SolutionPool&) override;

   virtual ~BoundSolution() = default;
};

#endif
