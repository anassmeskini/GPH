#ifndef INTSHIFTING_HPP
#define INTSHIFTING_HPP
#include "core/Heuristic.h"

class IntShifting : public HeuristicMethod
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

   virtual ~IntShifting() = default;
};

#endif
