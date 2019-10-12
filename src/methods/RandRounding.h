#ifndef RAND_ROUND_HPP
#define RAND_ROUND_HPP
#include "core/Heuristic.h"

#include <vector>

class RandRounding : public HeuristicMethod
{
 public:
   RandRounding() : HeuristicMethod("RandRounding") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               SolutionPool&) override;

   ~RandRounding() override = default;
};

#endif
