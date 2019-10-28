#ifndef MIN_FRAC_HPP
#define MIN_FRAC_HPP
#include "core/Heuristic.h"

#include <vector>

class MinFracRounding : public HeuristicMethod
{
 public:
   MinFracRounding() : HeuristicMethod("FracRounding") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               TimeLimit, SolutionPool&) override;

   ~MinFracRounding() override = default;
};

#endif
