#ifndef TRIVIAL_ROUND_HPP
#define TRIVIAL_ROUND_HPP
#include "core/Heuristic.h"

#include <vector>

class TrivialRounding : public FeasibilityHeuristic
{
 public:
   TrivialRounding() : FeasibilityHeuristic("TrivialRounding") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               TimeLimit, SolutionPool&) override;

   ~TrivialRounding() override = default;

   void setParam(const std::string&,
                 const std::variant<std::string, int, double>&) override
   {
   }
};

#endif
