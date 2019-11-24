#ifndef RAND_ROUND_HPP
#define RAND_ROUND_HPP
#include "core/Heuristic.h"

#include <vector>

class RandRounding : public FeasibilityHeuristic
{
 public:
   RandRounding() : FeasibilityHeuristic("RandRounding") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               TimeLimit, SolutionPool&) override;

   ~RandRounding() override = default;

   void
   setParam(const std::string& param,
            const std::variant<std::string, int, double>& value) override
   {
      if (param == "itermax")
      {
         // TODO make sure it is nonegative
         itermax = std::get<int>(value);
      }
   }

   int itermax = 10;
};

#endif
