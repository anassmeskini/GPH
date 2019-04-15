#ifndef HEURISTIC_HPP
#define HEURISTIC_HPP
#include "core/Heuristic.h"

#include <vector>

class TrivialSolutions : public HeuristicMethod
{
   public:
   virtual std::optional<std::vector<double>> search(
     const ProblemView& problem) override;

   virtual ~TrivialSolutions() = default;
};

#endif
