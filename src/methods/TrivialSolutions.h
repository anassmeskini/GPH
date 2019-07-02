#ifndef HEURISTIC_HPP
#define HEURISTIC_HPP
#include "core/Heuristic.h"

#include <vector>

class TrivialSolutions : public HeuristicMethod
{
   public:
   virtual std::optional<std::vector<double>> search(const MIP<double>&,
                                                     const std::vector<double>&,
                                                     const std::vector<double>&,
                                                     const LPResult&,
                                                     const std::vector<double>&,
                                                     const std::vector<size_t>&,
                                                     const LPFactory&) override;

   virtual ~TrivialSolutions() = default;
};

#endif
