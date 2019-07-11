#ifndef MINLOCK_HPP
#define MINLOCK_HPP
#include "core/Heuristic.h"

#include <vector>

class MinLockRounding : public HeuristicMethod
{
   public:
   virtual std::optional<Solution> search(const MIP&,
                                          const std::vector<double>&,
                                          const std::vector<double>&,
                                          const std::vector<Activity>&,
                                          const LPResult&,
                                          const std::vector<double>&,
                                          const std::vector<int>&,
                                          std::shared_ptr<LPSolver>) override;

   virtual ~MinLockRounding() = default;
};

#endif
