#ifndef PROPAGATION_HPP
#define PROPAGATION_HPP

#include "Common.h"
#include "MIP.h"

int
propagate(const MIP &problem, std::vector<double> &lb,
          std::vector<double> &ub, std::vector<Activity> &activities,
          int col, double oldlb, double oldub);

void
updateActivities(VectorView colview, double oldlb, double newlb,
                 double oldub, double newub,
                 std::vector<Activity> &activities) noexcept;

enum class ChangedBound
{
   LOWER,
   UPPER,
};

template <ChangedBound chgbd>
void
updateActivities(VectorView colview, double oldb, double newb,
                 std::vector<Activity> &activities) noexcept;
#endif
