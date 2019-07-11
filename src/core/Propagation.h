#ifndef PROPAGATION_HPP
#define PROPAGATION_HPP

#include "Common.h"
#include "MIP.h"

// assumes that the lb and ub reftlect the fixing, and the activities are up to
// date
int
propagate(const MIP& problem,
          std::vector<double>& lb,
          std::vector<double>& ub,
          std::vector<Activity>& activities,
          int col);

bool
updateActivities(VectorView colview,
                 double oldlb,
                 double newlb,
                 double oldub,
                 double newub,
                 std::vector<Activity>& activities,
                 const std::vector<double>& lhs,
                 const std::vector<double>& rhs);

enum class ChangedBound
{
   LOWER,
   UPPER,
};

template<ChangedBound chgbd>
bool
updateActivities(VectorView colview,
                 double oldb,
                 double newb,
                 std::vector<Activity>& activities,
                 const std::vector<double>& lhs,
                 const std::vector<double>& rhs);

#endif
