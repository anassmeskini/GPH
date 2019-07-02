#ifndef PROPAGATION_HPP
#define PROPAGATION_HPP

#include "Problem.h"

// assumes that the lb and ub reftlect the fixing, and the activities are up to
// date
int
propagate(const MIP<double>& problem,
          std::vector<double>& lb,
          std::vector<double>& ub,
          std::vector<Activity>& activities,
          size_t col);

bool
updateActivities(VectorView<double> colview,
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
updateActivities(VectorView<double> colview,
                 double oldb,
                 double newb,
                 std::vector<Activity>& activities,
                 const std::vector<double>& lhs,
                 const std::vector<double>& rhs);

#endif
