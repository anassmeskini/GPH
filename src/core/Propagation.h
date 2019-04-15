#ifndef PROPAGATION_HPP
#define PROPAGATION_HPP

#include "Problem.h"

// assumes that the lb and ub reftlect the fixing, and the activities are up to
// date
int
propagate(const ProblemView& problem,
          std::vector<double>& lb,
          std::vector<double>& ub,
          std::vector<Activity>& activities,
          size_t col);

void
updateActivities(VectorView<double> colview,
                 double oldlb,
                 double newlb,
                 double oldub,
                 double newub,
                 std::vector<Activity>& activities);

#endif
