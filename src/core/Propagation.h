#ifndef PROPAGATION_HPP
#define PROPAGATION_HPP

#include "Common.h"
#include "MIP.h"

int
propagate(const MIP& problem, std::vector<double>& lb,
          std::vector<double>& ub, std::vector<Activity>& activities,
          int col, double oldlb, double oldub);

bool
updateActivities(VectorView colview, double oldlb, double newlb,
                 double oldub, double newub,
                 std::vector<Activity>& activities,
                 const std::vector<double>& lhs,
                 const std::vector<double>& rhs) noexcept;

enum class ChangedBound
{
   LOWER,
   UPPER,
};

template <ChangedBound chgbd>
bool
updateActivities(VectorView colview, double oldb, double newb,
                 std::vector<Activity>& activities,
                 const std::vector<double>& lhs,
                 const std::vector<double>& rhs) noexcept;

bool
propagate_get_changed_cols(const MIP& mip, std::vector<double>& lb,
                           std::vector<double>& ub,
                           std::vector<Activity>& activities,
                           int changedcol, double oldlb, double oldub,
                           std::vector<int>& buffer);
#endif
