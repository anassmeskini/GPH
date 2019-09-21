#define CATCH_CONFIG_MAIN
#define UNIT_TEST

#include "catch2/catch.hpp"

#include "core/MIP.h"
#include "core/Propagation.h"
#include "make_mip.h"

TEST_CASE("Propagation test", "[Propagation]")
{
   // knapsack problem
   // 2x1 + 3x2 + 4x3 + x4 <= 3
   int ncols = 4;
   std::vector<double> coefs = {2.0, 3.0, 4.0, 1.0};
   std::vector<double> lhs = {-Num::infval};
   std::vector<double> rhs = {3.0};
   std::vector<double> lb = {0.0, 0.0, 0.0, 0.0};
   std::vector<double> ub = {1.0, 1.0, 1.0, 1.0};
   dynamic_bitset<> integer(4, true);
   std::vector<double> objective = {1.0, 1.0, 1.0, 1.0};
   MIP mip = make_mip(coefs, ncols, lhs, rhs, lb, ub, integer, objective);

   std::vector<Activity> activities = {{0.0, 10.0, 0, 0}};

   // fix x2 = 1 should result in fixing other
   // variables to 0
   mip.lb[1] = 1.0;
   REQUIRE(propagate(mip, mip.lb, mip.ub, activities, 1, 0.0, 1.0));

   REQUIRE(mip.ub[0] == 0.0);
   REQUIRE(mip.lb[1] == 1.0);
   REQUIRE(mip.ub[2] == 0.0);
   REQUIRE(mip.ub[3] == 0.0);

   // change the knapsack constraint tp
   // 2x1 - 3x2 + 4x3 + x4 <= 3
   // with x2 in {-1, 0}
   // and fix x2 to -1
   activities = {{0.0, 10.0, 0, 0}};
   mip.lb[1] = -1.0;
   mip.ub[1] = -1.0;
   change_coef(mip, {0, 1}, -3.0);

   REQUIRE(propagate(mip, mip.lb, mip.ub, activities, 1, -1.0, 0.0));

   REQUIRE(mip.ub[0] == 0.0);
   REQUIRE(mip.ub[1] == -1.0);
   REQUIRE(mip.ub[2] == 0.0);
   REQUIRE(mip.ub[3] == 0.0);
   REQUIRE(activities[0].min == 3.0);
   REQUIRE(activities[0].max == 10.0);
   REQUIRE(activities[0].ninfmin == 0);
   REQUIRE(activities[0].ninfmax == 0);

   // 2x1 + 3x2 + 4x3 + x4 <= 3
   // and fix x3 to 1
   change_coef(mip, {0, 1}, 3.0);
   mip.lb[1] = 0.0;
   mip.ub[1] = 1.0;
   mip.lb[2] = 1.0;

   activities = {{0.0, 10.0, 0, 0}};
   REQUIRE(!propagate(mip, mip.lb, mip.ub, activities, 2, 0.0, 1.0));
}
