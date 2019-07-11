#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <limits>

struct Num
{
   static double round(double val);

   static double floor(double val);

   static double ceil(double val);

   static double greater(double lhs, double rhs);

   static double less(double lhs, double rhs);

   static double equal(double lhs, double rhs);

   static int sign(double val);

   static bool isIntegral(double val);

   static bool isInf(double val);

   static bool isMinusInf(double val);

   static constexpr double infinity();

   static constexpr double infval = std::numeric_limits<double>::infinity();

   private:
   static constexpr double constol = 1e-6;
};

#endif
