#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <cmath>
#include <limits>

struct Num
{
   static double round(double val) { return std::floor(val + 0.5); }

   static double floor(double val) { return std::floor(val); }

   static double ceil(double val) { return std::ceil(val); }

   static bool isGE(double lhs, double rhs)
   {
      return lhs - rhs >= -epsilon;
   }

   static bool isLE(double lhs, double rhs)
   {
      return lhs - rhs <= epsilon;
   }

   static bool isEQ(double lhs, double rhs)
   {
      return std::fabs(rhs - lhs) <= epsilon;
   }

   static bool isFeasGE(double lhs, double rhs)
   {
      return lhs - rhs >= -constol;
   }

   static bool isFeasLE(double lhs, double rhs)
   {
      return lhs - rhs <= constol;
   }

   static bool isFeasEQ(double lhs, double rhs)
   {
      return std::fabs(rhs - lhs) <= constol;
   }

   static bool isFeasInt(double val)
   {
      return std::fabs(val - round(val)) < epsilon;
   }

   static int sign(double val)
   {
      if (val > 0.0)
         return 1;
      return -1;
   }

   static bool isIntegral(double val)
   {
      return std::floor(val) == std::ceil(val);
   }

   static bool isInf(double val) { return val == infval; }

   static bool isMinusInf(double val) { return val == -infval; }

   static constexpr double infinity() { return infval; }

   static constexpr double infval =
       std::numeric_limits<double>::infinity();

 private:
   static constexpr double constol = 1e-6;

   static constexpr double epsilon = 1e-9;
};

#endif
