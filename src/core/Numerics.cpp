#include "Numerics.h"
#include <cmath>

bool
Num::isIntegral(double val)
{
   return std::floor(val) == std::ceil(val);
}

double
Num::round(double val)
{
   return std::floor(val + 0.5);
}

double
Num::floor(double val)
{
   return std::floor(val);
}

double
Num::ceil(double val)
{
   return std::ceil(val);
}

double
Num::greater(double lhs, double rhs)
{
   return lhs >= rhs + constol;
}

double
Num::less(double lhs, double rhs)
{
   return lhs <= rhs - constol;
}

double
Num::equal(double lhs, double rhs)
{
   return lhs <= rhs + constol && lhs >= rhs - constol;
}

int
Num::sign(double val)
{
   if (val > 0.0)
      return 1;
   return -1;
}

bool
Num::isInf(double val)
{
   return val == infval;
}

bool
Num::isMinusInf(double val)
{
   return val == -infval;
}

constexpr double
Num::infinity()
{
   return infval;
}
