#include "Common.h"
#include "SparseMatrix.h"

Time::time_point
Time::now()
{
   return std::chrono::high_resolution_clock::now();
}

double
Time::seconds(time_point t1, time_point t0)
{
   using namespace std::chrono;
   return duration_cast<duration<double>>(t1 - t0).count();
}
