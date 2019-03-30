#include "Timer.h"

Timer::time_point
Timer::now()
{
   return tbb::tick_count::now();
}

double
Timer::seconds(time_point t1, time_point t0)
{
   return (t1 - t0).seconds();
}
