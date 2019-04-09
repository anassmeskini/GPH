#ifndef TIMER_HPP
#define TIMER_HPP

#include <tbb/tick_count.h>

struct Timer
{
  using time_point = tbb::tick_count;

  static time_point now();

  static double seconds(time_point, time_point);
};

#endif
