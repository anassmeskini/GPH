#ifndef TIMER_HPP
#define TIMER_HPP

#include <tbb/tick_count.h>
#include <tbb/atomic.h>

struct Timer
{
   using time_point = tbb::tick_count;

   static time_point now() {return tbb::tick_count::now();}

   static double seconds(time_point t1, time_point t0) {return (t1 - t0).seconds();}
};

class TimeLimit
{
 public:
   TimeLimit(Timer::time_point t, int l) : start(t), limit(l) {}

   bool reached(Timer::time_point t) const
   {
      return static_cast<int>(Timer::seconds(t, start)) > limit;
   }

 private:
   Timer::time_point start;
   int limit;
};

#endif
