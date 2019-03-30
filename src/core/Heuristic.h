#ifndef _HEURISTIC_HPP
#define _HEURISTIC_HPP

#include "AvaiLPSolver.h"
#include "MIP.h"

#include <memory>
#include <vector>

struct OptionalSolution
{
   std::unique_ptr<double[]> values = nullptr;
   size_t size = 0;

   OptionalSolution(OptionalSolution&& other)
   {
      values = std::move(other.values);
      size = other.size;
   }
};

// TODO templatize
// TODO use flags for inf
struct Activity
{
   double max;
   double min;
   double current;
};

// TODO templatize
class HeuristicMethod
{
 public:
   virtual OptionalSolution solve(const MIP<double>&) = 0;

   virtual ~HeuristicMethod() {}
};

// TODO templatize
class Heuristics
{
 public:
   OptionalSolution solve(const MIP<double>&);

 private:
   std::vector<std::unique_ptr<HeuristicMethod>> heuristics;

   std::vector<size_t> upLocks;
   std::vector<size_t> dowsLocks;

   std::vector<Activity> activities;
};

OptionalSolution
Heuristics::solve(const MIP<double>& mip)
{
   std::vector<OptionalSolution> solutions;
   for (size_t id = 0; id < heuristics.size(); ++id)
      solutions.push_back(heuristics[id]->solve(mip));

   return std::move(solutions[0]);
}
#endif
