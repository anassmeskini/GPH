#ifndef _HEURISTIC_HPP
#define _HEURISTIC_HPP

#include "Common.h"
#include "LPSolver.h"
#include "MIP.h"
#include "MySolver.h"
#include "Timer.h"

#include <memory>
#include <optional>
#include <vector>

#include <tbb/mutex.h>

class SolutionPool
{
 public:
   using value_type = std::pair<std::vector<double>, double>;

   void add(std::vector<double>&& sol, double obj)
   {
      solution_list.emplace_back(sol, obj);
   }

   size_t size() const { return solution_list.size(); }

   const value_type& operator[](size_t id) const
   {
      return solution_list[id];
   }

 private:
   std::vector<value_type> solution_list;
};

class HeuristicMethod
{
 public:
   HeuristicMethod(std::string_view heur_name)
       : name(heur_name), runtime(0.0)
   {
   }

   virtual ~HeuristicMethod() {}

   void execute(const MIP& mip,                   // original problem
                const std::vector<double>& lb,    // lb at the node
                const std::vector<double>& ub,    // ub at the node
                const std::vector<Activity>& act, // activities
                const LPResult& res, // LP solution at the current node
                const std::vector<double>& lpsol, // activities of the rows
                                                  // at the LP solution
                const std::vector<int>&
                    integer, // integer variables with fractional values
                std::shared_ptr<const LPSolver> lpsolver, // lp solver
                SolutionPool& pool)
   {
      auto t0 = Timer::now();
      search(mip, lb, ub, act, res, lpsol, integer, lpsolver, pool);
      auto t1 = Timer::now();

      runtime += Timer::seconds(t1, t0);
   }

   const std::string& getName() const { return name; }

   float getRunTime() const { return runtime; }

 protected:
   virtual void search(
       const MIP&,                   // original problem
       const std::vector<double>&,   // lb at the node
       const std::vector<double>&,   // ub at the node
       const std::vector<Activity>&, // activities
       const LPResult&,              // LP solution at the current node
       const std::vector<double>&,   // activities of the rows at the LP
                                     // solution
       const std::vector<int>&, // integer variables with fractional values
       std::shared_ptr<const LPSolver>, // lp solver
       SolutionPool&) = 0;              // solution pool

   std::string name;
   float runtime;
};

class Search
{
 public:
   Search() = default;
   Search(std::initializer_list<HeuristicMethod*>);

   void run(const MIP&);

 private:
   std::vector<std::unique_ptr<HeuristicMethod>> heuristics;
   std::vector<SolutionPool> heuristics_solutions;
};

#endif
