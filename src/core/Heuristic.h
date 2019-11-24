#ifndef _HEURISTIC_HPP
#define _HEURISTIC_HPP

#include "Common.h"
#include "LPSolver.h"
#include "MIP.h"
#include "MySolver.h"
#include "Timer.h"

#include "io/Config.h"
#include "io/Message.h"

#include <memory>
#include <optional>
#include <variant>
#include <vector>

#include <optional>
#include <tbb/mutex.h>

class Heuristic
{
 public:
   explicit Heuristic(std::string_view heur_name) : name(heur_name) {}

   const std::string& getName() const { return name; }

   virtual void setParam(const std::string& param,
                         const std::variant<std::string, int, double>&)
   {
      Message::warn("Parameter {} has been ignored", param);
   }

 private:
   std::string name;
};

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

class FeasibilityHeuristic : public Heuristic
{
 public:
   explicit FeasibilityHeuristic(std::string_view heur_name)
       : Heuristic(heur_name), runtime(0.0)
   {
   }

   virtual ~FeasibilityHeuristic() {}

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
                TimeLimit limit,                          // time limit
                SolutionPool& pool)
   {
      auto t0 = Timer::now();
      search(mip, lb, ub, act, res, lpsol, integer, lpsolver, limit, pool);
      auto t1 = Timer::now();

      runtime += Timer::seconds(t1, t0);
   }

   float getRunTime() const { return runtime; }

 private:
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
       TimeLimit limit,                 // time limit
       SolutionPool&) = 0;              // solution pool

   float runtime;
};

class ImprovementHeuristic : public Heuristic
{
 public:
   explicit ImprovementHeuristic(std::string_view heur_name)
       : Heuristic(heur_name), runtime(0.0)
   {
   }

   virtual ~ImprovementHeuristic() {}

   void execute(
       const MIP& mip,                   // original problem
       const std::vector<double>& lb,    // lb at the node
       const std::vector<double>& ub,    // ub at the node
       const std::vector<Activity>& act, // activities
       const LPResult& res,              // LP solution at the current node
       const std::vector<double>& lpsol, // activities of the rows
                                         // at the LP solution
       const std::vector<int>&
           integer, // integer variables with fractional values
       const std::vector<double>& best_int_sol, // best integer solution
       double best_int_cost, // cost of the best solution
       std::shared_ptr<const LPSolver> lpsolver, // lp solver
       TimeLimit limit,                          // time limit
       SolutionPool& pool)
   {
      auto t0 = Timer::now();
      improve(mip, lb, ub, act, res, lpsol, integer, best_int_sol,
              best_int_cost, lpsolver, limit, pool);
      auto t1 = Timer::now();

      runtime += Timer::seconds(t1, t0);
   }

   float getRunTime() const { return runtime; }

 private:
   virtual void improve(
       const MIP&,                   // original problem
       const std::vector<double>&,   // lb at the node
       const std::vector<double>&,   // ub at the node
       const std::vector<Activity>&, // activities
       const LPResult&,              // LP solution at the current node
       const std::vector<double>&,   // activities of the rows at the LP
                                     // solution
       const std::vector<int>&, // integer variables with fractional values
       const std::vector<double>& best_int_sol, // best integer solution
       double best_int_cost,            // cost of the best solution
       std::shared_ptr<const LPSolver>, // lp solver
       TimeLimit limit,                 // time limit
       SolutionPool&) = 0;              // solution pool

   float runtime;
};

class Search
{
 public:
   Search(std::initializer_list<FeasibilityHeuristic*>,
          std::initializer_list<ImprovementHeuristic*>, const Config&);

   std::optional<std::vector<double>> run(const MIP&, int);

 private:
   std::tuple<int, int, double, int> getFeasSolSummary() const;

   std::tuple<int, int, double, int> getImprSolSummary() const;

   bool checkSolFeas(const MIP&) const;

 private:
   std::vector<std::unique_ptr<FeasibilityHeuristic>> feas_heuristics;
   std::vector<std::unique_ptr<ImprovementHeuristic>> impr_heuristics;
   std::vector<SolutionPool> feas_solutions_pools;
   std::vector<SolutionPool> impr_solutions_pools;
};

#endif
