#ifndef _HEURISTIC_HPP
#define _HEURISTIC_HPP

#include "Common.h"
#include "LPSolver.h"
#include "MIP.h"
#include "MySolver.h"

#include <memory>
#include <optional>
#include <vector>

#include <tbb/mutex.h>
class SolutionPool
{
    struct Solution;

  public:
    SolutionPool() = default;

    void add(Solution &&sol)
    {
        std::unique_lock<tbb::mutex> lock(mut);
        pool.emplace_back(std::forward<Solution>(sol));
    }

    void add(std::vector<double> &&sol, double obj)
    {
        std::unique_lock<tbb::mutex> lock(mut);
        pool.emplace_back(std::forward<std::vector<double>>(sol), obj);
    }

  private:
    struct Solution
    {
        Solution(std::vector<double> &&sol, double obj)
            : solution(std::move(sol)), objective(obj)
        {
        }

        std::vector<double> solution;
        double objective;
    };

    tbb::mutex mut;
    std::vector<Solution> pool;

    friend class Search;
};

class HeuristicMethod
{
  public:
    virtual void
    search(const MIP &,                   // original problem
           const std::vector<double> &,   // lb at the node
           const std::vector<double> &,   // ub at the node
           const std::vector<Activity> &, // activities
           const LPResult &,              // LP solution at the current node
           const std::vector<double>
               &, // activities of the rows at the LP solution
           const std::vector<int> &, // integer variables with fractional values
           std::shared_ptr<const LPSolver>, // lp solver
           SolutionPool &) = 0;             // solution pool

    virtual ~HeuristicMethod() {}
};

class Search
{
  public:
    Search() = default;
    Search(std::initializer_list<HeuristicMethod *>);

    void run(const MIP &);

  private:
    std::vector<std::unique_ptr<HeuristicMethod>> heuristics;
    SolutionPool pool;
    std::vector<float> runtime;
};

#endif
