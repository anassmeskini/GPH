#ifndef _HEURISTIC_HPP
#define _HEURISTIC_HPP

#include "Common.h"
#include "LPSolver.h"
#include "MIP.h"
#include "MySolver.h"

#include <memory>
#include <optional>
#include <vector>

struct Node
{
   Node() = default;
   Node(Node&) = default;
   Node(Node&&) = default;

   std::vector<double> lb;
   std::vector<double> ub;
   int depth;
};

struct Solution
{
   std::vector<double> solution;
   double objective;
};

class HeuristicMethod
{
   public:
   virtual std::optional<Solution> search(
     const MIP&,                   // original problem
     const std::vector<double>&,   // lb at the node
     const std::vector<double>&,   // ub at the node
     const std::vector<Activity>&, // activities
     const LPResult&,              // LP solution at the current node
     const std::vector<double>&,   // activities of the rows at the LP solution
     const std::vector<int>&,      // integer variables with fractional values
     std::shared_ptr<LPSolver>) = 0; // lp solver

   virtual ~HeuristicMethod() {}
};

class Search
{
   public:
   Search() = default;
   Search(std::initializer_list<HeuristicMethod*>);

   void run(const MIP&);

   ~Search();

   private:
   std::vector<std::unique_ptr<HeuristicMethod>> heuristics;
};

#endif
