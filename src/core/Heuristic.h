#ifndef _HEURISTIC_HPP
#define _HEURISTIC_HPP

#include "Common.h"
#include "LPFactory.h"
#include "MIP.h"
#include "MySolver.h"
#include "Problem.h"

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

class HeuristicMethod
{
   public:
   virtual std::optional<std::vector<double>> search(
     const MIP<double>&,         // original problem
     const std::vector<double>&, // lb at the node
     const std::vector<double>&, // ub at the node
     const LPResult&,            // LP solution at the current node
     const std::vector<double>&, // activities of the rows at the LP solution
     const std::vector<size_t>&, // integer variables with fractional values
     const LPFactory&) = 0;      // factory to copy the lp model

   virtual ~HeuristicMethod() {}
};

class Heuristics
{
   public:
   Heuristics() = default;
   Heuristics(std::initializer_list<HeuristicMethod*>);

   void run(MIP<double>&&);

   ~Heuristics();

   private:
   void cleanNodes() {}

   std::vector<std::unique_ptr<HeuristicMethod>> heuristics;
};

#endif
