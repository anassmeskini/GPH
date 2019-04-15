#ifndef _HEURISTIC_HPP
#define _HEURISTIC_HPP

#include "AvaiLPSolver.h"
#include "Common.h"
#include "MIP.h"
#include "Problem.h"

#include <memory>
#include <optional>
#include <vector>

class HeuristicMethod
{
   public:
   virtual std::optional<std::vector<double>> search(const ProblemView&) = 0;

   virtual ~HeuristicMethod() {}
};

class Heuristics
{
   public:
   Heuristics() = default;

   Heuristics(const std::initializer_list<HeuristicMethod*>&);

   std::optional<std::vector<double>> run(MIP<double>&&);

   private:
   std::vector<std::unique_ptr<HeuristicMethod>> heuristics;
};

#endif
