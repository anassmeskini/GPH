#include <stdlib.h>
#include "AvaiLPSolver.h"
#include "CPXSolver.h"
#include "Common.h"
#include "LPSolver.h"
#include "MIP.h"
#include "mpsreader.h"

int main() {
   MIP<double> mip;

   try {
      mip = mpsreader::parse("mip.mps");
   } catch (const std::exception& ex) {
      std::cout << ex.what();
   }

   LPSolver<double>* solver;

   try {
      solver = new AvaiLPSolver(mip);
   } catch (...) {
      std::cout << "Solver raised an exception";
   }

   auto result = solver->solve();

   std::cout << "LP solver return status: " << to_str(result.status)
             << std::endl;

   if (result.status == LPResult::OPTIMAL) {
      std::cout << "obj: " << result.obj << std::endl;

      std::cout << "primal solution: ";
      for (auto val : result.primalSolution) std::cout << val << ", ";

      std::cout << "\ndual values: ";
      for (auto val : result.dualSolution) std::cout << val << ", ";

      bool feasible =
          checkFeasibility<double>(mip, result.primalSolution, 1e-9, 1e-6);

      std::cout << "\nfeasiblity check: " << feasible << std::endl;
   }

   return 0;
}
