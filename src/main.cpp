#include <memory>
#include "core/AvaiLPSolver.h"
#include "core/Common.h"
#include "core/LPSolver.h"
#include "core/MIP.h"
#include "core/mpsreader.h"

int main(int argc, char** argv) {
   MIP<double> mip;
   std::string filename("mip.mps");

   if (argc == 2) filename = std::string(argv[1]);

   try {
      mip = mpsreader::parse("mip.mps");
   } catch (const std::exception& ex) {
      std::cout << ex.what();
      return 1;
   }

   try {
      std::unique_ptr<LPSolver<double>> solver(new AvaiLPSolver(mip));
      LPResult result = solver->solve();

      std::cout << "LP solver return status: " << to_str(result.status)
                << std::endl;

      if (result.status == LPResult::OPTIMAL) {
         std::cout << "obj: " << result.obj << std::endl;
      }

   } catch (...) {
      std::cout << "Solver raised an exception";
   }
   return 0;
}
