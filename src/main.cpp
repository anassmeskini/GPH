#include <memory>
#include "core/AvaiLPSolver.h"
#include "core/Common.h"
#include "core/LPSolver.h"
#include "core/MIP.h"
#include "core/MPSReader.h"
#include "fmt/format.h"

int main(int argc, char** argv) {
   MIP<double> mip;
   std::string filename("mip.mps");

   if (argc == 2) filename = std::string(argv[1]);

   try {
      mip = mpsreader::parse(filename);
   } catch (const std::exception& ex) {
      std::cout << ex.what();
      return 1;
   }

   try {
      std::unique_ptr<LPSolver<double>> solver(new AvaiLPSolver(mip));
      LPResult result = solver->solve();

      fmt::print("LP solver return status:{}\n", to_str(result.status));

      if (result.status == LPResult::OPTIMAL) {
         fmt::print("obj: {}\n", result.obj);
      }

   } catch (...) {
      fmt::print("Solver raised an exception\n");
   }
   return 0;
}
