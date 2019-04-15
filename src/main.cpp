#include "core/AvaiLPSolver.h"
#include "core/Common.h"
#include "core/Heuristic.h"
#include "core/LPSolver.h"
#include "core/MIP.h"
#include "core/Problem.h"
#include "core/Timer.h"
#include "fmt/format.h"
#include "io/ArgParser.h"
#include "io/MPSReader.h"
#include "methods/TrivialSolutions.h"

#include <memory>

int
main(int argc, char** argv)
{
   auto optarginfo = parseArgs(argc, argv);

   if (!optarginfo)
      return 1;

   auto arginfo = optarginfo.value();

   MIP<double> mip;
   std::string filename(arginfo.probFile);

   try
   {
      auto t0 = Timer::now();
      mip = mpsreader::parse(filename);
      auto t1 = Timer::now();

      fmt::print("Reading the problem took: {}\n", Timer::seconds(t1, t0));
   }
   catch (const std::exception& ex)
   {
      std::cout << ex.what();
      return 1;
   }

   Heuristics heur({ new TrivialSolutions });

   heur.run(std::move(mip));

   return 0;
}
