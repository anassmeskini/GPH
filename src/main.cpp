#include "core/Common.h"
#include "core/Heuristic.h"
#include "core/LPSolver.h"
#include "core/MIP.h"
#include "core/MySolver.h"
#include "core/Timer.h"
#include "io/ArgParser.h"
#include "io/MPSReader.h"
#include "io/Message.h"
#include "methods/BoundSolution.h"
#include "methods/MinLockRounding.h"

#include <memory>

int
main(int argc, char** argv)
{
   std::optional<ArgInfo> optionalArgs = parseArgs(argc, argv);

   if (!optionalArgs)
      return 1;

   auto args = optionalArgs.value();

   MIP mip;

   try
   {
      auto t0 = Timer::now();
      mip = MPSReader::parse(args.probFile);
      auto t1 = Timer::now();

      Message::print("Reading the problem took: {:0.2f}s",
                     Timer::seconds(t1, t0));
   }
   catch (const std::exception& ex)
   {
      Message::print(ex.what());
      return 1;
   }

   printStats(mip.getStatistics());

   Search search({ new BoundSolution, new MinLockRounding });
   search.run(mip);

   return 0;
}
