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
#include "methods/CoefDiving.h"
#include "methods/DivingHeuristic.h"
#include "methods/FracDiving.h"
#include "methods/IntShifting.h"
#include "methods/MinLockRounding.h"
#include "methods/Shifting.h"
#include "methods/VecLengthDiving.h"

#include <cassert>
#include <exception>
#include <new>
#include <optional>
#include <tbb/task_scheduler_init.h>

int
main(int argc, char** argv)
{
   // read arguments
   std::optional<ArgInfo> optionalArgs = parseArgs(argc, argv);

   if (!optionalArgs)
      return 1;

   auto args = optionalArgs.value();

   // set verbosity level
   switch (args.verbosity)
   {
   case 1:
      Message::verbosity = Message::RELEASE;
      break;
   case 2:
      Message::verbosity = Message::DEBUG;
      break;
   case 3:
      Message::verbosity = Message::DEBUG_DETAILS;
      break;
   }

   // set number of threads
   assert(args.nthreads == -1 || args.nthreads >= 1);
   tbb::task_scheduler_init init(args.nthreads);

   MIP mip;

   // read mip
   try
   {
      auto t0 = Timer::now();
      mip = MPSReader::parse(args.probFile);
      auto t1 = Timer::now();

      Message::print("Reading the problem took: {:0.2f} sec.",
                     Timer::seconds(t1, t0));
   }
   catch (const std::exception& ex)
   {
      Message::print(ex.what());
      return 1;
   }

   printStats(mip.getStatistics());

   Search search{new BoundSolution,  new MinLockRounding, new Shifting,
                 new IntShifting,    new CoefDiving,      new FracDiving,
                 new VecLengthDiving};
   search.run(mip);

   return 0;
}
