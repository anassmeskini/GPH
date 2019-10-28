#include "core/Common.h"
#include "core/Heuristic.h"
#include "core/LPSolver.h"
#include "core/MIP.h"
#include "core/MySolver.h"
#include "core/Timer.h"

#include "io/ArgParser.h"
#include "io/MPSReader.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

#include "methods/BoundSolution.h"
#include "methods/CoefDiving.h"
#include "methods/DivingHeuristic.h"
#include "methods/FeasPump.h"
#include "methods/FracDiving.h"
#include "methods/IntShifting.h"
#include "methods/MinFracRounding.h"
#include "methods/MinLockRounding.h"
#include "methods/Octane.h"
#include "methods/RandRounding.h"
#include "methods/Shifting.h"
#include "methods/VecLengthDiving.h"

#include <cassert>
#include <exception>
#include <new>
#include <optional>
#include <filesystem>
#include <tbb/task_scheduler_init.h>

int
main(int argc, char** argv)
{
   // read arguments
   std::optional optionalArgs = parseArgs(argc, argv);

   if (!optionalArgs)
      return 1;

   const auto& args = optionalArgs.value();

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

   Search search{new BoundSolution,   new MinFracRounding,
                 new MinLockRounding, new Shifting,
                 new IntShifting,     new CoefDiving,
                 new FracDiving,      new VecLengthDiving,
                 new FeasPump,        new RandRounding};

   std::optional solution = search.run(mip, args.timelimit);

   // write the solution to disk
   if (solution)
   {
      std::vector solVec = solution.value();
      assert(solVec.size() == static_cast<size_t>(mip.getNCols()));

      std::string filename = std::filesystem::path(args.probFile).filename();
      SOLFormat::write(filename, solVec, mip.getVarNames());
   }

   return 0;
}
