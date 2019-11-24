#include "core/Common.h"
#include "core/Heuristic.h"
#include "core/LPSolver.h"
#include "core/MIP.h"
#include "core/MySolver.h"
#include "core/Timer.h"

#include "io/ArgParser.h"
#include "io/Config.h"
#include "io/MPSReader.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

#include "methods/BinaryLocalSearch.h"
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
#include <filesystem>
#include <new>
#include <optional>
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

   // read mip
   auto t0 = Timer::now();
   MIP mip = MPSReader::parse(args.probFile);
   auto t1 = Timer::now();

   Message::print("Reading the problem took: {:0.2f} sec.",
                  Timer::seconds(t1, t0));

   // takes ownership
   Search search(
       // feasibility heuristics
       {new BoundSolution, new IntShifting, new MinFracRounding,
        new MinLockRounding, new CoefDiving, new FracDiving,
        new RandRounding, new VecLengthDiving, new FeasPump, new Shifting},
       // improvement heuristics
       {new BinLocalSearch},
       // configuration
       Config(args.configFile));

   std::optional solution = search.run(mip, args.timelimit);

   if (solution)
   {
      // write the solution to disk
      std::vector solVec = solution.value();
      assert(solVec.size() == static_cast<size_t>(mip.getNCols()));

      std::string filename =
          std::filesystem::path(args.probFile).filename();
      SOLFormat::write(filename + ".sol", solVec, mip.getVarNames());
   }

   return 0;
}
