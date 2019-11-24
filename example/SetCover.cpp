
#include "GreedyHeuristic.h"

#include "io/ArgParser.h"
#include "io/MPSReader.h"
#include "io/SOLFormat.h"

#include "methods/BinaryLocalSearch.h"

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

   Search search({new GreedyHeuristic}, {new BinLocalSearch});
   search.run(mip, args.timelimit);
}
