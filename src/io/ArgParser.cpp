#include "ArgParser.h"
#include "Message.h"
#include "clipp/clipp.h"

#include <iostream>
#include <limits>

std::optional<ArgInfo>
parseArgs(int argc, char** argv)
{

   using namespace clipp;

   ArgInfo arginfo;

   arginfo.timelimit = std::numeric_limits<int>::max();
   arginfo.nthreads = -1;
   arginfo.probFile = "mip.mps";
   arginfo.writeSol = false;

#ifndef NDEBUG
   arginfo.verbosity = 2;
#endif

   auto cli =
       (value("input file", arginfo.probFile),
        option("-l", "--limit") & value("tlimit", arginfo.timelimit)
                                      .doc("time limit in seconds"),
#ifndef NDEBUG
        option("-v") & value("verbosity", arginfo.verbosity),
#endif
        option("-t", "--thread") & value("nthreads", arginfo.nthreads)
                                       .doc("number of threads to use"),
        option("-w").set(arginfo.writeSol).doc("write solution to disk"),
        option("-s", "--solution") &
            value("start_sol", arginfo.solutionFile)
                .doc("path to solution to improve"),
        option("-c", "--config") &
            value("config", arginfo.configFile).doc("configuration file"));

#ifndef NDEBUG
   // set verbosity level
   switch (arginfo.verbosity)
   {
   case 1:
      Message::verbosity = Message::QUIET;
      break;
   case 2:
      Message::verbosity = Message::DEBUG;
      break;
   case 3:
      Message::verbosity = Message::DEBUG_DETAILS;
      break;
   }
#endif

   if (!parse(argc, argv, cli))
   {
      std::cout << make_man_page(cli, argv[0]);
      return {};
   }

   return arginfo;
}
