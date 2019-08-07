#include "ArgParser.h"
#include "clipp/clipp.h"

#include <iostream>

std::optional<ArgInfo> parseArgs(int argc, char **argv)
{

   using namespace clipp;

   ArgInfo arginfo;

   arginfo.timelimit = -1;
   arginfo.nthreads = -1;
   arginfo.probFile = "mip.mps";
   arginfo.verbosity = 2;

   auto cli =
       (value("input file", arginfo.probFile),
        option("-l", "--limit") &
            value("tlimit", arginfo.timelimit).doc("time limit"),
        option("-o", "--output") & value("outfile", arginfo.inSolFile)
                                       .doc("output file for the solution"),
        option("-v") & value("verbosity", arginfo.verbosity),
        option("-t", "--thread") & value("nthreads", arginfo.nthreads)
                                       .doc("number of threads to use"));

   if (!parse(argc, argv, cli))
   {
      std::cout << make_man_page(cli, argv[0]);
      return {};
   }

   return arginfo;
}
