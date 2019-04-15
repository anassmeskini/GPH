#include "ArgParser.h"
#include <iostream>

std::optional<ArgInfo>
parseArgs(int argc, char** argv)
{

   using namespace clipp;

   ArgInfo arginfo;

   arginfo.timelimit = -1;
   arginfo.nthreads = -1;
   arginfo.probFile = "mip.mps";

   auto cli =
     (option("-i", "--in") &
        value("infile", arginfo.probFile).doc("input problem file"),
      option("-l", "--limit") &
        value("tlimit", arginfo.timelimit).doc("time limit"),
      option("-o", "--output") &
        value("outfile", arginfo.inSolFile).doc("output file for the solution"),
      option("-t", "--thread") &
        value("nthreads", arginfo.nthreads).doc("number of threads to use"));

   if (!parse(argc, argv, cli))
   {
      std::cout << make_man_page(cli, argv[0]);
      return {};
   }

   return arginfo;
}
