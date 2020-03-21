#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <optional>
#include <string>

struct ArgInfo
{
   // input problem
   std::string probFile;
   // configuration file
   std::string configFile;
   // input solution
   std::string solutionFile;

   int timelimit;
   int nthreads;
   // write solution to a file
   bool writeSol;
#ifndef NDEBUG
   // output level
   int verbosity;
#endif
};

std::optional<ArgInfo>
parseArgs(int argc, char** argv);

#endif
