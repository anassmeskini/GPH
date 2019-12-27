#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <optional>
#include <string>

struct ArgInfo
{
   std::string probFile;
   std::string configFile;
   // input solution
   std::string solutionFile;

   int timelimit;
   int nthreads;
   bool writeSol;
#ifndef NDEBUG
   int verbosity;
#endif
};

std::optional<ArgInfo>
parseArgs(int argc, char** argv);

#endif
