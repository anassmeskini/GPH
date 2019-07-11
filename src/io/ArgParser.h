#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <optional>
#include <string>

struct ArgInfo
{
   std::string probFile;
   std::string outSolFile;
   std::string inSolFile;

   int timelimit;
   int nthreads;
};

std::optional<ArgInfo>
parseArgs(int argc, char** argv);

#endif
