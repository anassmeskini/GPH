#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <optional>
#include <string>

struct ArgInfo
{
   std::string probFile;
   int timelimit;
   int nthreads;
#ifndef NDEBUG
   int verbosity;
#endif
};

std::optional<ArgInfo>
parseArgs(int argc, char** argv);

#endif
