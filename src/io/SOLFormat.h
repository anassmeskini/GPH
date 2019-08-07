#ifndef SOLFORMAT_HPP
#define SOLFORMAT_HPP

#include "Message.h"
#include "ska/Hash.hpp"

#include <cassert>
#include <exception>
#include <fstream>
#include <string>
#include <vector>

struct SOLFormat
{
   static std::vector<double> read(const std::string& file,
                                   const std::vector<std::string>& colNames)
   {
      std::vector<double> solution(colNames.size(), 0.0);

      HashMap<std::string, int> nameToId;
      for (size_t i = 0; i < colNames.size(); ++i)
         nameToId.emplace(colNames[i], i);

      std::ifstream in(file);
      if (!in.is_open())
         throw std::runtime_error("unable to open file");

      char line[256];
      const char* delim = " \t";
      char* field1 = nullptr;
      char* field2 = nullptr;
      int linenb = 0;
      bool format_error = false;
      while (in.getline(line, sizeof(line)))
      {
         ++linenb;
         assert(std::strlen(line));

         if (line[0] == '#')
            continue;

         field1 = std::strtok(line, delim);
         if (!field1)
         {
            format_error = true;
            break;
         }

         field2 = std::strtok(line, delim);
         if (!field2)
         {
            format_error = true;
            break;
         }

         auto iter = nameToId.find(field1);
         if (iter != nameToId.end())
            solution[iter->second] = std::stod(field2);
         else
            Message::warn(
              "Skipping unknown variable <{}> in line {}", field1, linenb);
      }

      if (format_error)
         throw std::runtime_error("unable to parse file, error in line " +
                                  std::to_string(linenb));

      return solution;
   }

   static void write(const std::string& file,
                     const std::vector<double>& solution,
                     const std::vector<std::string>& colNames)
   {
      FILE* out = fopen(file.c_str(), "w");
      if (out == nullptr)
         throw std::runtime_error("unable to open output file");

      size_t nbzero = 0;

      for (size_t i = 0; i < solution.size(); ++i)
      {
         if (solution[i] == 0.0)
         {
            ++nbzero;
            continue;
         }

         Message::print(out, "{}   {:<15}\n", colNames[i], solution[i]);
      }

      fclose(out);

      if (nbzero == solution.size())
         Message::warn("all solution values are zero");

      Message::debug("solution written to {}", file);
   }
};

#endif
