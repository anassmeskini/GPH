#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "Message.h"

#include <cstring>
#include <fstream>
#include <optional>
#include <string>
#include <variant>
#include <vector>

struct ParamValuePair
{
   std::string heuristic;
   std::string param;
   std::variant<std::string, int, double> value;
};

// format
// HeuristicName/Param = value

class Config
{
 public:
   Config(const std::string& config_file)
   {
      std::ifstream in;
      if (config_file.empty())
      {
         in.open("conf.txt");
         if (!in.is_open())
            return;

         Message::warn("using the configuration file conf.txt");
      }
      else
      {
         in.open(config_file);
         if (!in.is_open())
            throw std::runtime_error("unable to open config file");
      }

      bool format_error = false;
      char line[256];
      char* field1 = nullptr;
      char* field2 = nullptr;
      char* field3 = nullptr;
      while (in.getline(line, sizeof(line)))
      {
         if (line[0] == '\0' || line[0] == '#')
            continue;

         ParamValuePair pair;
         field1 = std::strtok(line, "/");
         if (!field1)
         {
            format_error = true;
            break;
         }

         field2 = std::strtok(nullptr, " =");
         if (!field2)
         {
            format_error = true;
            break;
         }

         field3 = std::strtok(nullptr, " =");
         if (!field3)
         {
            format_error = true;
            break;
         }

         pair.heuristic = std::string(field1);
         pair.param = std::string(field2);

         size_t value_size = 0;
         char* f3_copy = field3;
         while (*(f3_copy++) != '\0')
            ++value_size;

         size_t conv_size = 0;
         int int_param = std::stoi(field3, &conv_size);
         if (conv_size < value_size)
         {
            double double_param = std::stof(field3, &conv_size);
            if (conv_size < value_size)
               pair.value = std::string(field3);
            else
               pair.value = double_param;
         }
         else
            pair.value = int_param;

         config.push_back(pair);
      }

      if (format_error)
         throw std::runtime_error("format error in config file");
   }

   auto begin() const { return config.begin(); };
   auto end() const { return config.end(); };

 private:
   std::vector<ParamValuePair> config;
};

#endif
