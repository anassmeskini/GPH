#ifndef MESSAGE_HPP
#define MESSAGE_HPP

#include "fmt/format.h"
#include <string_view>

struct Message
{
   template<typename... Args>
   static void print(std::string_view str, Args... args)
   {
      fmt::print(str, std::forward<Args>(args)...);
      fmt::print("\n");
   }

   template<typename... Args>
   static void print(FILE* file, std::string_view str, Args... args)
   {
      fmt::print(file, str, std::forward<Args>(args)...);
      fmt::print("\n");
   }

   template<typename... Args>
   static void warn(std::string_view str, Args... args)
   {
      fmt::print("Warning: ");
      fmt::print(str, std::forward<Args>(args)...);
      fmt::print("\n");
   }

   template<typename... Args>
   static void error(std::string_view str, Args... args)
   {
      fmt::print(stderr, "Error: ");
      fmt::print(stderr, str, std::forward<Args>(args)...);
      fmt::print("\n");
   }

#ifndef NDEBUG
   template<typename... Args>
   static void debug(std::string_view str, Args... args)
   {
      fmt::print("[Debug]: ");
      fmt::print(str, std::forward<Args>(args)...);
      fmt::print("\n");
   }
#else
   // to prevent compiler warnings in release mode
   template<typename... Args>
   static void debug(std::string_view, Args...)
   {
   }
#endif
};

#endif
