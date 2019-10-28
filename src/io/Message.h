#ifndef MESSAGE_HPP
#define MESSAGE_HPP

#include "fmt/format.h"
#include <string_view>

#ifndef NDEBUG
#include <mutex>
#include <tbb/mutex.h>
#endif

struct Message
{
   template <typename... Args>
   static void print(std::string_view str, Args&&... args)
   {
      fmt::print(str, std::forward<Args>(args)...);
      fmt::print("\n");
   }

   template <typename... Args>
   static void print(std::FILE* file, std::string_view str, Args&&... args)
   {
      fmt::print(file, str, std::forward<Args>(args)...);
   }

   template <typename... Args>
   static void warn(std::string_view str, Args&&... args)
   {
      fmt::print("Warning: ");
      fmt::print(str, std::forward<Args>(args)...);
      fmt::print("\n");
   }

   template <typename... Args>
   static void error(std::string_view str, Args&&... args)
   {
      fmt::print(stderr, "Error: ");
      fmt::print(stderr, str, std::forward<Args>(args)...);
      fmt::print("\n");
   }

#ifndef NDEBUG
   inline static tbb::mutex mut;
   inline static enum Verbosity {
      QUIET,
      DEBUG,
      DEBUG_DETAILS
   } verbosity = QUIET;

   template <typename... Args>
   static void debug(std::string_view str, Args&&... args)
   {
      if (verbosity != QUIET)
      {
         std::lock_guard<tbb::mutex> lock(mut);
         fmt::print("[Debug]: ");
         fmt::print(str, std::forward<Args>(args)...);
         fmt::print("\n");
      }
   }

   template <typename... Args>
   static void debug_details(std::string_view str, Args&&... args)
   {
      std::lock_guard<tbb::mutex> lock(mut);
      if (verbosity == DEBUG_DETAILS)
      {
         fmt::print("[Debug]: ");
         fmt::print(str, std::forward<Args>(args)...);
         fmt::print("\n");
      }
   }
#else

   template <typename... Args>
   static void debug(std::string_view, Args&&...)
   {
   }

   template <typename... Args>
   static void debug_details(std::string_view, Args&&...)
   {
   }

#endif
};

#endif
