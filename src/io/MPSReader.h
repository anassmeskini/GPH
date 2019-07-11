#ifndef _MPS_READER_HPP
#define _MPS_READER_HPP

#include "core/MIP.h"
#include "dynamic_bitset/dynamic_bitset.hpp"
#include "ska/Hash.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <numeric>
#include <string>
#include <string_view>

class MPSWrapper
{
   public:
   MPSWrapper(std::istream& _is)
     : is(_is)
     , integer_section(false)
     , linenb(0)
   {
      buf[0] = '\0';
   }

   bool readLine() noexcept;

   const char* field1() const { return field_1; }
   const char* field2() const { return field_2; }
   const char* field3() const { return field_3; }
   const char* field4() const { return field_4; }
   const char* field5() const { return field_5; }
   const char* field6() const { return field_6; }

   bool isIntSection() const { return integer_section; }

   int getLineNb() const { return linenb; }

   private:
   std::istream& is;

   char buf[256];
   bool integer_section;

   char* field_1 = nullptr;
   char* field_2 = nullptr;
   char* field_3 = nullptr;
   char* field_4 = nullptr;
   char* field_5 = nullptr;
   char* field_6 = nullptr;

   int linenb = 0;

   static constexpr char blank = ' ';
   static constexpr char tab = '\t';
};

class MPSReader
{
   public:
   static MIP parse(const std::string&);

   private:
   enum Section : uint8_t
   {
      NAME,
      ROWS,
      COLUMNS,
      RHS,
      BOUNDS,
      RANGES,
      END,
      FORMAT_ERROR,
   };
   static Section error_section;

   static Section parseName(MPSWrapper&, std::string& name);

   static Section parseRows(MPSWrapper&, Rows& rows, std::string& objName);

   static Section parseColumns(MPSWrapper&,
                               const Rows& rows,
                               Cols& cols,
                               std::vector<double>& coefs,
                               std::vector<int>& idxT,
                               std::vector<int>& rstart,
                               std::vector<double>& obj,
                               const std::string& objName,
                               dynamic_bitset<>&,
                               std::vector<int>&,
                               std::vector<std::string>&);

   static Section parseRhs(MPSWrapper&,
                           const Rows& rows,
                           std::vector<double>& lhs,
                           std::vector<double>& rhs);

   static Section parseBounds(MPSWrapper&,
                              const Cols& cols,
                              std::vector<double>& lbs,
                              std::vector<double>& ubs,
                              dynamic_bitset<>& integer);

   static Section parseRanges(MPSWrapper&,
                              const Rows& rows,
                              std::vector<double>& lhs,
                              std::vector<double>& rbs);

   static SparseMatrix transpose(const SparseMatrix&, const std::vector<int>&);
};

#endif
