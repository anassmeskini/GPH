#include "MPSReader.h"
#include "Message.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Timer.h"
#include <cstring>
#include <exception>

#ifdef ZLIB_FOUND
#include <zlib.h>

static std::string
decompress_string(const std::string& str)
{
   z_stream zs;
   std::memset(&zs, 0, sizeof(zs));

   if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK)
      throw std::runtime_error("inflateInit failed while decompressing.");

   zs.next_in = (unsigned char*)str.data();
   zs.avail_in = str.size();

   int ret;
   char outbuffer[32768];
   std::string outstring;

   do
   {
      zs.next_out = reinterpret_cast<unsigned char*>(outbuffer);
      zs.avail_out = sizeof(outbuffer);

      ret = inflate(&zs, 0);

      if (outstring.size() < zs.total_out)
         outstring.append(outbuffer, zs.total_out - outstring.size());

   } while (ret == Z_OK);

   inflateEnd(&zs);

   if (ret != Z_STREAM_END)
      throw std::runtime_error("Exception during zlib decompression");

   return outstring;
}

#else

static std::string
decompress_string(const std::string&)
{
   throw std::runtime_error("The library was built without zlib");

   return {};
}

#endif

MPSWrapper::MPSWrapper(const std::string& filename)
    : ptr(0), integer_section(false), linenb(0)
{
   std::fstream in(filename);
   in.seekg(0, std::ios::end);
   std::streampos size = in.tellg();
   in.seekg(0, std::ios::beg);
   content.resize(size);
   in.read(&content[0], size);
   in.close();

   auto pos = filename.find_last_of('.');

   if (filename.substr(pos + 1) == "gz")
   {
      std::string temp = decompress_string(content);
      content = std::move(temp);
   }

   buf = nullptr;
   next = &content[0];

   for (size_t i = 0; i < content.size(); ++i)
   {
      if (content[i] == tab)
         content[i] = blank;
   }
}

bool
MPSWrapper::getLine() noexcept
{
   buf = next;

   if (ptr == content.size())
      return false;

   do
   {
      ++next;
      ++ptr;
   } while (ptr < content.size() && *next != '\n');

   if (ptr < content.size())
   {
      *next = '\0';
      ++next;
      ++ptr;
   }

   return true;
}

bool
MPSWrapper::readLine() noexcept
{
   bool marker;
   do
   {
      marker = false;

      field_1 = nullptr;
      field_2 = nullptr;
      field_3 = nullptr;
      field_4 = nullptr;
      field_5 = nullptr;

      do
      {
         ++linenb;
         if (!getLine())
            return false;

      } while (buf[0] == '*' || buf[0] == '\n');

      // new section
      if (buf[0] != blank)
      {
         field_1 = std::strtok(buf, " ");
         if (!std::strcmp(field_1, "NAME"))
            field_2 = std::strtok(nullptr, " ");
         return true;
      }

      // not a comment, empty line or a new section
      // if it's a marker change integer_section and
      // move to the next line
      // else get the tokens and return
      do
      {
         field_1 = std::strtok(buf + 1, " ");
         assert(field_1);

         if (!(field_2 = std::strtok(nullptr, " ")))
            break;

         if (!std::strcmp(field_2, "'MARKER'"))
            marker = true;

         if (marker)
         {
            if (!(field_3 = std::strtok(nullptr, " ")))
               return false;

            if (!std::strcmp(field_3, "'INTORG'"))
               integer_section = true;
            else if (!std::strcmp(field_3, "'INTEND'"))
               integer_section = false;
            else
               return false;

            break;
         }

         if (!(field_3 = std::strtok(nullptr, " ")))
            break;

         if (!(field_4 = std::strtok(nullptr, " ")))
            break;

         field_5 = std::strtok(nullptr, " ");

         break;
      } while (false);
   } while (marker);

   return true;
}

// read the mps, fill a transposed constraint matrix
// construct the sparse constraint matrix from the transposed
MIP
MPSReader::parse(const std::string& file)
{
   std::fstream in(file);
   if (!in.is_open())
      throw std::runtime_error("unable to open file");

   MPSWrapper mps(file);

   std::string name;

   Rows rows;
   Cols cols;

   // transposed sparse matrix
   std::vector<double> coefsT;
   std::vector<int> idxT;
   std::vector<int> rstatrtT;

   std::vector<double> rhs;
   std::vector<double> lhs;

   std::vector<double> lbs;
   std::vector<double> ubs;

   std::vector<std::string> varNames;
   std::vector<std::string> consNames;

   std::vector<double> objective;
   std::string objName;

   dynamic_bitset<> integer;

   std::vector<int> rowSize;

   Timer::time_point t0;
   Timer::time_point t1;

   Section nextsection = NAME;
   while (nextsection != FORMAT_ERROR && nextsection != END)
   {
      switch (nextsection)
      {
      case NAME:
         Message::debug("section name");
         nextsection = parseName(mps, name);
         break;

      case ROWS:
         t0 = Timer::now();
         nextsection = parseRows(mps, rows, objName);
         t1 = Timer::now();
         Message::debug("Section ROWS parsed in {:0.2f}s",
                        Timer::seconds(t1, t0));
         break;

      case COLUMNS:
         t0 = Timer::now();
         nextsection =
             parseColumns(mps, rows, cols, coefsT, idxT, rstatrtT,
                          objective, objName, integer, rowSize, varNames);
         t1 = Timer::now();
         Message::debug("Section COLUMNS parsed in {:0.2f}s",
                        Timer::seconds(t1, t0));
         break;

      case RHS:
         t0 = Timer::now();
         nextsection = parseRhs(mps, rows, lhs, rhs);
         t1 = Timer::now();
         Message::debug("Section RHS parsed in {:0.2f}s",
                        Timer::seconds(t1, t0));
         break;

      case BOUNDS:
         t0 = Timer::now();
         nextsection = parseBounds(mps, cols, lbs, ubs, integer);
         t1 = Timer::now();
         Message::debug("Section BOUNDS parsed in {:0.2f}s",
                        Timer::seconds(t1, t0));
         break;

      case RANGES:
         t0 = Timer::now();
         nextsection = parseRanges(mps, rows, lhs, rhs);
         t1 = Timer::now();
         Message::debug("Section BOUNDS parsed in {:0.2f}s",
                        Timer::seconds(t1, t0));
         break;

      case FORMAT_ERROR:
      case END:
         assert(0);
      }
   }

   if (nextsection == FORMAT_ERROR)
      throw std::runtime_error("unable to parse MPS file (error in line " +
                               std::to_string(mps.getLineNb()) + ")");

   return MIP(rows, cols, std::move(coefsT), std::move(idxT),
              std::move(rstatrtT), std::move(rhs), std::move(lhs),
              std::move(lbs), std::move(ubs), std::move(objective),
              std::move(integer), rowSize, std::move(varNames));
}

MPSReader::Section
MPSReader::parseName(MPSWrapper& mps, std::string& name)
{
   if (!mps.readLine() || std::strcmp(mps.field1(), "NAME"))
      return FORMAT_ERROR;

   name = std::string(mps.field2());

   if (!mps.readLine() || std::strcmp(mps.field1(), "ROWS"))
      return FORMAT_ERROR;

   return ROWS;
}

MPSReader::Section
MPSReader::parseRows(MPSWrapper& mps, Rows& rows, std::string& objname)
{
   int rowcounter = 0;
   ConsType type;

   while (true)
   {
      if (!mps.readLine())
         return FORMAT_ERROR;

      if (!mps.field2())
         break;

      if (std::strlen(mps.field1()) > 1)
         return FORMAT_ERROR;

      switch (mps.field1()[0])
      {
      case 'N':
         type = OBJECTIVE;
         objname = std::string(mps.field2());
         break;
      case 'L':
         type = LESS;
         break;
      case 'G':
         type = GREATER;
         break;
      case 'E':
         type = EQUAL;
         break;
      default:
         return FORMAT_ERROR;
      }

      if (type != OBJECTIVE)
      {
         auto pair =
             rows.emplace(mps.field2(), std::make_pair(type, rowcounter));
         ++rowcounter;

         // duplicate rows
         if (!pair.second)
            return FORMAT_ERROR;
      }
   }

   if (std::strcmp(mps.field1(), "COLUMNS"))
      return FORMAT_ERROR;
   return COLUMNS;
}

MPSReader::Section
MPSReader::parseColumns(MPSWrapper& mps, const Rows& rows, Cols& cols,
                        std::vector<double>& coefs, std::vector<int>& idxT,
                        std::vector<int>& rstart,
                        std::vector<double>& objective,
                        const std::string& objname,
                        dynamic_bitset<>& integer,
                        std::vector<int>& rowSize,
                        std::vector<std::string>& varNames)
{
   int colid = -1;
   std::string row_name;
   objective.clear();

   rowSize = std::vector<int>(rows.size(), 0);

   while (true)
   {
      if (!mps.readLine())
         return FORMAT_ERROR;

      if (!mps.field2())
         break;

      if (!mps.field3())
         return FORMAT_ERROR;

      if (varNames.empty() ||
          std::strcmp(mps.field1(), varNames.back().c_str()))
      {
         ++colid;
         std::string current_col(mps.field1());

         auto insertion = cols.emplace(current_col, colid);
         if (!insertion.second)
            return FORMAT_ERROR;

         // the last row of the transposed matrix ends here
         rstart.push_back(coefs.size());
         integer.push_back(mps.isIntSection());

         varNames.push_back(std::move(current_col));

         assert(varNames.size() == static_cast<size_t>(colid + 1));
      }
      assert(coefs.size() == idxT.size());

      double coef;
      try
      {
         coef = std::stod(mps.field3());
      }
      catch (const std::exception& ex)
      {
         return FORMAT_ERROR;
      }

      if (!std::strcmp(objname.c_str(), mps.field2()))
      {
         int beforesize = objective.size();
         objective.resize(colid + 1);

         assert(colid >= beforesize);
         assert(objective.size() == static_cast<size_t>(colid) + 1);

         std::memset(objective.data() + beforesize, 0,
                     sizeof(double) * (colid - beforesize));

         objective[colid] = coef;
      }
      else
      {
         auto iter = rows.find(mps.field2());
         // row not declared in the ROWS section
         if (iter == rows.end())
         {
            return FORMAT_ERROR;
         }

         auto rowid = iter->second.second;

         if (coef != 0.0)
         {
            coefs.push_back(coef);
            idxT.push_back(rowid);
            ++rowSize[rowid];
         }

         assert(coefs.size() == idxT.size());
      }

      if (!mps.field4())
         continue;
      if (!mps.field5())
         return FORMAT_ERROR;

      try
      {
         coef = std::stod(mps.field5());
      }
      catch (const std::exception& ex)
      {
         return FORMAT_ERROR;
      }

      if (!std::strcmp(objname.c_str(), mps.field4()))
      {
         int beforesize = objective.size();
         objective.resize(colid + 1);

         assert(colid >= beforesize);
         assert(objective.size() == static_cast<size_t>(colid) + 1);

         std::memset(objective.data() + beforesize, 0,
                     sizeof(double) * (colid - beforesize));

         objective[colid] = coef;
      }
      else
      {
         auto iter = rows.find(mps.field4());
         // row not declared in the ROWS section
         if (iter == rows.end())
            return FORMAT_ERROR;

         auto rowid = iter->second.second;
         if (coef != 0.0)
         {
            coefs.push_back(coef);
            idxT.push_back(rowid);
            ++rowSize[rowid];
         }

         assert(coefs.size() == idxT.size());
      }
   }

   // end the last row
   rstart.push_back(coefs.size());

   int ncols = static_cast<int>(colid + 1);
   int beforesize = objective.size();
   objective.resize(ncols);

   assert(ncols >= beforesize);
   std::memset(objective.data() + beforesize, 0,
               sizeof(double) * (ncols - beforesize));

   if (!std::strcmp(mps.field1(), "RHS"))
      return RHS;

   return FORMAT_ERROR;
}

MPSReader::Section
MPSReader::parseRhs(MPSWrapper& mps, const Rows& rows,
                    std::vector<double>& lhs, std::vector<double>& rhs)
{
   const double inf = std::numeric_limits<double>::infinity();

   int nrows = rows.size();
   lhs = std::vector<double>(nrows);
   rhs = std::vector<double>(nrows);

   // set default bounds
   for (auto row : rows)
   {
      int id = row.second.second;
      ConsType type = row.second.first;

      assert(id < static_cast<int>(lhs.size()));
      switch (type)
      {
      case LESS:
         lhs[id] = -inf;
         rhs[id] = 0.0;
         break;
      case GREATER:
         lhs[id] = 0.0;
         rhs[id] = inf;
         break;
      case EQUAL:
         lhs[id] = 0.0;
         rhs[id] = 0.0;
         break;
      case OBJECTIVE:
         break;
      }
   }

   do
   {
      if (!mps.readLine())
         return FORMAT_ERROR;

      if (!mps.field2())
         break;
      if (!mps.field3())
         return FORMAT_ERROR;

      auto iter = rows.find(mps.field2());
      if (iter == rows.end())
         return FORMAT_ERROR;

      int rowid = iter->second.second;

      double side;
      try
      {
         side = std::stod(mps.field3());
      }
      catch (const std::exception& ex)
      {
         return FORMAT_ERROR;
      }

      switch (iter->second.first)
      {
      case LESS:
         lhs[rowid] = -inf;
         rhs[rowid] = side;
         break;
      case GREATER:
         lhs[rowid] = side;
         rhs[rowid] = inf;
         break;
      case EQUAL:
         lhs[rowid] = side;
         rhs[rowid] = side;
         break;
      case OBJECTIVE:
         assert(0);
         break;
      }

      if (!mps.field4())
         continue;
      if (!mps.field5())
         return FORMAT_ERROR;

      iter = rows.find(mps.field4());
      if (iter == rows.end())
         return FORMAT_ERROR;

      rowid = iter->second.second;

      try
      {
         side = std::stod(mps.field5());
      }
      catch (const std::exception& ex)
      {
         return FORMAT_ERROR;
      }

      switch (iter->second.first)
      {
      case LESS:
         lhs[rowid] = -inf;
         rhs[rowid] = side;
         break;
      case GREATER:
         lhs[rowid] = side;
         rhs[rowid] = inf;
         break;
      case EQUAL:
         lhs[rowid] = side;
         rhs[rowid] = side;
         break;
      case OBJECTIVE:
         assert(0);
         break;
      }

   } while (true);

   if (!std::strcmp(mps.field1(), "BOUNDS"))
      return BOUNDS;
   if (!std::strcmp(mps.field1(), "RANGES"))
      return RANGES;

   return FORMAT_ERROR;
}

MPSReader::Section
MPSReader::parseBounds(MPSWrapper& mps, const Cols& cols,
                       std::vector<double>& lbs, std::vector<double>& ubs,
                       dynamic_bitset<>& integer)
{
   lbs = std::vector<double>(cols.size(), 0);
   ubs = std::vector<double>(cols.size(),
                             std::numeric_limits<double>::infinity());

   dynamic_bitset<> lb_changed(cols.size(), false);

   // TODO
   constexpr double inf = std::numeric_limits<double>::infinity();

   do
   {
      if (!mps.readLine())
         return FORMAT_ERROR;
      if (!mps.field2())
         break;
      if (!mps.field3())
         return FORMAT_ERROR;

      auto iter = cols.find(mps.field3());
      if (iter == cols.end())
         return FORMAT_ERROR;
      int colid = iter->second;
      if (mps.field4())
      {
         // TODO catch exceptions
         double bound;
         try
         {
            bound = std::stod(mps.field4());
         }
         catch (const std::exception& ex)
         {
            return FORMAT_ERROR;
         }

         if (!std::strcmp(mps.field1(), "UP"))
         {
            ubs[colid] = bound;
            if (bound < 0.0 && !lb_changed[colid])
               lbs[colid] = -inf;
         }
         else if (!std::strcmp(mps.field1(), "LO"))
         {
            lbs[colid] = bound;
            lb_changed[colid] = true;
         }
         else if (!std::strcmp(mps.field1(), "FX"))
         {
            lbs[colid] = bound;
            ubs[colid] = bound;
         }
         else if (!std::strcmp(mps.field1(), "MI"))
            lbs[colid] = -inf;
         else if (!std::strcmp(mps.field1(), "PL"))
            ubs[colid] = inf;
         else
            return FORMAT_ERROR;
      }
      else
      {
         if (!std::strcmp(mps.field1(), "FX"))
         {
            lbs[colid] = -inf;
            ubs[colid] = inf;
         }
         else if (!std::strcmp(mps.field1(), "BV"))
         {
            lbs[colid] = 0.0;
            ubs[colid] = 1.0;
            integer[colid] = true;
         }
      }

   } while (true);

   if (std::strcmp(mps.field1(), "ENDATA"))
      return FORMAT_ERROR;

   return END;
}

MPSReader::Section
MPSReader::parseRanges(MPSWrapper& mps, const Rows& rows,
                       std::vector<double>& lhs, std::vector<double>& rhs)
{
   std::string rangeVectorName;

   do
   {
      if (!mps.readLine())
         return FORMAT_ERROR;
      if (!mps.field2())
         break;
      if (!mps.field3())
         return FORMAT_ERROR;
      if (rangeVectorName.empty())
         rangeVectorName = std::string(mps.field1());
      else if (std::strcmp(rangeVectorName.c_str(), mps.field1()))
         // TODO warning
         continue;

      auto iter = rows.find(mps.field2());

      if (iter == rows.end())
         return FORMAT_ERROR;

      ConsType type = iter->second.first;
      int rowid = iter->second.second;
      double range;
      try
      {
         range = std::stod(mps.field3());
      }
      catch (const std::exception& ex)
      {
         return FORMAT_ERROR;
      }

      switch (type)
      {
      case GREATER:
         rhs[rowid] = lhs[rowid] + range;
         break;
      case LESS:
         lhs[rowid] = rhs[rowid] - range;
         break;
      case EQUAL:
         if (range > 0.0)
            rhs[rowid] = lhs[rowid] + range;
         else
            lhs[rowid] = rhs[rowid] + range;
         break;
      default:
         assert(0);
         break;
      }

   } while (true);

   if (!std::strcmp(mps.field1(), "BOUNDS"))
      return BOUNDS;

   return FORMAT_ERROR;
}
