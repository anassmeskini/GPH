#include "MPSReader.h"
#include "Timer.h"
#include "fmt/format.h"
#include <cstring>

std::string
string(const std::string& line, std::pair<int, int> word)
{
   return std::string(line.cbegin() + word.first, line.cbegin() + word.second);
}

// TODO improve decompression time

static Split
split(const std::string& line)
{
   constexpr char space = ' ';
   constexpr char tab = '\t';

   Split split;

   for (size_t id = 0; id < line.size(); ++id) {
      if (line[id] == space || line[id] == tab)
         continue;

      split.words[split.size_].first = id;

      while (true) {
         if (split.size_ >= 5)
            throw std::runtime_error("too many words");

         ++id;
         assert(id <= line.size());

         if (id >= line.size() || line[id] == space || line[id] == tab) {
            split.words[split.size_].second = id;
            ++split.size_;
            break;
         }
      }
   }

   return split;
}

bool
strcompare(const char* line,
           std::pair<size_t, size_t> word,
           const char* str,
           size_t size)
{
   if (word.second != word.first + size)
      return false;

   return !std::memcmp(line + word.first, str, size * sizeof(char));
}

mpsreader::Section
mpsreader::parseName(boost::iostreams::filtering_istream& file,
                     std::string& name)
{
   error_section = NAME;

   std::string line;
   while (std::getline(file, line) &&
          (line.empty() || line[line.find_first_not_of(" ")] == '*')) {
   }

   auto tokens = split(line);

   if (tokens.size() != 2 ||
       !strcompare(line.c_str(), tokens.words[0], "NAME", 4))
      return FAIL;

   name = string(line, tokens.words[1]);

   while (std::getline(file, line) &&
          (line.empty() || line[line.find_first_not_of(" ")] == '*')) {
   }

   tokens = split(line);
   if (tokens.size() > 1 ||
       !strcompare(line.c_str(), tokens.words[0], "ROWS", 4))
      return FAIL;

   error_section = NONE;
   return ROWS;
}

mpsreader::Section
mpsreader::parseRows(boost::iostreams::filtering_istream& file,
                     Rows& rows,
                     std::string& objName)
{
   error_section = ROWS;
   std::string line;
   Split tokens;

   objName.clear();
   size_t rowcounter = 0;
   while (std::getline(file, line)) {
      if (line.empty() || line[line.find_first_not_of(" ")] == '*')
         continue;

      tokens = split(line);
      if (tokens.size() == 1)
         break;
      if (tokens.size() != 2 || tokens.size(0) != 1)
         return FAIL;

      ConsType type;
      char typechar = line[tokens.words[0].first];
      switch (typechar) {
         case 'N':
            type = OBJECTIVE;
            objName = string(line, tokens.words[1]);
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
            return FAIL;
      }

      if (type != OBJECTIVE) {
         auto pair = rows.emplace(string(line, tokens.words[1]),
                                  std::make_pair(type, rowcounter));
         ++rowcounter;

         // duplicate rows
         if (!pair.second)
            return FAIL;
      }
   }

   assert(tokens.size() == 1);
   if (!strcompare(line.c_str(), tokens.words[0], "COLUMNS", 7))
      return FAIL;

   error_section = NONE;
   return COLUMNS;
}

mpsreader::Section
mpsreader::parseColumns(boost::iostreams::filtering_istream& file,
                        const Rows& rows,
                        Cols& cols,
                        std::vector<double>& coefs,
                        std::vector<size_t>& idxT,
                        std::vector<size_t>& rstart,
                        std::vector<double>& objective,
                        const std::string& objName,
                        bitset& integer,
                        std::vector<size_t>& rowSize,
                        std::vector<std::string>& varNames)
{
   error_section = COLUMNS;
   std::string line;
   Split tokens;

   int colId = -1;
   std::string prevCol("");

   bool integerSection = false;

   rowSize = std::vector<size_t>(rows.size(), 0);

   while (std::getline(file, line)) {
      assert(coefs.size() == idxT.size());

      if (line.empty() || line[line.find_first_not_of(" ")] == '*')
         continue;

      tokens = split(line);
      if (tokens.size() == 1)
         break;

      // check if it's the start or end of an integer section
      if (tokens.size() == 3) {
         if (strcompare(line.c_str(), tokens.words[1], "'MARKER'", 8)) {
            if (strcompare(line.c_str(), tokens.words[2], "'INTORG'", 8)) {
               if (integerSection)
                  return FAIL;
               integerSection = true;
               continue;
            } else if (strcompare(
                         line.c_str(), tokens.words[2], "'INTEND'", 8)) {
               if (!integerSection)
                  return FAIL;
               integerSection = false;
               continue;
            } else
               return FAIL;
         }
      } else if (!(tokens.size() % 2))
         return FAIL;

      auto curCol = string(line, tokens.words[0]);

      // if it's a new column
      if (curCol != prevCol) {
         ++colId;

         auto insertion = cols.insert(std::make_pair(curCol, colId));
         if (!insertion.second)
            return FAIL;

         // the last row of the transposed matrix ends here
         rstart.push_back(coefs.size());
         integer.push_back(integerSection);

         varNames.push_back(curCol);
         assert(varNames.size() == colId + 1);

         prevCol = std::move(curCol);
      }

      for (size_t i = 1; i < tokens.size(); i += 2) {
         auto rowname = string(line, tokens.words[i]);
         double coef = std::stod(string(line, tokens.words[i + 1]));

         // if the entry is for the objective
         if (rowname == objName) {
            if (objective.size() < static_cast<size_t>(colId)) {
               int objsize = objective.size();
               objective.resize(colId + 1);
               std::memset(objective.data() + objsize,
                           0.0,
                           (colId - objsize) * sizeof(double));
            }

            objective.resize(colId + 1);
            objective[colId] = coef;
            assert(objective.size() == static_cast<size_t>(colId) + 1);
         } else {
            auto iter = rows.find(rowname);
            // row not declared in the ROWS section
            if (iter == rows.end())
               return FAIL;

            auto rowId = iter->second.second;
            coefs.push_back(coef);
            idxT.push_back(rowId);
            ++rowSize[rowId];
         }
      }
   }

   while (objective.size() < static_cast<size_t>(colId) + 1)
      objective.push_back(0.0);

   // end the last row
   rstart.push_back(coefs.size());

   assert(tokens.size() == 1);
   if (!strcompare(line.c_str(), tokens.words[0], "RHS", 3))
      return FAIL;

   assert(objective.size() == cols.size());

   error_section = NONE;
   return RHS;
}

mpsreader::Section
mpsreader::parseRhs(boost::iostreams::filtering_istream& file,
                    const Rows& rows,
                    std::vector<double>& lhs,
                    std::vector<double>& rhs)
{
   error_section = RHS;
   std::string line;
   Split tokens;

   std::string prevCol("");
   std::set<std::string> colset;

   const double inf = std::numeric_limits<double>::infinity();

   assert(rows.size() > 1);
   size_t nrows = rows.size();
   lhs = std::vector<double>(nrows);
   rhs = std::vector<double>(nrows);

   for (auto row : rows) {
      size_t id = row.second.second;
      ConsType type = row.second.first;
      switch (type) {
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

   // set default bounds
   for (auto row : rows) {
      auto constype = row.second.first;
      auto id = row.second.second;

      switch (constype) {
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

   while (std::getline(file, line)) {
      if (line.empty() || line[line.find_first_not_of(" ")] == '*')
         continue;

      tokens = split(line);
      if (tokens.size() == 1)
         break;
      if (!(tokens.size() % 2))
         return FAIL;

      for (size_t i = 1; i < tokens.size(); i += 2) {
         auto rowname = string(line, tokens.words[i]);
         double side = std::stof(string(line, tokens.words[i + 1]));

         auto iter = rows.find(rowname);
         if (iter == rows.end())
            return FAIL;

         size_t rowid = iter->second.second;
         switch (iter->second.first) {
            case LESS:
               lhs[rowid] = -std::numeric_limits<double>::infinity();
               rhs[rowid] = side;
               break;
            case GREATER:
               lhs[rowid] = side;
               rhs[rowid] = std::numeric_limits<double>::infinity();
               break;
            case EQUAL:
               lhs[rowid] = side;
               rhs[rowid] = side;
               break;
            case OBJECTIVE:
               assert(0);
               break;
         }
      }
   }

   assert(tokens.size() == 1);
   if (!strcompare(line.c_str(), tokens.words[0], "BOUNDS", 6))
      return FAIL;

   error_section = NONE;
   return BOUNDS;
}

mpsreader::Section
mpsreader::parseBounds(boost::iostreams::filtering_istream& file,
                       const Cols& cols,
                       std::vector<double>& lbs,
                       std::vector<double>& ubs,
                       bitset& integer)
{
   error_section = BOUNDS;
   std::string line;
   Split tokens;

   lbs = std::vector<double>(cols.size(), 0);
   ubs =
     std::vector<double>(cols.size(), std::numeric_limits<double>::infinity());

   bitset lb_changed(cols.size(), false);

   // TODO
   constexpr double inf = std::numeric_limits<double>::infinity();

   while (std::getline(file, line)) {
      if (line.empty() || line[line.find_first_not_of(" ")] == '*')
         continue;

      tokens = split(line);
      if (tokens.size() == 1)
         break;
      if (tokens.size() < 3)
         return FAIL;

      auto colname = string(line, tokens.words[2]);
      auto iter = cols.find(colname);
      if (iter == cols.end())
         return FAIL;
      size_t colid = iter->second;

      if (tokens.size() == 4) {
         double bound = std::stof(string(line, tokens.words[3]));

         if (strcompare(line.c_str(), tokens.words[0], "UP", 2)) {
            ubs[colid] = bound;
            if (bound < 0.0 && !lb_changed[colid])
               lbs[colid] = -inf;
         } else if (strcompare(line.c_str(), tokens.words[0], "LO", 2)) {
            lbs[colid] = bound;
            lb_changed[colid] = true;
         } else if (strcompare(line.c_str(), tokens.words[0], "FX", 2)) {
            lbs[colid] = bound;
            ubs[colid] = bound;
         } else if (strcompare(line.c_str(), tokens.words[0], "MI", 2)) {
            lbs[colid] = -inf;
         } else if (strcompare(line.c_str(), tokens.words[0], "PI", 2)) {
            ubs[colid] = inf;
         } else
            return FAIL;
      } else if (tokens.size() == 3) {
         if (strcompare(line.c_str(), tokens.words[0], "FR", 2)) {
            lbs[colid] = -inf;
            ubs[colid] = inf;
         } else if (strcompare(line.c_str(), tokens.words[0], "BV", 2)) {
            integer[colid] = true;
            lbs[colid] = 0.0;
            ubs[colid] = 1.0;
         } else
            return FAIL;
      } else
         return FAIL;
   }

   assert(tokens.size() == 1);
   if (!strcompare(line.c_str(), tokens.words[0], "ENDATA", 6))
      return FAIL;

   error_section = NONE;
   return END;
}

MIP<double>
mpsreader::makeMip(const Rows& rows,
                   const Cols& cols,
                   std::vector<double>&& coefsT,
                   std::vector<size_t>&& idxT,
                   std::vector<size_t>&& rstartT,
                   std::vector<double>&& rhs,
                   std::vector<double>&& lhs,
                   std::vector<double>&& lbs,
                   std::vector<double>&& ubs,
                   std::vector<double>&& objective,
                   bitset&& integer,
                   const std::vector<size_t>& rowSize,
                   std::vector<std::string>&& varNames)
{
   assert(coefsT.size() == idxT.size());

   MIP<double> mip;

   size_t ncols = cols.size();
   size_t nrows = rows.size();

   // fill the column major matrix
   mip.constMatrixT.ncols = nrows;
   mip.constMatrixT.nrows = ncols;
   mip.constMatrixT.coefficients = coefsT;
   mip.constMatrixT.indices = idxT;
   mip.constMatrixT.rowStart = rstartT;

   mip.constMatrix = transpose(mip.constMatrixT, rowSize);

   // fill the rest
   mip.lhs = lhs;
   mip.rhs = rhs;

   mip.lb = lbs;
   mip.ub = ubs;

   mip.objective = objective;

   mip.integer = integer;

   mip.varNames = varNames;

   using RowInfo = std::pair<std::string, size_t>;
   std::vector<RowInfo> rowinfo;
   rowinfo.reserve(rows.size());

   for (auto& element : rows)
      rowinfo.emplace_back(element.first, element.second.second);

   std::sort(std::begin(rowinfo),
             std::end(rowinfo),
             [](const RowInfo& lhs, const RowInfo& rhs) {
                return lhs.second < rhs.second;
             });

   mip.consNames.reserve(rowinfo.size());

   for (auto& element : rowinfo)
      mip.consNames.push_back(std::move(element.first));

   return mip;
}

// read the mps, fill a transposed constraint matrix
// construct the sparse constraint matrix from the transposed
MIP<double>
mpsreader::parse(const std::string& filename)
{
   std::ifstream file(filename);
   boost::iostreams::filtering_istream in;

   if (!file.is_open())
      throw std::runtime_error("unable to open file: " + filename);

   if (boost::algorithm::ends_with(filename, ".gz"))
      in.push(boost::iostreams::gzip_decompressor());
   else if (boost::algorithm::ends_with(filename, ".bz2"))
      in.push(boost::iostreams::bzip2_decompressor());

   in.push(file);

   std::string name;

   Rows rows;
   Cols cols;

   // transposed sparse matrix
   std::vector<double> coefsT;
   std::vector<size_t> idxT;
   std::vector<size_t> rstatrtT;

   std::vector<double> rhs;
   std::vector<double> lhs;

   std::vector<double> lbs;
   std::vector<double> ubs;

   std::vector<std::string> varNames;
   std::vector<std::string> consNames;

   std::vector<double> objective;
   std::string objName;

   bitset integer;

   std::vector<size_t> rowSize;

   Section nextsection = NAME;

   Timer::time_point t0;
   Timer::time_point t1;

   try {
      while (nextsection != FAIL && nextsection != END) {
         switch (nextsection) {
            case NAME:
               nextsection = parseName(in, name);
               break;
            case ROWS:
               nextsection = parseRows(in, rows, objName);
               break;
            case COLUMNS:
               t0 = Timer::now();
               nextsection = parseColumns(in,
                                          rows,
                                          cols,
                                          coefsT,
                                          idxT,
                                          rstatrtT,
                                          objective,
                                          objName,
                                          integer,
                                          rowSize,
                                          varNames);
               t1 = Timer::now();
               fmt::print("Section COLUMNS parsed in {}\n",
                          Timer::seconds(t1, t0));
               break;
            case RHS:
               nextsection = parseRhs(in, rows, lhs, rhs);
               break;
            case BOUNDS:
               nextsection = parseBounds(in, cols, lbs, ubs, integer);
               break;
            case FAIL:
            case END:
            case NONE:
               assert(0);
         }
      }
   } catch (const std::exception& exp) {
      throw std::runtime_error("unable to parse file: " + filename + "\n" +
                               getErrorStr());
   }

   if (nextsection == FAIL)
      throw std::runtime_error("unable to parse file: " + filename + "\n" +
                               getErrorStr());

   return makeMip(rows,
                  cols,
                  std::move(coefsT),
                  std::move(idxT),
                  std::move(rstatrtT),
                  std::move(rhs),
                  std::move(lhs),
                  std::move(lbs),
                  std::move(ubs),
                  std::move(objective),
                  std::move(integer),
                  rowSize,
                  std::move(varNames));
}

std::string
mpsreader::getErrorStr()
{
   std::string errorstr = "Failed to parse mps file, error in section: ";
   switch (error_section) {
      case NAME:
         errorstr += "NAME";
         break;
      case ROWS:
         errorstr += "ROWS";
         break;
      case COLUMNS:
         errorstr += "COLUMNS";
         break;
      case RHS:
         errorstr += "RHS";
         break;
      case BOUNDS:
         errorstr += "BOUNDS";
         break;
      default:
         assert(0);
   }

   errorstr += "\n";

   return errorstr;
}

SparseMatrix<double>
mpsreader::transpose(const SparseMatrix<double>& matrix,
                     const std::vector<size_t>& rowSize)
{
   size_t nnz = matrix.coefficients.size();
   size_t ncols = matrix.nrows;
   size_t nrows = matrix.ncols;

   assert(std::accumulate(rowSize.begin(), rowSize.end(), 0) == nnz);

   SparseMatrix<double> transposed;
   transposed.nrows = nrows;
   transposed.ncols = ncols;
   transposed.coefficients.resize(nnz);
   transposed.indices.resize(nnz);

   auto& rowStart = transposed.rowStart;
   rowStart.reserve(nrows + 1);
   rowStart.push_back(0);

   std::vector<size_t> offset(nrows, 0);

   for (size_t i = 0; i < nrows; ++i)
      rowStart.push_back(rowSize[i] + rowStart[i]);

   for (size_t col = 0; col < ncols; ++col) {
      for (size_t rowid = matrix.rowStart[col];
           rowid < matrix.rowStart[col + 1];
           ++rowid) {
         size_t row = matrix.indices[rowid];
         double coef = matrix.coefficients[rowid];

         size_t rowstart = rowStart[row];
         transposed.coefficients[rowstart + offset[row]] = coef;
         transposed.indices[rowstart + offset[row]] = col;
         ++offset[row];
         assert(offset[row] <= rowSize[row]);
      }
   }

   assert(std::all_of(transposed.coefficients.begin(),
                      transposed.coefficients.end(),
                      [](double coef) { return coef != 0.0; }));

   return transposed;
}

mpsreader::Section mpsreader::error_section = mpsreader::NONE;
