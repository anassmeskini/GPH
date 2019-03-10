#ifndef MPS_READER_HPP
#define MPS_READER_HPP

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "MIP.h"
#include "unordered_map"

// TODO display meaningfull error messages
// TODO handle all section types
class mpsreader {
  public:
   static MIP<double> parse(const std::string& filename);

  private:
   enum Section : int {
      NAME,
      ROWS,
      COLUMNS,
      RHS,
      BOUNDS,
      END,
      FAIL,
      NONE,
   };

   enum ConsType {
      LESS,
      GREATER,
      EQUAL,
      OBJECTIVE,
   };

   static Section error_section;

   using Rows = std::unordered_map<std::string, std::pair<ConsType, size_t>>;

   using Cols = std::unordered_map<std::string, size_t>;

   static std::vector<std::string> split(const std::string& line);

   static std::vector<std::string> getNextElements(std::ifstream& file);

   static Section parseName(std::ifstream& file, std::string& name);

   static Section parseRows(std::ifstream& file, Rows& rows);

   static Section parseColumns(std::ifstream& file, const Rows& rows,
                               Cols& cols, std::vector<double>& coefs,
                               std::vector<size_t>& idxT,
                               std::vector<size_t>& rstart,
                               std::vector<double>& obj, bitset&,
                               std::vector<std::string>&);

   static Section parseRhs(std::ifstream& file, const Rows& rows,
                           std::vector<double>& lhs, std::vector<double>& rhs);

   static Section parseBounds(std::ifstream& file, const Cols& cols,
                              std::vector<double>& lbs,
                              std::vector<double>& ubs);

   static SparseMatrix<double> compress(const std::vector<double>&, size_t);

   static MIP<double> makeMip(
       const Rows& rows, const Cols& cols, std::vector<double>&& coefsT,
       std::vector<size_t>&& idxT, std::vector<size_t>&& rstartT,
       std::vector<double>&& rhs, std::vector<double>&& lhs,
       std::vector<double>&& lbs, std::vector<double>&& ubs,
       std::vector<double>&& obj, bitset&& integer, std::vector<std::string>&&);

   static std::string getErrorStr();
};

std::vector<std::string> mpsreader::split(const std::string& line) {
   constexpr char space = ' ';
   constexpr char tab = '\t';

   size_t beginId;
   std::vector<std::string> tokens;

   for (size_t id = 0; id < line.size(); ++id) {
      if (line[id] == space || line[id] == tab) continue;

      beginId = id;
      while (true) {
         ++id;
         assert(id <= line.size());

         if (id >= line.size() || line[id] == space || line[id] == tab) {
            tokens.emplace_back(line.begin() + beginId, line.begin() + id);
            break;
         }
      }
   }

   return tokens;
}

std::vector<std::string> mpsreader::getNextElements(std::ifstream& file) {
   std::string line;
   while (std::getline(file, line) && (line.empty() || line[0] == '*')) {
   }

   return split(line);
}

mpsreader::Section mpsreader::parseName(std::ifstream& file,
                                        std::string& name) {
   auto tokens = getNextElements(file);
   error_section = NAME;

   if (tokens.size() != 2 || tokens[0] != "NAME") return FAIL;

   name = std::move(tokens[1]);

   tokens = getNextElements(file);

   if (tokens.size() != 1 || tokens[0] != "ROWS") return FAIL;

   error_section = NONE;
   return ROWS;
}

mpsreader::Section mpsreader::parseRows(std::ifstream& file, Rows& rows) {
   error_section = ROWS;
   std::string line;
   std::vector<std::string> tokens;

   size_t rowcounter = 0;
   while (std::getline(file, line)) {
      if (line.empty() || line[line.find_last_not_of(" ")] == '*') continue;

      tokens = split(line);
      if (tokens.size() == 1) break;
      if (tokens.size() != 2 || tokens[0].size() != 1) return FAIL;

      ConsType type;
      switch (tokens[0][0]) {
         case 'N':
            type = OBJECTIVE;
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
      }

      auto pair =
          rows.emplace(std::move(tokens[1]), std::make_pair(type, rowcounter));
      if (!pair.second) return FAIL;

      if (type != OBJECTIVE) ++rowcounter;
   }

   assert(tokens.size() == 1);
   if (tokens[0] != "COLUMNS") return FAIL;

   error_section = NONE;
   return COLUMNS;
}

mpsreader::Section mpsreader::parseColumns(
    std::ifstream& file, const Rows& rows, Cols& cols,
    std::vector<double>& coefs, std::vector<size_t>& idxT,
    std::vector<size_t>& rstart, std::vector<double>& objective,
    bitset& integer, std::vector<std::string>& varNames) {
   error_section = COLUMNS;
   std::string line;
   std::vector<std::string> tokens;

   int colId = -1;
   std::string prevCol("");
   std::set<std::string> colset;

   bool integerSection = false;

   while (std::getline(file, line)) {
      assert(coefs.size() == idxT.size());

      if (line.empty() || line[line.find_last_not_of(" ")] == '*') continue;

      tokens = split(line);
      if (tokens.size() == 1) break;

      // check if it's the start or end of an integer section
      if (tokens.size() == 3) {
         if (tokens[1] == "'MARKER'") {
            if (tokens[2] == "'INTORG'") {
               if (integerSection) return FAIL;
               integerSection = true;
               continue;
            } else if (tokens[2] == "'INTEND'") {
               if (!integerSection) return FAIL;
               integerSection = false;
               continue;
            } else
               return FAIL;
         }
      } else if (!(tokens.size() % 2))
         return FAIL;

      auto& curCol = tokens[0];

      // if it's a new column
      if (curCol != prevCol) {
         if (colset.count(curCol)) return FAIL;
         ++colId;
         colset.insert(curCol);
         cols.emplace(tokens[0], colId);

         rstart.push_back(coefs.size());
         integer.push_back(integerSection);

         varNames.push_back(tokens[0]);
         assert(varNames.size() == colId + 1);

         prevCol = curCol;
      }

      for (size_t i = 1; i < tokens.size(); i += 2) {
         auto& rowname = tokens[i];
         double coef = std::stod(tokens[i + 1]);

         auto iter = rows.find(rowname);
         // row not declared in the ROWS section
         if (iter == rows.end()) return FAIL;

         // if the entry is for the objective
         if (iter->second.first == OBJECTIVE) {
            // TODO
            assert(colId >= 0);
            // account for the columns that didn't have an objective entry
            // colId = 2 size = 0
            // iter 1: size = 1
            // iter 2: size = 2
            while (objective.size() < static_cast<size_t>(colId))
               objective.push_back(0.0);

            objective.push_back(coef);
            assert(colId && objective.size() == static_cast<size_t>(colId) + 1);
         } else {
            auto rowId = iter->second.second;
            coefs.push_back(coef);
            idxT.push_back(rowId);
         }
      }
   }

   while (objective.size() < static_cast<size_t>(colId) + 1)
      objective.push_back(0.0);

   rstart.push_back(coefs.size());

   assert(tokens.size() == 1);
   if (tokens[0] != "RHS") return FAIL;

   assert(objective.size() == cols.size());

   error_section = NONE;
   return RHS;
}

mpsreader::Section mpsreader::parseRhs(std::ifstream& file, const Rows& rows,
                                       std::vector<double>& lhs,
                                       std::vector<double>& rhs) {
   error_section = RHS;
   std::string line;
   std::vector<std::string> tokens;

   std::string prevCol("");
   std::set<std::string> colset;

   // TODO handle case where the objective is missing
   lhs = std::vector<double>(rows.size() - 1);

   rhs = std::vector<double>(rows.size() - 1);

   // set default bounds
   for (auto row : rows) {
      auto constype = row.second.first;
      auto id = row.second.second;
      double inf = std::numeric_limits<double>::infinity();

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
      if (line.empty() || line[line.find_last_not_of(" ")] == '*') continue;

      tokens = split(line);
      if (tokens.size() == 1) break;
      if (!(tokens.size() % 2)) return FAIL;

      for (size_t i = 1; i < tokens.size(); i += 2) {
         auto& rowname = tokens[i];
         double side = std::stof(tokens[i + 1]);

         auto iter = rows.find(rowname);
         if (iter == rows.end()) return FAIL;

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
   if (tokens[0] != "BOUNDS") return FAIL;

   error_section = NONE;
   return BOUNDS;
}

// TODO handle bounds properly
mpsreader::Section mpsreader::parseBounds(std::ifstream& file, const Cols& cols,
                                          std::vector<double>& lbs,
                                          std::vector<double>& ubs) {
   error_section = BOUNDS;
   std::string line;
   std::vector<std::string> tokens;

   lbs = std::vector<double>(cols.size(), 0);
   ubs = std::vector<double>(cols.size(),
                             std::numeric_limits<double>::infinity());

   bitset lb_changed(cols.size(), false);

   while (std::getline(file, line)) {
      if (line.empty() || line[line.find_last_not_of(" ")] == '*') continue;

      tokens = split(line);
      if (tokens.size() == 1) break;
      if (tokens.size() != 4 && tokens.size() != 3) return FAIL;

      auto& colname = tokens[2];
      double bound = std::stof(tokens[3]);

      auto iter = cols.find(colname);
      if (iter == cols.end()) return FAIL;

      size_t colid = iter->second;

      if (tokens.size() == 4) {
         if (tokens[0] == "UP") {
            ubs[colid] = bound;
            if (bound < 0.0 && !lb_changed[colid])
               lbs[colid] = -std::numeric_limits<double>::infinity();
         } else if (tokens[0] == "LO") {
            lbs[colid] = bound;
            lb_changed[colid] = true;
         } else if (tokens[0] == "FX") {
            lbs[colid] = bound;
            ubs[colid] = bound;
         } else
            return FAIL;

      } else if (tokens.size() == 3) {
         if (tokens[0] == "FR") {
            lbs[colid] = -std::numeric_limits<double>::infinity();
            ubs[colid] = std::numeric_limits<double>::infinity();
         } else
            return FAIL;
      }
   }

   assert(tokens.size() == 1);
   if (tokens[0] != "ENDATA") return FAIL;

   error_section = NONE;
   return END;
}

SparseMatrix<double> mpsreader::compress(const std::vector<double>& denseCoefs,
                                         size_t ncols) {
   SparseMatrix<double> matrix;

   assert(!(denseCoefs.size() % ncols));

   size_t nrows = static_cast<size_t>(denseCoefs.size() / ncols);
   size_t nnz = 0;

   for (size_t row = 0; row < nrows; ++row) {
      matrix.rowStart.push_back(nnz);
      for (size_t col = 0; col < ncols; ++col) {
         double coef = denseCoefs[ncols * row + col];
         if (coef != 0.0) {
            matrix.coefficients.push_back(coef);
            matrix.indices.push_back(col);
            ++nnz;
         }
      }
   }

   matrix.rowStart.push_back(nnz);

   matrix.ncols = ncols;
   matrix.nrows = nrows;

   return matrix;
}

MIP<double> mpsreader::makeMip(
    const Rows& rows, const Cols& cols, std::vector<double>&& coefsT,
    std::vector<size_t>&& idxT, std::vector<size_t>&& rstartT,
    std::vector<double>&& rhs, std::vector<double>&& lhs,
    std::vector<double>&& lbs, std::vector<double>&& ubs,
    std::vector<double>&& objective, bitset&& integer,
    std::vector<std::string>&& varNames) {
   assert(coefsT.size() == idxT.size());

   MIP<double> mip;

   // TODO
   size_t ncols = cols.size();
   size_t nrows = rows.size() - 1;

   // fill the column major matrix
   mip.constMatrixT.ncols = nrows;
   mip.constMatrixT.nrows = ncols;
   mip.constMatrixT.coefficients = coefsT;
   mip.constMatrixT.indices = idxT;
   mip.constMatrixT.rowStart = rstartT;

   std::vector<double> rowMajorCoefs(ncols * nrows, 0.0);

   for (size_t colid = 0; colid < rstartT.size() - 1; ++colid) {
      for (size_t rowid = rstartT[colid]; rowid < rstartT[colid + 1]; ++rowid) {
         size_t slot = ncols * idxT[rowid] + colid;
         assert(slot < rowMajorCoefs.size());

         rowMajorCoefs[slot] = coefsT[rowid];
      }
   }

   // get the row-major matrix
   mip.constMatrix = compress(rowMajorCoefs, ncols);

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

   std::sort(std::begin(rowinfo), std::end(rowinfo),
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
MIP<double> mpsreader::parse(const std::string& filename) {
   std::ifstream file(filename);

   if (!file.is_open())
      throw std::runtime_error("unable to open file: " + filename);

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

   bitset integer;

   Section nextsection = NAME;
   while (nextsection != FAIL && nextsection != END) {
      switch (nextsection) {
         case NAME:
            nextsection = parseName(file, name);
            break;
         case ROWS:
            nextsection = parseRows(file, rows);
            break;
         case COLUMNS:
            nextsection = parseColumns(file, rows, cols, coefsT, idxT, rstatrtT,
                                       objective, integer, varNames);
            break;
         case RHS:
            nextsection = parseRhs(file, rows, lhs, rhs);
            break;
         case BOUNDS:
            nextsection = parseBounds(file, cols, lbs, ubs);
            break;
         case FAIL:
         case END:
         case NONE:
            assert(0);
      }
   }

   if (nextsection == FAIL)
      throw std::runtime_error("unable to parse file: " + filename + "\n" +
                               getErrorStr());

   return makeMip(rows, cols, std::move(coefsT), std::move(idxT),
                  std::move(rstatrtT), std::move(rhs), std::move(lhs),
                  std::move(lbs), std::move(ubs), std::move(objective),
                  std::move(integer), std::move(varNames));
}

std::string mpsreader::getErrorStr() {
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

mpsreader::Section mpsreader::error_section = mpsreader::NONE;

#endif
