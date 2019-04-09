#ifndef MPS_READER_HPP
#define MPS_READER_HPP

#include "MIP.h"

#include "ska/Hash.hpp"
#include <algorithm>
#include <array>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>

struct Split
{
  std::array<std::pair<int, int>, 5> words;
  size_t size_ = 0;

  int size(size_t i) { return words[i].second - words[i].first; }
  size_t size() { return size_; }
};

// TODO handle RANGES section
class mpsreader
{
public:
  static MIP<double> parse(const std::string& filename);

private:
  enum Section : int
  {
    NAME,
    ROWS,
    COLUMNS,
    RHS,
    BOUNDS,
    END,
    FAIL,
    NONE,
  };

  enum ConsType
  {
    LESS,
    GREATER,
    EQUAL,
    OBJECTIVE,
  };

  static Section error_section;

  // stores info about the rows when reading the ROWS section
  // name -> <contraint type, row id>
  using Rows = HashMap<std::string, std::pair<ConsType, size_t>>;

  // name -> column id
  using Cols = HashMap<std::string, size_t>;

  static Section parseName(boost::iostreams::filtering_istream& file,
                           std::string& name);

  static Section parseRows(boost::iostreams::filtering_istream& file,
                           Rows& rows,
                           std::string& objName);

  static Section parseColumns(boost::iostreams::filtering_istream& file,
                              const Rows& rows,
                              Cols& cols,
                              std::vector<double>& coefs,
                              std::vector<size_t>& idxT,
                              std::vector<size_t>& rstart,
                              std::vector<double>& obj,
                              const std::string& objName,
                              bitset&,
                              std::vector<size_t>&,
                              std::vector<std::string>&);

  static Section parseRhs(boost::iostreams::filtering_istream& file,
                          const Rows& rows,
                          std::vector<double>& lhs,
                          std::vector<double>& rhs);

  static Section parseBounds(boost::iostreams::filtering_istream& file,
                             const Cols& cols,
                             std::vector<double>& lbs,
                             std::vector<double>& ubs,
                             bitset& integer);

  static SparseMatrix<double> compress(const std::vector<double>&, size_t);

  static SparseMatrix<double> transpose(const SparseMatrix<double>&,
                                        const std::vector<size_t>&);

  static MIP<double> makeMip(const Rows& rows,
                             const Cols& cols,
                             std::vector<double>&& coefsT,
                             std::vector<size_t>&& idxT,
                             std::vector<size_t>&& rstartT,
                             std::vector<double>&& rhs,
                             std::vector<double>&& lhs,
                             std::vector<double>&& lbs,
                             std::vector<double>&& ubs,
                             std::vector<double>&& obj,
                             bitset&& integer,
                             const std::vector<size_t>& rowSize,
                             std::vector<std::string>&&);

  static std::string getErrorStr();
};

#endif
