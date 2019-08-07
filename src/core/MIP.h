#ifndef LP_HPP
#define LP_HPP

#include "SparseMatrix.h"
#include "dynamic_bitset/dynamic_bitset.hpp"
#include "ska/Hash.hpp"

#include <string>
#include <vector>

enum ConsType
{
   LESS,
   GREATER,
   EQUAL,
   OBJECTIVE,
};

// name -> <contraint type, row id>
using Rows = HashMap<std::string, std::pair<ConsType, int>>;

// name -> column id
using Cols = HashMap<std::string, int>;

struct VectorView
{
   VectorView(const double* _array, const int* _indices, int _size)
     : coefs(_array)
     , indices(_indices)
     , size(_size)
   {
   }

   VectorView(const VectorView&) = default;

   const double* coefs;
   const int* indices;
   int size;
};

struct Statistics
{
   int ncols = -1;
   int nrows = -1;

   int nbin = 0;
   int nint = 0;
   int ncont = 0;

   int nequality = 0;
   int ninequality = 0;

   int avgRowSupport = 0;
   int avgColSupport = 0;
   int maxRowSupport = 0;
   int maxColSupport = 0;
   int minRowSupport = std::numeric_limits<int>::max();
   int minColSupport = std::numeric_limits<int>::max();

   int maxLocks = 0;
   int minLocks = std::numeric_limits<int>::max();

   int nnzmat = 0;
   int nnzobj = 0;
};

void printStats(Statistics);

class MIP
{
   public:
   MIP() = default;

   MIP(const Rows& rows,
       const Cols& cols,
       std::vector<double>&& coefsT,
       std::vector<int>&& idxT,
       std::vector<int>&& rstartT,
       std::vector<double>&& rhs,
       std::vector<double>&& lhs,
       std::vector<double>&& lbs,
       std::vector<double>&& ubs,
       std::vector<double>&& obj,
       dynamic_bitset<>&& integer,
       std::vector<int>& rowSize,
       std::vector<std::string>&& colNames);

   MIP(MIP&&) = default;

   MIP& operator=(MIP&&) noexcept;

   int getNCols() const { return constMatrixT.nrows; }

   int getNRows() const { return constMatrix.nrows; }

   const std::vector<double>& getObj() const { return objective; }

   const std::vector<double>& getLB() const { return lb; }

   const std::vector<double>& getUB() const { return ub; }

   const std::vector<double>& getLHS() const { return lhs; }

   const std::vector<double>& getRHS() const { return rhs; }

   const dynamic_bitset<>& getInteger() const { return integer; }

   const std::vector<std::string>& getVarNames() const { return varNames; }

   const std::vector<std::string>& getConsNames() const { return consNames; }

   VectorView getRow(int row) const noexcept;

   VectorView getCol(int col) const noexcept;

   int getRowSize(int row) const noexcept;

   int getColSize(int row) const noexcept;

   const std::vector<int>& getUpLocks() const { return upLocks; }

   const std::vector<int>& getDownLocks() const { return downLocks; }

   Statistics getStatistics() const { return stats; }

   const std::vector<int>& getBinaryVars() const { return binary; }

   const std::vector<int>& getIntVars() const { return generalInt; }

   const std::vector<int>& getContVars() const { return continuous; }

   private:
   static SparseMatrix transpose(const SparseMatrix&, const std::vector<int>&);

   // min {obj*x}
   std::vector<double> objective;
   double objoffset = 0.0;

   // lhs <= Ax <= rhs
   std::vector<double> lhs;
   std::vector<double> rhs;

   // lb <= x <= ub
   std::vector<double> lb;
   std::vector<double> ub;

   std::vector<std::string> varNames;
   std::vector<std::string> consNames;

   // row-major sparse
   SparseMatrix constMatrix;

   // column-major sparse
   SparseMatrix constMatrixT;

   dynamic_bitset<> integer;

   std::vector<int> binary;
   std::vector<int> generalInt;
   std::vector<int> continuous;

   std::vector<int> downLocks;
   std::vector<int> upLocks;

   Statistics stats;
};

#endif
