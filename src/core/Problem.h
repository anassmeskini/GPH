#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include "Common.h"
#include "MIP.h"

#include <vector>

class ProblemView
{
   public:
   ProblemView(MIP<double>&&, std::vector<double>&&);

   size_t getNCols() const;
   size_t getNRows() const;

   VectorView<double> getRow(size_t) const;
   VectorView<double> getCol(size_t) const;

   const std::vector<double>& getLB() const;
   const std::vector<double>& getUB() const;

   const std::vector<double>& getLHS() const;
   const std::vector<double>& getRHS() const;

   const std::vector<size_t>& getUpLocks() const;
   const std::vector<size_t>& getDownLocks() const;

   const bitset& integer() const;

   const std::vector<Activity>& getActivities() const;

   const std::vector<double>& getLPSol() const;

   const std::vector<double>& getSolActivity() const;

   private:
   const MIP<double> mip;

   size_t ncols;
   size_t nrows;

   std::vector<Activity> activities;

   std::vector<double> lpSol;
   std::vector<double> lpSolActivity;

   std::vector<size_t> fractional;

   std::vector<size_t> downLocks;
   std::vector<size_t> upLocks;
};

#endif
