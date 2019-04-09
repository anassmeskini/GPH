#ifndef _COMMON_HPP_
#define _COMMON_HPP_

#include <cassert>
#include <vector>

template<typename REAL>
class MIP;

template<typename REAL>
bool
checkFeasibility(const MIP<REAL>& mip,
                 const std::vector<REAL>& sol,
                 REAL boundtol,
                 REAL constol)
{
  const auto& ub = mip.getUB();
  const auto& lb = mip.getLB();

  assert(sol.size() == mip.getNCols());

  for (size_t colid = 0; colid < mip.getNCols(); ++colid) {
    if (sol[colid] > ub[colid] + boundtol || sol[colid] < lb[colid] - boundtol)
      return false;
  }

  const auto& lhs = mip.getLHS();
  const auto& rhs = mip.getRHS();

  for (size_t rowid = 0; rowid < mip.getNRows(); ++rowid) {
    const auto& row = mip.getRow(rowid);
    REAL activity{ 0.0 };

    for (size_t id = 0; id < row.size; ++id)
      activity += row.coefs[id] * sol[row.indices[id]];

    if (activity > rhs[rowid] + constol || activity < lhs[rowid] - constol)
      return false;
  }

  // TODO integrality check

  return true;
}

template<typename REAL>
bool
checkExactFeasibility(const MIP<REAL>& mip, const std::vector<REAL>& sol)
{
  const auto& ub = mip.getUB();
  const auto& lb = mip.getLB();

  assert(sol.size() == mip.getNCols());

  for (size_t colid = 0; colid < mip.getNCols(); ++colid) {
    if (sol[colid] > ub[colid] || sol[colid] < lb[colid])
      return false;
  }

  const auto& lhs = mip.getLHS();
  const auto& rhs = mip.getRHS();

  for (size_t rowid = 0; rowid < mip.getNRows(); ++rowid) {
    auto& row = mip.getRow(rowid);
    REAL activity{ 0.0 };

    for (size_t id = 0; id < row.size(); ++id)
      activity += row.coefs[id] * sol[row.indices[id]];

    if (activity > rhs[rowid] || activity < lhs[rowid])
      return false;
  }

  // TODO integrality check

  return true;
}

#endif
