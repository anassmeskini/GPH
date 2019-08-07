#include "GLPKSolver.h"
#include "core/Common.h"
#include <algorithm>
#include <numeric>

#ifdef GLPK_FOUND
GLPKSolver::GLPKSolver(const MIP& mip)
  : problem(nullptr)
  , ncols(mip.getNCols())
  , nrows(mip.getNRows())
{
   const auto& lb = mip.getLB();
   const auto& ub = mip.getUB();
   const auto& obj = mip.getObj();
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();

   problem = glp_create_prob();
   glp_set_prob_name(problem, "NONAME");
   glp_set_obj_dir(problem, GLP_MIN);
   glp_add_rows(problem, nrows);
   glp_add_cols(problem, ncols);

   constexpr double inf = std::numeric_limits<double>::infinity();

   std::vector<int> ind_buffer;
   ind_buffer.reserve(ncols);

   for (int col = 0; col < ncols; ++col)
   {
      // set objective
      glp_set_obj_coef(problem, col + 1, obj[col]);

      int boundtype;
      if (lb[col] == -inf && ub[col] == inf)
         boundtype = GLP_FR;
      else if (lb[col] == -inf)
         boundtype = GLP_UP;
      else if (ub[col] == inf)
         boundtype = GLP_LO;
      else if (lb[col] == ub[col])
         boundtype = GLP_FX;
      else
         boundtype = GLP_DB;

      glp_set_col_bnds(problem, col + 1, boundtype, lb[col], ub[col]);
   }

   for (int row = 0; row < nrows; ++row)
   {
      int constype;
      if (lhs[row] == -inf && rhs[row] == inf)
         constype = GLP_FR;
      else if (lhs[row] == -inf)
         constype = GLP_UP;
      else if (rhs[row] == inf)
         constype = GLP_LO;
      else if (lhs[row] == rhs[row])
         constype = GLP_FX;
      else
         constype = GLP_DB;

      glp_set_row_bnds(problem, row + 1, constype, lhs[row], rhs[row]);

      // set coefficients
      auto rowview = mip.getRow(row);

      std::transform(rowview.indices,
                     rowview.indices + rowview.size,
                     ind_buffer.begin(),
                     [&](int id) { return ++id; });

      glp_set_mat_row(problem,
                      row + 1,
                      rowview.size,
                      ind_buffer.data() - 1,
                      rowview.coefs - 1);
   }
}

LPResult
GLPKSolver::solve()
{
   LPResult result;
   int ret = glp_simplex(problem, nullptr);

   if (!ret)
   {
      int st = glp_get_status(problem);
      switch (st)
      {
         case GLP_OPT:
            result.status = LPResult::OPTIMAL;
            result.primalSolution.resize(ncols);
            for (int i = 0; i < ncols; ++i)
               result.primalSolution[i] = glp_get_col_prim(problem, i + 1);

            result.dualSolution.resize(nrows);
            for (int j = 0; j < nrows; ++j)
               result.dualSolution[j] = glp_get_row_prim(problem, j + 1);

            result.obj = glp_get_obj_val(problem);
            break;
         case GLP_INFEAS:
            result.status = LPResult::INFEASIBLE;
            break;
         case GLP_UNBND:
            result.status = LPResult::UNBOUNDED;
            break;
         case GLP_UNDEF:
            result.status = LPResult::OTHER;
            break;
      }
   }
   else
      result.status = LPResult::OTHER;
   // TODO

   return result;
}

GLPKSolver::GLPKSolver(const GLPKSolver& glpksolver)
  : problem(nullptr)
  , ncols(glpksolver.ncols)
  , nrows(glpksolver.nrows)
{
   problem = glp_create_prob();
   glp_copy_prob(problem, glpksolver.problem, GLP_OFF);
   assert(problem);
}

std::unique_ptr<LPSolver>
GLPKSolver::makeCopy() const
{
   return std::make_unique<GLPKSolver>(*this);
}

GLPKSolver::~GLPKSolver()
{
   glp_delete_prob(problem);
}

void
GLPKSolver::branch(int col, double val, Direction direction)
{
   glp_add_rows(problem, 1);
   ++nrows;

   const double one = 1.0;
   const int glpCol = col + 1;

   glp_set_mat_row(problem, nrows, 1, &glpCol - 1, &one - 1);

   if (direction == Direction::UP)
      glp_set_row_bnds(problem, nrows, GLP_LO, val, 0.0);
   else
      glp_set_row_bnds(problem, nrows, GLP_UP, 0.0, val);
}

void
GLPKSolver::changeBounds(int column, double lb, double ub)
{
   constexpr double inf = std::numeric_limits<double>::infinity();

   int boundtype;
   if (lb == -inf && ub == inf)
      boundtype = GLP_FR;
   else if (lb == -inf)
      boundtype = GLP_UP;
   else if (ub == inf)
      boundtype = GLP_LO;
   else if (lb == ub)
      boundtype = GLP_FX;
   else
      boundtype = GLP_DB;

   glp_set_col_bnds(problem, column + 1, boundtype, lb, ub);
}

void
GLPKSolver::changeBounds(const std::vector<double>& lb,
                         const std::vector<double>& ub)
{
   constexpr double inf = std::numeric_limits<double>::infinity();

   for (int col = 0; col < ncols; ++col)
   {
      int boundtype;
      if (lb[col] == -inf && ub[col] == inf)
         boundtype = GLP_FR;
      else if (lb[col] == -inf)
         boundtype = GLP_UP;
      else if (ub[col] == inf)
         boundtype = GLP_LO;
      else if (lb[col] == ub[col])
         boundtype = GLP_FX;
      else
         boundtype = GLP_DB;

      glp_set_col_bnds(problem, col + 1, boundtype, lb[col], ub[col]);
   }
}
#endif // GLPK_FOUND
