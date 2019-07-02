#include "GLPKSolver.h"
#include <algorithm>
#include <numeric>

#ifdef GLPK_FOUND
GLPKSolver::GLPKSolver(const MIP<double>& mip)
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

   for (size_t col = 0; col < ncols; ++col)
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

   for (size_t row = 0; row < nrows; ++row)
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
      // TODO remove this when int is used instead of size_t
      std::vector<int> intbuffer(rowview.size);
      std::transform(rowview.indices,
                     rowview.indices + rowview.size,
                     intbuffer.begin(),
                     [&](size_t id) { return static_cast<int>(id) + 1; });

      glp_set_mat_row(problem,
                      row + 1,
                      rowview.size,
                      intbuffer.data() - 1,
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
            for (size_t i = 0; i < ncols; ++i)
               result.primalSolution[i] = glp_get_col_prim(problem, i + 1);

            result.dualSolution.resize(nrows);
            for (size_t j = 0; j < nrows; ++j)
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

LPResult
GLPKSolver::solve(LPAlgorithm alg)
{
   glp_smcp param;
   glp_init_smcp(&param);

   param.msg_lev = GLP_MSG_ERR;

   bool barrier = false;
   switch (alg)
   {
      case LPAlgorithm::PRIMAL_SPX:
         param.meth = GLP_PRIMAL;
         break;
      case LPAlgorithm::AUTO:
      case LPAlgorithm::DUAL_SPX:
         param.meth = GLP_DUAL;
         break;
      case LPAlgorithm::BARRIER:
         barrier = true;
         break;
   }

   int ret;
   if (!barrier)
      ret = glp_simplex(problem, &param);
   else
      ret = glp_interior(problem, nullptr);
   // TODO

   LPResult result;
   if (!ret)
   {
      int st = glp_get_status(problem);
      switch (st)
      {
         case GLP_OPT:
            result.status = LPResult::OPTIMAL;
            result.primalSolution.resize(ncols);
            for (size_t i = 0; i < ncols; ++i)
               result.primalSolution[i] = glp_get_col_prim(problem, i + 1);

            result.dualSolution.resize(nrows);
            for (size_t j = 0; j < nrows; ++j)
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

std::unique_ptr<LPSolver<double>>
GLPKSolver::makeCopy() const
{
   return std::make_unique<GLPKSolver>(*this);
}

GLPKSolver::~GLPKSolver()
{
   glp_delete_prob(problem);
}

void GLPKSolver::branch(int col, double val, Direction direction){
   glp_add_rows(problem, 1);
   ++nrows;

   const double one = 1.0;
   const int glpCol = col + 1;

   glp_set_mat_row(problem,
                  nrows,
                  1,
                  &glpCol - 1,
                  &one -1);

   if(direction == Direction::UP)
       glp_set_row_bnds(problem, nrows, GLP_LO, val, 0.0);
   else
       glp_set_row_bnds(problem, nrows, GLP_UP, 0.0, val);
}

#endif // GLPK_FOUND
