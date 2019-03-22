#include "interface/GLPKSolver.h"
#include <algorithm>
#include <numeric>

#ifdef GLPK_FOUND
GLPKSolver::GLPKSolver(const MIP<double>& mip)
    : LPSolver(mip), problem(nullptr) {
   size_t ncols = mip.getNCols();
   size_t nrows = mip.getNRows();

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

   //

   for (size_t col = 0; col < ncols; ++col) {
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

   for (size_t row = 0; row < nrows; ++row) {
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
      std::transform(rowview.indices, rowview.indices + rowview.size,
                     intbuffer.begin(),
                     [&](size_t id) { return static_cast<int>(id) + 1; });

      glp_set_mat_row(problem, row + 1, rowview.size, intbuffer.data() - 1,
                      rowview.coefs - 1);
   }
}

LPResult GLPKSolver::solve() {
   LPResult result;
   int ret = glp_simplex(problem, NULL);

   size_t ncols = mip.getNCols();
   size_t nrows = mip.getNRows();

   if (!ret) {
      int st = glp_get_status(problem);
      switch (st) {
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
   } else
      result.status = LPResult::OTHER;

   return result;
}

GLPKSolver::~GLPKSolver() {}

#endif  // GLPK_FOUND
