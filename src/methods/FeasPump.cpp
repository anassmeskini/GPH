#include "FeasPump.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "core/Propagation.h"
#include "io/Message.h"
#include "io/SOLFormat.h"

#include <cmath>

void
FeasPump::search(const MIP& mip, const std::vector<double>& lb,
                 const std::vector<double>& ub,
                 const std::vector<Activity>& activities,
                 const LPResult& result, const std::vector<double>&,
                 const std::vector<int>&,
                 std::shared_ptr<const LPSolver> solver, TimeLimit tlimit,
                 SolutionPool& pool)
{
   int ncols = mip.getNCols();
   int nrows = mip.getNRows();
   auto st = mip.getStats();
   const auto& objective = mip.getObj();

   std::unique_ptr<LPSolver> localsolver = solver->clone();

   if (st.nnzmat >= 1.5e5)
   {
      // TODO : disable propagation in this case
      return;
   }

   // the factor in the weight of the old objective in
   // the new objective
   // modified objective:
   // (1 - aplha) * d(x, x_int) + objfactor * aplha * c * x
   double obj_norm_square = 0.0;
   for (int col = 0; col < ncols; ++col)
      obj_norm_square += objective[col] * objective[col];
   double obj_factor = std::sqrt((st.nbin + st.nint) / obj_norm_square);

   Message::debug_details("FeasPump: objective factor {:0.4f}",
                          obj_factor);

   // copies
   std::vector<double> locallb(ncols);
   std::vector<double> localub(ncols);
   std::vector<double> rounded_sol(ncols);
   std::vector<Activity> local_activities(nrows);

   // the current lp sol
   auto lp_sol = result.primalSolution;

   // the history of visited integer solutions
   std::array<std::vector<double>, cycle_length> int_sols;

   int iter = 0;
   double total_frac = -1.0;
   double last_frac = static_cast<double>(ncols);
   double alpha = 1.0;
   int stall_iter = 0;
   int restarts = 0;

   int nlp_iter = 0;
   bool limit_reached = false;
   do
   {
      locallb = lb;
      localub = ub;
      local_activities = activities;

      // fix and propagate
      bool propagation_feas = true;
      for (int col = 0; col < st.nbin + st.nint; ++col)
      {
         if (Num::isFeasEQ(locallb[col], localub[col]))
            continue;

         double oldlb = locallb[col];
         double oldub = localub[col];

         double int_val = Num::round(lp_sol[col]);
         if (Num::isFeasGE(int_val, locallb[col]) &&
             Num::isFeasLE(int_val, localub[col]))
         {
            locallb[col] = int_val;
            localub[col] = int_val;
         }
         else if (objective[col] > 0.0)
            locallb[col] = localub[col];
         else
            localub[col] = locallb[col];

         if (propagation_feas)
            propagation_feas =
                propagate(mip, locallb, localub, local_activities, col,
                          oldlb, oldub);
      }

      // the rounded solution is feasible
      if (propagation_feas)
      {
         Message::debug("feasPump: Rounded solution is feasible, iter {} "
                        "restarts {} stalls {}",
                        iter, restarts, stall_iter);

         assert(!locallb.empty());
         std::vector<double> sol = std::move(locallb);

         if (st.ncont > 0)
         {
            localsolver.reset();
            localsolver = solver->clone();
            for (int col = 0; col < st.nbin + st.nint; ++col)
            {
               localsolver->changeBounds(col, sol[col], sol[col]);
            }

            auto local_result = localsolver->solve(Algorithm::DUAL);

            if (local_result.status != LPResult::OPTIMAL)
            {
               Message::debug("FeasPump: lp is infeasible");
               break;
            }

            sol = std::move(local_result.primalSolution);
         }

         assert(!sol.empty());
         double cost = 0;
         for (int col = 0; col < ncols; ++col)
            cost += objective[col] * sol[col];
         pool.add(std::move(sol), cost);
         break;
      }

      assert(!locallb.empty());
      rounded_sol = std::move(locallb);

      // random flips and cycle detection
      if (iter > min_iter_pert && iter % pert_freq == 0)
      {
         Message::debug(
             "feasPump: making periodic perturbation, iter {} stalls {}",
             iter, stall_iter);

         make_rand_perturbation(rounded_sol, lp_sol, st.nbin + st.nint,
                                lb);
      }
      else if (alpha < alpha_cycle_detection_threshold)
      {
         assert(iter < 3 || std::all_of(int_sols.begin(), int_sols.end(),
                                        [](const auto& vec) -> bool {
                                           return !vec.empty();
                                        }));

         int cycle =
             is_already_visited<cycle_length>(int_sols, rounded_sol);

         if (cycle)
         {
            Message::debug_details("feasPump: found cycle of length {}",
                                   cycle);
            if (++restarts >= max_restarts)
            {
               Message::debug_details("feasPump: reached maximum restarts "
                                      "after {} iterations",
                                      iter);
               break;
            }

            if (cycle == 1)
               handle_one_cycle(rounded_sol, lp_sol, lb, ub,
                                st.nbin + st.nint);
            else
               make_rand_perturbation(rounded_sol, lp_sol,
                                      st.nbin + st.nint, lb);
         }
      }

      // set up the current iteration's objective : integers
      int nlbvar = 0;
      int nubvar = 0;
      for (int col = 0; col < st.nbin + st.nint; ++col)
      {
         if (Num::isFeasEQ(rounded_sol[col], lb[col]))
         {
            ++nlbvar;
            localsolver->changeObjective(
                col, 1.0 + alpha * (obj_factor * objective[col] - 1.0));
         }
         else if (Num::isFeasEQ(rounded_sol[col], ub[col]))
         {
            ++nubvar;
            localsolver->changeObjective(
                col, -1.0 + alpha * (obj_factor * objective[col] + 1.0));
         }
         else
            localsolver->changeObjective(col, obj_factor * alpha *
                                                  objective[col]);
      }

      // set up the current iteration's objective : continuous
      for (int col = st.nbin + st.nint; col < st.ncols; ++col)
      {
         localsolver->changeObjective(col,
                                      obj_factor * alpha * objective[col]);
      }

      // TODO handle case where all integer variables are not binary
      // and no variable is rounded to its bound
      assert(nlbvar || nubvar);

      // solve the lp
      auto local_result = localsolver->solve(Algorithm::PRIMAL);
      assert(local_result.status == LPResult::OPTIMAL);
      assert(!local_result.primalSolution.empty());
      lp_sol = std::move(local_result.primalSolution);

      // compute the fractionality
      last_frac = total_frac;
      total_frac = get_frac(lp_sol, st.nbin + st.nint);

      // check if the lp solution is integer
      if (total_frac < zero_frac)
      {
         Message::debug("feasPump: lp solution is feasible");

         assert(!lp_sol.empty());
         // roundFeasIntegers(lp_sol, integer);

         double cost = 0;
         for (int col = 0; col < ncols; ++col)
            cost += objective[col] * lp_sol[col];

         // TODO add feasibility asserts
         pool.add(std::move(lp_sol), cost);
         break;
      }

      // check for stalling
      if (iter > 5 && total_frac > last_frac - min_frac_improv)
         ++stall_iter;
      else
         stall_iter = 0;

#ifndef NDEBUG
      double l1dist = 0.0;
      for (int col = 0; col < ncols; ++col)
         l1dist += std::fabs(rounded_sol[col] - lp_sol[col]);
      Message::debug_details(
          "iter {}, restarts {}, stalls {}, alpha {:0.2f}"
          ", frac {:0.2f}, l1dist {:0.2f}, nvarlb {}, nvarub {}, "
          "spx_pivots {}",
          iter, restarts, stall_iter, alpha, total_frac, l1dist, nlbvar,
          nubvar, local_result.niter);
#endif

      // update the solution history
      for (int i = cycle_length - 1; i > 0; --i)
         int_sols[i] = std::move(int_sols[i - 1]);
      int_sols[0] = std::move(rounded_sol);

      assert(iter < 3 || std::all_of(int_sols.begin(), int_sols.end(),
                                     [](const auto& vec) -> bool {
                                        return !vec.empty();
                                     }));

      // update alpha
      alpha *= alpha_update_factor;

      nlp_iter = local_result.niter;

      limit_reached = (nlp_iter > max_lp_iter_ratio * result.niter &&
                       iter % pert_freq >= free_lp_iter_rounds) ||
                      (stall_iter >= max_stall_iter);

      if (limit_reached)
      {
         Message::debug("FeasPump: limit reached nlp_iter {}, nstalls {}",
                        nlp_iter, stall_iter);
         break;
      }

   } while (++iter < max_iter && !tlimit.reached(Timer::now()));
}

double
FeasPump::get_frac(const std::vector<double>& lp_sol, int ninteger)
{
   double frac = 0.0;

   for (int col = 0; col < ninteger; ++col)
   {
      double floor_dist = lp_sol[col] - Num::floor(lp_sol[col]);
      frac += std::min(floor_dist, 1.0 - floor_dist);
   }

   return frac;
}

void
FeasPump::make_rand_perturbation(std::vector<double>& rounded_sol,
                                 const std::vector<double>& lp_sol,
                                 int ninteger,
                                 const std::vector<double>& lb)
{
   for (int col = 0; col < ninteger; ++col)
   {
      double rnd = static_cast<double>(rand()) / RAND_MAX;
      double floor_dist = lp_sol[col] - Num::floor(lp_sol[col]);
      double frac = std::min(floor_dist, 1.0 - floor_dist);

      if (frac + std::max(rnd, 0.3) > 0.7)
      {
         if (Num::isFeasEQ(rounded_sol[col], lb[col]))
            rounded_sol[col] += 1.0;
         else
            rounded_sol[col] -= 1.0;
      }
   }
}

void
FeasPump::handle_one_cycle(std::vector<double>& rounded_sol,
                           const std::vector<double>& lp_sol,
                           const std::vector<double>& lb,
                           const std::vector<double>& ub, int ninteger)
{
   int ncols = rounded_sol.size();
   std::vector<int> perm(ncols);
   std::vector<int> score(ncols);

   // compute score for sorting
   for (int col = 0; col < ncols; ++col)
   {
      perm[col] = col;
      if (col < ninteger)
      {
         double is_bin = lb[col] == 0.0 && ub[col] == 1.0;
         double floor_dist = lp_sol[col] - Num::floor(lp_sol[col]);
         double frac = std::min(floor_dist, 1.0 - floor_dist);
         score[col] = -(10.0 * is_bin + frac);
      }
      else
         score[col] = 1.0;
   }

   // sort
   std::sort(std::begin(perm), std::end(perm),
             [&score](int left, int right) -> bool {
                return score[left] < score[right];
             });

   int real_min_flips = std::min(min_flips, ninteger);
   int real_max_flips = std::min(max_flips, ninteger);

   int nflips =
       rand() % (real_max_flips - real_min_flips) + real_min_flips;

   // flip the most fractional variables
   for (int i = 0; i < nflips; ++i)
   {
      int col = perm[i];

      if (Num::isFeasEQ(lb[col], rounded_sol[col]))
         rounded_sol[col] += 1.0;
      else
         rounded_sol[col] -= 1.0;
   }
}
