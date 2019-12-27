#include "BinaryLocalSearch.h"
#include "io/Message.h"

#include <array>
#include <random>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>

void
BinLocalSearch::improve(const MIP& mip, const std::vector<double>&,
                        const std::vector<double>&,
                        const std::vector<Activity>&,
                        const std::vector<double>& incumbent,
                        double incumbent_cost,
                        std::shared_ptr<const LPSolver>, TimeLimit tlimit,
                        SolutionPool& pool)
{
   auto st = mip.getStats();
   if (st.nbin < st.ncols || tlimit.reached(Timer::now()))
      return;

   std::shared_ptr<uint32_t[]> hash_coefs(new uint32_t[st.nrows]);

   {
      std::default_random_engine gen;
      std::uniform_int_distribution<uint32_t> dist(
          1, std::numeric_limits<uint32_t>::max());

      for (int row = 0; row < st.nrows; ++row)
         hash_coefs[row] = dist(gen);
   }

   double best_cost = incumbent_cost;
   dynamic_bitset<> best_sol(st.ncols);
   for (int col = 0; col < st.ncols; ++col)
   {
      assert(Num::round(incumbent[col]) == 1.0 ||
             Num::round(incumbent[col] == 0.0));
      best_sol[col] = static_cast<bool>(Num::round(incumbent[col]));
   }

   size_t run = 0;
   int iter;
   int stalliter;
   Solution adj_sol;
   do
   {
      iter = 0;
      stalliter = 0;

      SolutionList list(st.nrows, region_size[run]);
      list.add(Solution(mip, best_sol, hash_coefs));

      Message::debug("new run Q {} reg_size {}", Q[run], region_size[run]);

      do
      {
         auto cur_sol = list.pop();
         Message::debug("BinLs: iter {} size {} cost {}, nviolated {}",
                        iter, list.size(), cur_sol.cost,
                        cur_sol.nviolated);

         bool found_improvement = false;
         for (int col = 0; col < st.ncols; ++col)
         {
            adj_sol = cur_sol;

            auto [viol_diff, loos_diff] = adj_sol.flip(mip, col);
            if (adj_sol.large_violation)
               continue;

            if (adj_sol.nviolated == 0 && adj_sol.cost < best_cost)
            {
               // TODO improve this
               Message::debug("BinLs: new best solution {} -> {}",
                              best_cost, adj_sol.cost);
               best_sol = adj_sol.values;
               best_cost = adj_sol.cost;
               list.clear();
               list.add(std::move(adj_sol));
               stalliter = 0;
               found_improvement = true;
            }
            else if (!found_improvement &&
                     viol_diff + loos_diff < st.nrows * Q[run] / 100)
               list.add(std::move(adj_sol));
         }

         ++stalliter;
      } while (stalliter < stallitermax && !list.empty() &&
               !tlimit.reached(Timer::now()));
      ++iter;
   } while (!tlimit.reached(Timer::now()) && ++run < Q.size());

   if (best_cost < incumbent_cost)
   {
      std::vector<double> solution(st.ncols);
      for (int col = 0; col < st.ncols; ++col)
         solution[col] = static_cast<double>(best_sol[col]);
      pool.add(std::move(solution), best_cost);
   }
}
