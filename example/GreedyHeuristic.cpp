#include "GreedyHeuristic.h"
#include "io/Message.h"

#include <numeric>
#include <random>

void
GreedyHeuristic::search(const MIP& mip, const std::vector<double>&,
                        const std::vector<double>&,
                        const std::vector<Activity>&, const LPResult&,
                        const std::vector<double>&,
                        const std::vector<int>&,
                        std::shared_ptr<const LPSolver>, TimeLimit,
                        SolutionPool& pool)
{
   auto st = mip.getStats();
   int ncols = st.ncols;
   int nrows = st.nrows;
   const auto& lhs = mip.getLHS();
   const auto& rhs = mip.getRHS();
   const auto& objective = mip.getObj();

#ifndef NDEBUG
   for (int row = 0; row < nrows; ++row)
   {
      assert(lhs[row] == 1.0 && Num::isInf(rhs[row]));

      auto [coefs, _, size] = mip.getRow(row);

      for (int i = 0; i < size; ++i)
         assert(coefs[i] == 1.0);
   }

   for (int col = 0; col < ncols; ++col)
      assert(objective[col] > 0.0);
#endif

   dynamic_bitset<> best_cover;
   double best_cost = Num::infval;
   int iter = 0;

   static thread_local std::default_random_engine gen(0);
   std::uniform_real_distribution<double> rdist(0, 1.0);

   std::vector<int> nelements_left_original(ncols);
   for (int col = 0; col < ncols; ++col)
      nelements_left_original[col] = mip.getColSize(col);

   std::vector<int> nelements_left;
   dynamic_bitset<> in_cover;
   dynamic_bitset<> covered;
   std::vector<int> cover;
   std::vector<int> ncoverage;
   std::vector<int> candidates;

   cover.reserve(ncols / 10);
   candidates.reserve(10);
   do
   {
      nelements_left = nelements_left_original;
      in_cover = dynamic_bitset<>(ncols, false);
      covered = dynamic_bitset<>(nrows, false);
      cover.clear();

      int uncovered = nrows;
      double cost = 0.0;

      ncoverage = std::vector<int>(nrows, 0);
      do
      {
         double max_ratio = -1.0;
         int best_col = -1;
         for (int col = 0; col < ncols; ++col)
         {
            if (in_cover[col])
               continue;

            assert(objective[col] > 0.0);
            double col_ratio = nelements_left[col] / objective[col];

            if (col_ratio > max_ratio)
            {
               max_ratio = col_ratio;
               best_col = col;
            }
         }

         double p = rdist(gen);
         if (p > priority)
         {
            candidates.clear();
            for (int col = 0;
                 col < ncols && candidates.size() < cand_list_max; ++col)
            {
               double col_ratio = nelements_left[col] / objective[col];

               if (col_ratio * (1.0 + improvement) > max_ratio)
                  candidates.push_back(col);
            }

            std::uniform_int_distribution<int> idist(1, candidates.size());

            best_col = candidates[idist(gen) - 1];
         }

         assert(!in_cover[best_col]);
         assert(best_col >= 0);
         in_cover[best_col] = true;
         cover.push_back(best_col);
         cost += objective[best_col];

         auto [_, colindices, colsize] = mip.getCol(best_col);

         for (int i = 0; i < colsize; ++i)
         {
            int row = colindices[i];
            ++ncoverage[row];
            if (!covered[row])
            {
               --uncovered;
               covered[row] = true;
               auto [_, rowindices, rowsize] = mip.getRow(row);
               for (int j = 0; j < rowsize; ++j)
                  --nelements_left[rowindices[j]];
            }
         }

      } while (uncovered > 0);
      Message::debug("SCPGreedy: cover has {} columns ({:0.1f}%)",
                     in_cover.count(),
                     in_cover.count() / static_cast<double>(ncols) * 100);

      std::vector<int> mincoverage(ncols, ncols + 1);

      for (size_t i = 0; i < cover.size(); ++i)
      {
         int col = cover[i];
         auto [_, indices, size] = mip.getCol(col);

         for (int i = 0; i < size; ++i)
            mincoverage[col] =
                std::min(mincoverage[col], ncoverage[indices[i]]);
      }

      int nredundant = 0;
      for (int col = 0; col < ncols; ++col)
      {
         if (in_cover[col] && mincoverage[col] > 1)
         {
            in_cover[col] = false;
            ++nredundant;
            cost -= objective[col];

            assert(cost > 0.0);

            auto [_, colindices, colsize] = mip.getCol(col);

            for (int i = 0; i < colsize; ++i)
            {
               int row = colindices[i];
               --ncoverage[row];

               auto [_, rowindices, rowsize] = mip.getRow(row);
               for (int j = 0; j < rowsize; ++j)
               {
                  int other_col = rowindices[j];
                  mincoverage[other_col] =
                      std::min(mincoverage[other_col], ncoverage[row]);
               }
            }
         }
      }

      Message::debug(
          "SCPGreedy: removed {} redundant columns pop count {}",
          nredundant, in_cover.count());

      if (cost < best_cost)
      {
         best_cost = cost;
         best_cover = in_cover;
      }

   } while (++iter < itermax);

   std::vector<double> solution(ncols);
   for (int col = 0; col < ncols; ++col)
      solution[col] = static_cast<double>(best_cover[col]);
   pool.add(std::move(solution), best_cost);
}
