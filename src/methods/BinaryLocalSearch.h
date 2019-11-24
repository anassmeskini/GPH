#ifndef BIN_LS_HPP
#define BIN_LS_HPP

#include "core/Heuristic.h"

#include <queue>
#include <random>
#include <tuple>

class BinLocalSearch final : public ImprovementHeuristic
{
 public:
   BinLocalSearch() : ImprovementHeuristic("BinLS") {}

   void improve(const MIP&,                   // original problem
                const std::vector<double>&,   // lb at the node
                const std::vector<double>&,   // ub at the node
                const std::vector<Activity>&, // activities
                const LPResult&, // LP solution at the current node
                const std::vector<double>&, // activities of the rows at
                                            // the LP solution
                const std::vector<int>&,    // integer variables with
                                            // fractional values
                const std::vector<double>&, // best integer solution
                double,                     // cost of the best solution
                std::shared_ptr<const LPSolver>, // lp solver
                TimeLimit limit,                 // time limit
                SolutionPool&) override;

   void
   setParam(const std::string& param,
            const std::variant<std::string, int, double>& value) override
   {
      if (param == "Q0")
      {
         // TODO make sure it is nonegative
         Q[0] = std::get<int>(value);
      }
      else if (param == "Q1")
      {
         Q[1] = std::get<int>(value);
      }
      else if (param == "Q2")
      {
         Q[2] = std::get<int>(value);
      }
      else if (param == "stallitermax")
      {
         stallitermax = std::get<int>(value);
      }
   }

 private:
   struct Solution
   {
      Solution(const MIP& mip, const dynamic_bitset<>& sol,
               std::shared_ptr<uint32_t[]> coefs)
          : cost(0.0), values(mip.getNCols()),
            violation(2 * mip.getNRows(), false),
            looseness(2 * mip.getNRows(), false),
            activities(mip.getNRows(), 0), nviolated(0), nloose(0),
            large_violation(false), violated_row(-1), loose_row(-1),
            violation_hash(0), looseness_hash(0), hash_coefs(coefs)
      {
         auto st = mip.getStats();
         const auto& rhs = mip.getRHS();
         const auto& lhs = mip.getLHS();
         const auto& objective = mip.getObj();

         for (int col = 0; col < st.ncols; ++col)
         {
            assert(sol[col] == 1.0 || sol[col] == 0.0);
            values[col] = sol[col];
            cost += objective[col] * sol[col];
         }

         for (int row = 0; row < st.nrows; ++row)
         {
            auto [coefs, indices, size] = mip.getRow(row);

            for (int i = 0; i < size; ++i)
            {
               activities[row] += coefs[i] * sol[indices[i]];
            }

            if (!Num::isFeasGE(activities[row], lhs[row]))
            {
               violation[2 * row] = true;
               violation_hash += hash_coefs[row];
               ++nviolated;
            }
            else if (!Num::isFeasLE(activities[row], rhs[row]))
            {
               violation[2 * row + 1] = true;
               violation_hash += hash_coefs[row];
               ++nviolated;
            }

            if (!Num::isMinusInf(lhs[row]) &&
                Num::isFeasGE(activities[row], lhs[row] + 1))
            {
               looseness[2 * row] = true;
               looseness_hash += hash_coefs[row];
               loose_row = row;
               ++nloose;
            }
            else if (!Num::isInf(rhs[row]) &&
                     Num::isFeasLE(activities[row], rhs[row] - 1))
            {
               looseness[2 * row + 1] = true;
               looseness_hash += hash_coefs[row];
               loose_row = row;
               ++nloose;
            }
         }
      }

      Solution() = default;
      // Solution& operator=(Solution&&) = default;
      Solution& operator=(const Solution&) = default;
      Solution(Solution&&) = default;
      Solution(const Solution&) = default;

      std::pair<int, int> flip(const MIP& mip, int col)
      {
         const auto& lhs = mip.getLHS();
         const auto& rhs = mip.getRHS();
         const auto& objective = mip.getObj();
         auto [coefs, indices, size] = mip.getCol(col);

         double oldval = values[col];
         values[col] = !values[col];
         double newval = values[col];
         cost += objective[col] * (newval - oldval);

         int viol_diff = 0;
         int loose_diff = 0;

         for (int i = 0; i < size; ++i)
         {
            int row = indices[i];
            double coef = coefs[i];
            activities[row] += coef * (newval - oldval);

            bool lhs_violated = !Num::isFeasGE(activities[row], lhs[row]);
            bool rhs_violated = !Num::isFeasLE(activities[row], rhs[row]);
            bool lhs_loose =
                Num::isFeasGE(activities[row], lhs[row] + 1.0);
            bool rhs_loose =
                Num::isFeasLE(activities[row], rhs[row] - 1.0);

            large_violation |=
                !Num::isFeasGE(activities[row], lhs[row] - 2.0) ||
                !Num::isFeasLE(activities[row], rhs[row] + 2.0);

            // check if the lhs changed from feasible -> violated ot
            // violated -> feasible
            if ((lhs_violated && !violation[2 * row]) ||
                (rhs_violated && !violation[2 * row + 1]))
            {
               ++viol_diff;
               ++nviolated;
               violation_hash += hash_coefs[row];
               violated_row = row;
            }

            if ((!lhs_violated && violation[2 * row]) ||
                (!rhs_violated && violation[2 * row + 1]))
            {
               ++viol_diff;
               --nviolated;
               violation_hash -= hash_coefs[row];
               violated_row = -1;
               assert(nviolated >= 0);
            }

            if ((lhs_loose && !looseness[2 * row] &&
                 !looseness[2 * row + 1]) ||
                (rhs_loose && !looseness[2 * row + 1] &&
                 !looseness[2 * row]))
            {
               ++loose_diff;
               ++nloose;
               looseness_hash += hash_coefs[row];
               loose_row = row;
            }

            if ((!lhs_loose && looseness[2 * row] &&
                 !looseness[2 * row + 1]) ||
                (!rhs_loose && looseness[2 * row + 1] &&
                 !looseness[2 * row]))
            {
               ++loose_diff;
               --nloose;
               looseness_hash -= hash_coefs[row];
               loose_row = -1;
               assert(nviolated >= 0);
            }

            violation[2 * row] = lhs_violated;
            violation[2 * row + 1] = rhs_violated;
            looseness[2 * row] = lhs_loose;
            looseness[2 * row + 1] = rhs_loose;
         }

         return {viol_diff, loose_diff};
      }

      int get_violated_row() const
      {
         for (size_t row = 0; row < violation.size(); ++row)
         {
            if (violation[2 * row] || violation[2 * row + 1])
               return row;
         }

         return -1;
      }

      int get_loose_row() const
      {
         for (size_t row = 0; row < violation.size(); ++row)
         {
            if (looseness[2 * row] || looseness[2 * row + 1])
               return row;
         }

         return -1;
      }

      double cost;

      dynamic_bitset<> values;
      dynamic_bitset<> violation;
      dynamic_bitset<> looseness;
      std::vector<double> activities;

      int nviolated;
      int nloose;
      bool large_violation;

      int violated_row;
      int loose_row;

      uint32_t violation_hash;
      uint32_t looseness_hash;
      std::shared_ptr<uint32_t[]> hash_coefs;
   };

   class SolutionList
   {
    public:
      // TODO better sorting criteria
      SolutionList(int nrows_p, int subregion_size_p)
          : nrows(nrows_p), subregion_size(subregion_size_p),
            boxes(2 * (nrows + (nrows - 1) * subregion_size), Num::infval),
            region1_offset(0), region2_offset(nrows),
            region3_offset(nrows + (nrows - 1) * subregion_size),
            region4_offset(2 * nrows + (nrows - 1) * subregion_size),
            cmp([](const Solution& l, const Solution& r) -> bool {
               return l.nviolated + 2 * l.cost > r.nviolated + 2 * r.cost;
            }),
            solutions(cmp)
      {
      }

      bool add(Solution&& sol)
      {
         if (sol.nviolated == 1)
         {
            int violated_row = sol.get_violated_row();
            if (sol.cost < boxes[region1_offset + violated_row])
            {
               boxes[region1_offset + violated_row] = sol.cost;
               solutions.emplace(std::move(sol));
            }
         }
         else if (sol.nviolated > 1)
         {
            int boxid = sol.violation_hash % subregion_size;

            if (sol.cost <
                boxes[region2_offset + sol.nviolated - 2 + boxid])
            {
               boxes[region2_offset + sol.nviolated - 2 + boxid] =
                   sol.cost;
               solutions.emplace(std::move(sol));

               return true;
            }
         }
         else if (sol.nloose == 1)
         {
            int loose_row = sol.get_loose_row();
            if (sol.cost < boxes[region3_offset + loose_row])
            {
               boxes[region3_offset + loose_row] = sol.cost;
               solutions.emplace(std::move(sol));
            }
         }
         else if (sol.nloose > 1)
         {
            int boxid = sol.looseness_hash % subregion_size;

            if (sol.cost < boxes[region4_offset + sol.nloose - 2 + boxid])
            {
               boxes[region4_offset + sol.nloose - 2 + boxid] = sol.cost;
               solutions.emplace(std::move(sol));

               return true;
            }
         }
         else
         {
            // TODO
            assert(0);
         }

         return false;
      }

      void clear()
      {
         std::fill(std::begin(boxes), std::end(boxes), Num::infval);
         solutions = Queue(cmp);
      }

      auto pop()
      {
         auto sol = solutions.top();
         solutions.pop();
         return sol;
      }

      size_t size() const { return solutions.size(); }

      bool empty() const { return solutions.empty(); }

    private:
      int nrows;
      int subregion_size;

      std::vector<double> boxes;
      int region1_offset;
      int region2_offset;
      int region3_offset;
      int region4_offset;

      using Queue = std::priority_queue<
          Solution, std::vector<Solution>,
          std::function<bool(const Solution&, const Solution&)>>;
      std::function<bool(const Solution&, const Solution&)> cmp;
      Queue solutions;
   };

   std::array<int, 3> Q = {6, 8, 10};
   std::array<int, 3> region_size = {3, 8, 30};

   int stallitermax = 1000;
};

#endif
