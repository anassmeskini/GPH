#ifndef FEAS_PUMP_HPP
#define FEAS_PUMP_HPP
#include "core/Heuristic.h"
#include "core/Numerics.h"

#include <vector>

class FeasPump : public HeuristicMethod
{
 public:
   FeasPump() : HeuristicMethod("FeasPump") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               SolutionPool&) override;

 private:
   constexpr static int pert_freq = 100;
   constexpr static int max_stall_iter = 10;
   constexpr static int min_flips = 10;
   constexpr static int max_flips = 30;
   constexpr static int cycle_length = 3;
   constexpr static int max_iter = 1000;
   constexpr static int max_restarts = 50;
   constexpr static int min_iter_pert = 101;

   constexpr static double alpha_update_factor = 0.9;
   constexpr static double max_lp_iter_ratio = 0.01;
   constexpr static double min_frac_improv = 0.001;
   constexpr static double zero_frac = 1e-4;
   constexpr static double alpha_cycle_detection_threshold = 0.5;

   template <int cyclelength>
   int is_already_visited(
       const std::array<std::vector<double>, cyclelength>& history,
       const std::vector<double>& int_sol)
   {
      assert(cyclelength);
      int ncols = int_sol.size();
      std::array<bool, cyclelength> is_same;
      is_same.fill(true);

      for (int col = 0; col < ncols; ++col)
      {
         for (int j = 0; j < cyclelength; ++j)
            is_same[j] =
                is_same[j] && Num::isFeasEQ(int_sol[col], history[j][col]);
      }

      for (int j = 0; j < cyclelength; ++j)
      {
         if (is_same[j])
            return j;
      }

      return 0;
   }

   void handle_one_cycle(std::vector<double>&, const std::vector<double>&,
                         const std::vector<double>&,
                         const std::vector<double>&, int);

   void make_rand_perturbation(std::vector<double>&,
                               const std::vector<double>&, int,
                               const std::vector<double>&);

   double get_frac(const std::vector<double>&, int);
};

#endif
