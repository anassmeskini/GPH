#include "Octane.h"
#include "core/Common.h"
#include "core/Numerics.h"
#include "io/Message.h"

void
Octane::search(const MIP& mip, const std::vector<double>& lb,
               const std::vector<double>& ub, const std::vector<Activity>&,
               const LPResult& result, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               TimeLimit, SolutionPool& pool)

{
   int ncols = mip.getNCols();
   const auto& objective = mip.getObj();
   if (mip.getStats().nbin < mip.getNCols())
      return;

   // get the active columns in Octane
   std::vector<int> active_columns;
   std::vector<int> fixed_columns;
   if (use_frac_subspace)
   {
      Message::debug_details("Octane: using fractional space");
      for (int col = 0; col < ncols; ++col)
      {
         if (Num::isFeasInt(result.primalSolution[col]))
            fixed_columns.push_back(col);
         else
            active_columns.push_back(col);
      }
   }
   else
   {
      Message::debug_details("Octane: using the entire space");
      for (int col = 0; col < ncols; ++col)
      {
         if (!Num::isFeasEQ(lb[col], ub[col]))
            active_columns.push_back(col);
         else
            fixed_columns.push_back(col);
      }
   }

   Message::debug_details("Octane: active cols: {}, fixed cols: {}",
                          active_columns.size(), fixed_columns.size());

   assert(static_cast<int>(fixed_columns.size() + active_columns.size()) ==
          ncols);

   if (std::all_of(std::begin(active_columns), std::end(active_columns),
                   [&objective](int col) -> bool {
                      return objective[col] == 0.0;
                   }))
   {
      Message::debug("Octane: Active ray = zero");
      return;
   }

   // the origin is the shifted optimal lp vertex
   std::vector<double> origin(active_columns.size());
   for (size_t i = 0; i < active_columns.size(); ++i)
   {
      int col = active_columns[i];
      origin[i] = result.primalSolution[col] - 0.5;
   }

   // create the shooting ray
   Ray raytype = Ray::OBJECTIVE;
   std::vector<double> ray = makeRay(raytype, active_columns, mip, result);

   // scale the ray
   if (scale_ray)
   {
      double orig_norm =
          0.5 * std::sqrt(static_cast<double>(origin.size()));
      double ray_norm = 0.0;
      for (double val : ray)
         ray_norm += val * val;

      ray_norm = std::sqrt(ray_norm);

      for (double& val : ray)
         val = orig_norm / ray_norm * val;

      Message::debug_details("Octane: origin norm: {}, ray norm {}",
                             orig_norm, ray_norm);
   }

   // what columns had changed signs to make
   // the ray nonnegative
   std::vector<int> flipped_signs;
   // flip signs and compute column ratios
   std::vector<double> column_ratios(active_columns.size());
   for (size_t i = 0; i < active_columns.size(); ++i)
   {
      if (ray[i] < 0)
      {
         ray[i] *= -1.0;
         origin[i] *= -1.0;
         flipped_signs.push_back(i);
      }

      if (Num::isFeasEQ(ray[i], 0.0))
      {
         if (origin[i] > 0)
            column_ratios[i] = -Num::infval;
         else
            column_ratios[i] = Num::infval;
      }
      else
         column_ratios[i] = -origin[i] / ray[i];
   }

   Message::debug_details("Octane: flipped signs for {} cols",
                          flipped_signs.size());

   // columns are permuted to have a nonincreasing order of column_ratios
   // (v_i)
   std::vector<int> active_cols_perm = getIdentity(active_columns.size());
   std::sort(std::begin(active_cols_perm), std::end(active_cols_perm),
             [&column_ratios](int left, int right) {
                return column_ratios[left] > column_ratios[right];
             });

   // permute the column ratios
   {
      std::vector<double> buffer(column_ratios.size());
      for (size_t i = 0; i < column_ratios.size(); ++i)
      {
         buffer[i] = column_ratios[active_cols_perm[i]];
      }

      column_ratios = std::move(buffer);

      assert(std::is_sorted(std::rbegin(column_ratios),
                            std::rend(column_ratios)));
   }

   // permute the origin
   {
      std::vector<double> buffer(column_ratios.size());
      for (size_t i = 0; i < column_ratios.size(); ++i)
      {
         buffer[i] = origin[active_cols_perm[i]];
      }

      origin = std::move(buffer);
   }

   // permute the ray
   {
      std::vector<double> buffer(column_ratios.size());
      for (size_t i = 0; i < column_ratios.size(); ++i)
      {
         buffer[i] = ray[active_cols_perm[i]];
      }

      ray = std::move(buffer);
   }

   assert(std::all_of(std::begin(ray), std::end(ray),
                      [](double v) { return v >= 0; }));
   assert(std::all_of(std::begin(origin), std::end(origin),
                      [](double v) { return v >= -0.5 && v <= 0.5; }));

   auto firstFacet = getFirstFacet(origin, ray, column_ratios);

   if (Num::isFeasEQ(firstFacet.lambda_denom, 0.0))
   {
      Message::debug("Octane: lambda denom numerically = 0");
      return;
   }

   // get the k closest facets
   auto k_closest_facets =
       getKFacets(std::move(firstFacet), origin, ray, column_ratios);

#ifndef NDEBUG
   // check that we did not generate the same facets twice
   auto equal = [](const Bitset& l, const Bitset& r) -> bool {
      assert(l.size() == r.size());
      for (size_t k = 0; k < l.size(); ++k)
      {
         if (l[k] != r[k])
            return false;
      }
      return true;
   };

   for (size_t i = 0; i < k_closest_facets.size(); ++i)
   {
      for (size_t j = i + 1; j < k_closest_facets.size(); ++j)
         if (equal(k_closest_facets[i].facet, k_closest_facets[j].facet))
         {
            Message::debug_details("Octane: two similar facets {} and {}",
                                   i, j);
         }
   }
#endif

   Message::debug_details("Octane: generated {} facets",
                          k_closest_facets.size());

   // check feasibility of the solutions corresponding to the facets
   for (auto& [facet, ln, ld, id] : k_closest_facets)
   {

#ifndef NDEBUG
      std::vector<double> solution(ncols);
#else
      std::vector<double> solution(ncols, Num::infval);
#endif

      double cost = 0.0;
      // get the solution in the original search space
      std::vector<double> partial_solution =
          getSolFromFacet(facet, flipped_signs, active_cols_perm);

      assert(partial_solution.size() == active_columns.size());

      // set the active columns
      for (size_t i = 0; i < active_columns.size(); ++i)
      {
         int col = active_columns[i];
         solution[col] = partial_solution[i];
         cost += solution[col] * objective[col];
      }

      // set the fixed columns
      for (size_t i = 0; i < fixed_columns.size(); ++i)
      {
         int col = fixed_columns[i];
         solution[col] = result.primalSolution[col];
         cost += solution[col] * objective[col];
      }

      assert(std::all_of(
          std::begin(solution), std::end(solution), [](double val) {
             return val != Num::infval && Num::isIntegral(val);
          }));

      static auto lpFeas = checkFeasibility<double, true>;
      if (lpFeas(mip, solution, 1e-6, 1e-6))
      {
         Message::debug_details("Octane: found a feasible solution");
         pool.add(std::move(solution), cost);
      }
   }
}

Octane::FacetInfo
Octane::getFirstFacet(const std::vector<double>& origin,
                      const std::vector<double>& ray,
                      const std::vector<double>& column_ratios)
{
   int space_size = origin.size();

   assert(ray.size() == origin.size());
   assert(std::all_of(std::begin(ray), std::end(ray),
                      [](double val) { return val >= 0.0; }));
   assert(std::is_sorted(std::rbegin(column_ratios),
                         std::rend(column_ratios)));

   double lambda_num = 0.5 * space_size;
   double lambda_denom = 0.0;

   // we start from the [1, 1, ... 1, 1] facet
   Bitset facet;
   for (int i = 0; i < space_size; ++i)
   {
      facet.push_back(true);
      lambda_num += -origin[i];
      lambda_denom += ray[i];
   }

   assert(lambda_num >= 0);
   assert(lambda_denom > 0);

   for (int i = 0; i < space_size; ++i)
   {
      if (column_ratios[i] * lambda_denom > lambda_num)
      {
         assert(facet[i]);
         facet[i] = 0;
         lambda_num += 2.0 * origin[i];
         lambda_denom += -2.0 * ray[i];

         assert(lambda_num >= 0);
         assert(lambda_denom > 0);
      }
      else
         break;
   }

   Message::debug_details("Octane: First facet: lambda {}",
                          lambda_num / lambda_denom);

   int min_plus_index = -1;
   for (int i = 0; i < space_size; ++i)
   {
      if (facet[i])
      {
         min_plus_index = i;
         break;
      }
   }

   assert(min_plus_index >= 0);

   return {std::move(facet), lambda_num, lambda_denom, min_plus_index};
}

std::vector<Octane::FacetInfo>
Octane::getKFacets(FacetInfo&& firstFacet,
                   const std::vector<double>& origin,
                   const std::vector<double>& ray,
                   const std::vector<double>& column_ratios)
{
   assert(ray.size() == origin.size());
   assert(std::all_of(std::begin(ray), std::end(ray),
                      [](double val) { return val >= 0.0; }));

   auto insert = [](std::vector<FacetInfo>& unscanned, FacetInfo&& facet) {
      double facet_lambda = facet.lambda_num / facet.lambda_denom;

      int i = 0;
      while (i < static_cast<int>(unscanned.size()) &&
             unscanned[i].lambda_num <
                 facet_lambda * unscanned[i].lambda_denom)
      {
         ++i;
      }

      unscanned.resize(unscanned.size() + 1);
      int size = static_cast<int>(unscanned.size());

      for (int j = size - 2; j >= i; --j)
      {
         assert(j + 1 < static_cast<int>(unscanned.size()));
         assert(j >= 0);
         assert(j >= i);

         unscanned[j + 1] = std::move(unscanned[j]);
      }

      unscanned[i] = std::move(facet);
   };

   auto pop_front = [](std::vector<FacetInfo>& unscanned) {
      auto tmp = std::move(unscanned[0]);

      int size = unscanned.size();
      for (int i = 0; i < size - 1; ++i)
         unscanned[i] = std::move(unscanned[i + 1]);

      assert(!unscanned.empty());
      unscanned.resize(unscanned.size() - 1);

      return tmp;
   };

   // array of unscanned facets
   std::vector<FacetInfo> unscanned_facets = {firstFacet};
   std::vector<FacetInfo> scanned_facets;

   unscanned_facets.reserve(5 * fmax);
   scanned_facets.reserve(fmax);

   Message::debug_details("Octane: generating neighbor facets");
   int count = 0;
   do
   {
      assert(!unscanned_facets.empty());

      auto min_lambda_facet = pop_front(unscanned_facets);

      assert(!min_lambda_facet.facet.empty());

      auto neighbors =
          getReverseF(min_lambda_facet, origin, ray, column_ratios);

      Message::debug_details("Octane: current lambda {} node degree {}",
                             min_lambda_facet.lambda_num /
                                 min_lambda_facet.lambda_denom,
                             neighbors.size());

      for (auto& facet : neighbors)
         insert(unscanned_facets, std::move(facet));

      scanned_facets.emplace_back(std::move(min_lambda_facet));

      if (unscanned_facets.size() >= candidate_max)
      {
         Message::debug_details(
             "Octane: max number of unscanned candidates reached");

         size_t i = 0;
         while (scanned_facets.size() < fmax &&
                i < unscanned_facets.size())
         {
            scanned_facets.emplace_back(std::move(unscanned_facets[i++]));
         }

         break;
      }
   } while (++count < fmax && !unscanned_facets.empty());

   return scanned_facets;
}

std::vector<Octane::FacetInfo>
Octane::getReverseF(const FacetInfo& facetInfo,
                    const std::vector<double>& origin,
                    const std::vector<double>& ray,
                    const std::vector<double>& column_ratios)
{
   int space_size = origin.size();
   assert(ray.size() == origin.size());

   std::vector<FacetInfo> reverseF;

   auto& [facet, lambda_num, lambda_denom, min_plus_index] = facetInfo;
   double lambda = lambda_num / lambda_denom;

#ifndef NDEBUG
   int min_plus_dbg = -1;
   for (int i = 0; i < space_size; ++i)
   {
      if (facet[i])
      {
         min_plus_dbg = i;
         break;
      }
   }
   assert(min_plus_dbg >= 0);
   assert(min_plus_dbg == min_plus_index);
#endif

   assert(facet.size() == origin.size());

   int i = 0;
   while (i < space_size && static_cast<bool>(!facet[i]) &&
          column_ratios[i] > lambda)
   {
      double new_num = lambda_num - 2.0 * origin[i];
      double new_denom = lambda_denom + 2.0 * ray[i];
      assert(new_num / new_denom >= lambda);

      reverseF.emplace_back(facet, new_num, new_denom, i);
      reverseF.back().facet[i] = 1;

      ++i;
   }

   i = space_size - 1;

   double ratio_min_plus = -origin[min_plus_index] / ray[min_plus_index];
   while (i >= 0 && facet[i] && column_ratios[i] <= lambda)
   {
      double new_num = lambda_num + 2.0 * origin[i];
      double new_denom = lambda_denom - 2.0 * ray[i];
      double new_lambda = new_num / new_denom;

      if (new_denom > 0 && ratio_min_plus <= new_lambda)
      {
         assert(new_num / new_denom >= lambda);
         int new_min_plus = min_plus_index + (i == min_plus_index);

         reverseF.emplace_back(facet, new_num, new_denom, new_min_plus);
         reverseF.back().facet[i] = 0;

         assert(min_plus_index < space_size && min_plus_index >= 0);
      }

      --i;
   }

   return reverseF;
}

// the changes made to the original space:
// 1- restrict the search space to a subset of columns
// 2- flip signs of negative ray components
// 3- permute the columns for nonincreasing columns ratios
std::vector<double>
Octane::getSolFromFacet(const Bitset& facet,
                        const std::vector<int>& flipped_signs,
                        const std::vector<int>& permutation)
{
   int space_size = facet.size();
   std::vector<double> sol(space_size);
   std::vector<double> buffer(space_size);

   assert(static_cast<size_t>(space_size) == permutation.size());

   for (int i = 0; i < space_size; ++i)
      sol[i] = static_cast<double>(facet[i]);

   // undo the permutation
   for (int i = 0; i < space_size; ++i)
      buffer[permutation[i]] = sol[i];

   sol = std::move(buffer);

   // undo the flipped signs
   for (int col : flipped_signs)
      sol[col] = 1.0 - sol[col];

   return sol;
}

std::vector<double>
Octane::makeRay(Ray raytype, const std::vector<int>& columns,
                const MIP& mip, const LPResult&)
{
   int space_size = columns.size();
   std::vector<double> ray(space_size);

   if (raytype == Ray::OBJECTIVE)
   {
      const auto& objective = mip.getObj();

      for (int i = 0; i < space_size; ++i)
         ray[i] = -objective[columns[i]];
   }
   else
   {
      assert(0);
      return {};
   }

   return ray;
}
