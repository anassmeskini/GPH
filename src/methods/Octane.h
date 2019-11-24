#ifndef OCTANE_HPP
#define OCTANE_HPP

#include "core/Heuristic.h"
#include "core/Numerics.h"
#include "dynamic_bitset/dynamic_bitset.hpp"

#include <tuple>
#include <vector>

class Octane : public FeasibilityHeuristic
{
 public:
   Octane() : FeasibilityHeuristic("Octane") {}

   void search(const MIP&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<Activity>&,
               const LPResult&, const std::vector<double>&,
               const std::vector<int>&, std::shared_ptr<const LPSolver>,
               TimeLimit, SolutionPool&) override;

 private:
   enum class Ray
   {
      OBJECTIVE,
      AVERAGE,
   };

#ifndef NDEBUG
   using Bitset = std::vector<bool>;
#else
   using Bitset = dynamic_bitset<>;
#endif

   struct FacetInfo
   {
      FacetInfo() = default;

      FacetInfo(const Bitset& f, double n, double d, int m)
          : facet(f), lambda_num(n), lambda_denom(d), min_plus(m)
      {
      }

      Bitset facet;
      double lambda_num;
      double lambda_denom;
      int min_plus;
   };

   static FacetInfo
   getFirstFacet(const std::vector<double>& origin,
                 const std::vector<double>& ray,
                 const std::vector<double>& column_ratios);

   static std::vector<FacetInfo>
   getKFacets(FacetInfo&& firstfacet, const std::vector<double>& origin,
              const std::vector<double>& ray,
              const std::vector<double>& column_ratios);

   static std::vector<FacetInfo>
   getReverseF(const FacetInfo& facet, const std::vector<double>& origin,
               const std::vector<double>& ray,
               const std::vector<double>& column_ratios);

   static std::vector<double>
   getSolFromFacet(const Bitset& facet,
                   const std::vector<int>& flipped_signs,
                   const std::vector<int>& permutation);

   static std::vector<double> makeRay(Ray raytype,
                                      const std::vector<int>& columns,
                                      const MIP& mip, const LPResult& res);

   constexpr static int fmax = 100;
   constexpr static int candidate_max = 5000;
   constexpr static bool use_frac_subspace = true;
   constexpr static bool scale_ray = true;
};

#endif
