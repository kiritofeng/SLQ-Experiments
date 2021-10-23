#include <algorithm>
#include <utility>

#include "pagerank.h"
#include "util.h"

std::tuple<std::unordered_set<size_t>, double, std::tuple<int, int, int, int>>
weighted_local_sweep_cut(const Graph<double>& G, const std::unordered_map<int, double> &x,
    int Gvol) {

  int n = G.size();
  std::vector<std::pair<int, double>> sx(x.begin(), x.end());
  std::sort(sx.begin(), sx.end(),
      [](const std::pair<int, double> &p1, const std::pair<int, double> &p2){
        return p1.second > p2.second;
    });

  std::unordered_set<size_t> S;
  double volS = 0.0;
  double cutS = 0.0;
  double bestcond = 1.0;
  std::tuple<int, int, int, int> beststats = std::make_tuple(1, 1, 1, Gvol - 1);
  std::unordered_set<size_t> bestset;

  for (const auto &p : sx) {
    if (S.size() == n - 1)
      break;

    size_t u = p.first;
    volS += G.deg(u);

    for (const auto &a : G.neighbours(u)) {
      int v;
      double ew;
      std::tie(v, ew) = a;
      if (S.count(v)) {
        cutS -= ew;
      } else {
        cutS += ew;
      }
    }
    S.insert(u);
    // TODO: will I die to FPE?
    if (cutS / std::min(volS, Gvol - volS) <= bestcond) {
      bestcond = cutS / std::min(volS, Gvol - volS);
      bestset = S;
      beststats = std::make_tuple(cutS, std::min(volS, Gvol - volS), volS, Gvol - volS);
    }
  }
  return std::make_tuple(std::move(bestset), bestcond, std::move(beststats));
}

std::tuple<std::unordered_set<size_t>, double>
round_to_cluster(const Graph<double> &G, std::vector<double> x) {
  std::unordered_map<int, double> nnz_dict;
  for (const auto &a : enumerate(x)) {
    int i; double v;
    std::tie(i, v) = a;
    if (v > 0) {
      nnz_dict[i] = v;
    }
  }
  std::unordered_set<size_t> bestset;
  double bestcond;
  std::tuple<int, int, int, int> beststats;
  int Gvol = 0;
  for (int d : G.deg()) {
    Gvol += d;
  }
  std::tie(bestset, bestcond, beststats) = weighted_local_sweep_cut(G, nnz_dict, Gvol);
  return std::make_tuple(std::move(bestset), bestcond);
}

