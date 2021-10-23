#ifndef PAGERANK_H
#define PAGERANK_H
#include <queue>
#include <tuple>
#include <unordered_set>
#include <unordered_map>

#include "graph.h"

std::tuple<std::unordered_set<size_t>, double, std::tuple<int, int, int, int>>
weighted_local_sweep_cut(const Graph<double>& G,
    const std::unordered_map<int, double> &x, int Gvol);

std::tuple<std::unordered_set<size_t>, double>
round_to_cluster(const Graph<double> &G, std::vector<double> x);
#endif
