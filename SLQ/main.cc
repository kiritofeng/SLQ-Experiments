#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>
#include <vector>

#include "graph.h"
#include "pagerank.h"
#include "SLQ.h"
#include "util.h"

const int NTRIALS = 10;
const int LABEL = 2008;
const double DELTA = 0.0;

int main(int argc, char **argv) {
  // argv[1] is the graph
  // argv[2] is the labels
  std::ifstream graph_in(argv[1]), labels_in(argv[2]);
  int n, m;
  graph_in >> n >> m;
  Graph<double> G{n};
  for (int i = 0; i < m; ++i) {
    int u, v;
    graph_in >> u >> v;
    G.add_edge(u, v, 1.0);
  }
  Graph<double> GCC;
  std::vector<bool> p;
  std::tie(GCC, p) = largest_component(G);

  std::vector<std::tuple<double, double, double, double, double, double, double>> labels{n};
  for (int i = 0; i < n; ++i) {
    labels_in >> std::get<0>(labels[i]) >> std::get<1>(labels[i])
              >> std::get<2>(labels[i]) >> std::get<3>(labels[i])
              >> std::get<4>(labels[i]) >> std::get<5>(labels[i])
              >> std::get<6>(labels[i]);
  }
  std::vector<size_t> truth;
  for (int i = 0; i < n; ++i) {
    if (p[i] && std::get<5>(labels[i]) == LABEL) {
      truth.push_back(i + 1);
    }
  }
  std::mt19937 gen;
  for (int seed = 0; seed < NTRIALS; ++seed) {
    gen.seed(seed);
    std::shuffle(truth.begin(), truth.end(), gen);
    std::vector<size_t> S(truth.begin(), truth.begin() + std::max(1, (int)round(0.01 * truth.size())));
    std::shared_ptr<QHuberLoss> L = std::make_shared<QHuberLoss>(1.2, DELTA);

    std::vector<double> x_slq, r;
    int iter;
    std::tie(x_slq, r, iter) = slq_diffusion(GCC, S, 0.05, 0.005, 0.5, L, 100000, 1e-8);

    std::unordered_set<size_t> cluster_slq;
    double cond_slq;
    std::tie(cluster_slq, cond_slq) = round_to_cluster(GCC, x_slq);

    double pr_slq, rc_slq;
    std::tie(pr_slq, rc_slq) = compute_pr_rc(cluster_slq, truth);
    if (!(iter < 100000)) {
      cond_slq = 1.0;
    }
    std::cout << seed << " " << cond_slq << " " << pr_slq << " " << rc_slq;
    std::cout << " " << 2*pr_slq*rc_slq/(pr_slq+rc_slq) << std::endl;
  }
}
