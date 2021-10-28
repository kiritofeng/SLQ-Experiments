#ifndef GRAPH_H
#define GRAPH_H
#include <cstddef>
#include <vector>
#include <utility>

template<typename W>
class Graph {
  size_t n;
  std::vector<std::vector<std::pair<size_t, W>>> adj;

public:
  Graph(size_t n=0) : n{n}, adj(n) {}

  //Graph(const Graph<W> &G): n{G.n}, adj{G.adj} {}

  //Graph(Graph<W> &&G): n{G.n}, adj{std::move(G.adj)} {}

  //Graph &operator=(const Graph<W> &G) {
  //  n = G.n;
  //  adj = G.adj;
  //  return *this;
  //}
  //Graph &operator=(Graph<W> &&G) {
  //  n = G.n;
  //  adj = std::move(G.adj);
  //  return *this;
  //}

  size_t size() const { return n; }

  void add_edge(size_t u, size_t v, const W &w) {
    adj[u-1].emplace_back(v, w);
  }

  const std::vector<size_t> deg() const {
    std::vector<size_t> ret(n, 0);
    for (size_t i = 0; i < n; ++i) {
      ret[i] = adj[i].size();
    }
    return ret;
  }

  size_t deg(size_t u) const { return adj[u - 1].size(); }

  const std::vector<std::pair<size_t, W>> &neighbours(size_t u) const {
    return adj[u - 1];
  }
};
#endif // GRAPH_H
