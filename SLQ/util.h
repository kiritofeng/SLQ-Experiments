#ifndef UTIL_H
#define UTIL_H
#include <cstddef>
#include <iostream>
#include <iterator>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph.h"

template<typename InputIt> class enumerate_iterable;
template<typename InputIt> class enumerate_iterator;
template<typename T> enumerate_iterable<typename T::const_iterator> enumerate(const T &container);
template<typename InputIt> enumerate_iterable<InputIt> enumerate(const InputIt &begin, const InputIt &end);

template<typename InputIt> class enumerate_iterable {
private:
  const InputIt first, last;

  enumerate_iterable(InputIt first, InputIt last) : first{first}, last{last} {}

public:
  class iterator {
  private:
    size_t cnt;
    InputIt it;
    explicit iterator(InputIt it): cnt{0}, it{it} {}
    iterator(size_t cnt, InputIt it): cnt{cnt}, it{it} {}

    class postinc_return {
      size_t cnt;
      typename InputIt::value_type val;
      postinc_return(size_t cnt, typename InputIt::value_type val) : cnt{cnt}, val{val} {}

    public:
      std::pair<size_t, typename InputIt::value_type> operator *() {
        return std::make_pair(cnt, val);
      }
    };

  public:
      typedef std::pair<size_t, typename InputIt::value_type>                     value_type;
      typedef typename InputIt::difference_type          difference_type;
      typedef value_type*                    pointer;
      typedef value_type&                    reference;
      typedef std::input_iterator_tag iterator_category;

      value_type operator*() const { return std::make_pair(cnt, *it); }

      bool operator==(const iterator &other) const {
        return it == other.it;
      }

      bool operator!=(const iterator &other) const {
        return !(*this == other);
      }

      iterator &operator++() {
        cnt++;
        it++;
        return *this;
      }

      postinc_return operator++(int) {
        postinc_return ret(cnt, *it);
        ++*this;
        return ret;
      }
      friend class enumerate_iterable;
  };

  iterator begin() const { return iterator(first); }

  iterator end() const { return iterator(last); }

  friend enumerate_iterable<InputIt> enumerate<>(const InputIt &begin, const InputIt &end);
};

template<typename InputIt> enumerate_iterable<InputIt> enumerate(const InputIt &begin, const InputIt &end) {
  return enumerate_iterable<InputIt>(begin, end);
}

template<typename T> enumerate_iterable<typename T::const_iterator>
enumerate(const T &container) {
  return enumerate(container.begin(), container.end());
}

template<typename T>
std::tuple<double, double> compute_pr_rc(const std::unordered_set<T> &prediction,
    const std::unordered_set<T> &truth) {

    size_t common = 0;
  if (prediction.size() > truth.size()) {
    for (const T &t : truth) {
      if (prediction.count(t) > 0) {
        common += 1;
      }
    }
  } else {
    for (const T &t : prediction) {
      if (truth.count(t) > 0) {
        common += 1;
      }
    }
  }

  return std::make_tuple(1 - 1.0 * (prediction.size() - common) / prediction.size(),
      1 - 1.0 * (truth.size() - common) / truth.size());
}

template<typename T>
std::tuple<double, double> compute_pr_rc(const std::unordered_set<T> &prediction,
    const std::vector<T> &truth) {

  size_t common = 0;
  for (const T &t : truth) {
    if (prediction.count(t) > 0) {
      common += 1;
    }
  }

  return std::make_tuple(1 - 1.0 * (prediction.size() - common) / prediction.size(),
      1 - 1.0 * (truth.size() - common) / truth.size());
}

template<typename W>
size_t _component_dfs(size_t u, const Graph<W> &G, std::vector<bool> &vis) {
  size_t ret = 1;
  vis[u - 1] = 1;
  for (const auto &a : G.neighbours(u)) {
    if (!vis[a.first - 1]) {
      ret += _component_dfs(a.first, G, vis);
    }
  }
  return ret;
}

template<typename W>
void _component_dfs(size_t u, const Graph<W> &G, Graph<W> &cc, std::vector<bool> &vis) {
  vis[u - 1] = 1;
  for (const auto &a : G.neighbours(u)) {
    if (!vis[a.first - 1]) {
      cc.add_edge(u, a.first, a.second);
      _component_dfs(a.first, G, cc, vis);
    }
  }
}

template<typename W>
std::tuple<Graph<W>, std::vector<bool>> largest_component(const Graph<W> &G) {
  Graph<W> cc{G.size()};
  std::vector<bool> vis(G.size());
  size_t max_ind = 1, max_sz = 0;
  for (size_t u = 1; u <= G.size(); ++u) {
    if (!vis[u - 1]) {
      size_t sz = _component_dfs(u, G, vis);
      if (sz > max_sz) {
        max_ind = u;
        max_sz = sz;
      }
    }
  }
  for (size_t u = 0; u < G.size(); ++u) {
    vis[u] = 0;
  }
  _component_dfs(max_ind, G, cc, vis);
  return std::make_tuple(std::move(cc), std::move(vis));
}
#endif // UTIL_H
