#include "SLQ.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_set>

#include "util.h"

QHuberLoss::QHuberLoss(double q, double delta) : EdgeLoss{}, q{q}, delta{delta} {}

double QHuberLoss::operator()(double d) const { return pow(d, 1 / (q - 1)); }

double sign(double d) {
  if (d == 0) {
    return 0;
  } else {
    return copysign(1.0, d);
  }
}

double QHuberLoss::loss_gradient(double x) const {
  if (std::abs(x) < delta) {
    return pow(delta, q - 2) * x;
  } else {
    return sign(x) * pow(std::abs(x), q - 1);
  }
}

double QHuberLoss::loss_function(double x) const {
  if (std::abs(x) < delta) {
    return 0.5 * pow(std::abs(x), q - 2) * pow(x, 2);
  } else {
    return pow(std::abs(x), q) / q + (0.5 - 1 / q) * pow(delta, q);
  }
}

double TwoNormLoss::operator()(double d) const { return sqrt(d); }

double TwoNormLoss::loss_gradient(double x) const { return x; }

double TwoNormLoss::loss_function(double x) const { return 0.5 * pow(x, 2); }

std::shared_ptr<EdgeLoss> lossType(double q, double delta) {
  if (q == 2.0) {
    return std::make_shared<TwoNormLoss>();
  } else {
    return std::make_shared<QHuberLoss>(q, delta);
  }
}
size_t _buffer_neighbours(const std::vector<double> &x, const Graph<double> &G,
    int i, std::vector<double> &buf_x, std::vector<double> &buf_vals) {

  size_t nneighs = G.neighbours(i).size();
  for (const auto &a : enumerate(G.neighbours(i))) {
    int iter = a.first;
    int j = a.second.first;
    buf_x[iter] = x[j - 1];
    buf_vals[iter] = a.second.second;
  }
  return nneighs;
}

double _eval_residual_i(double xi, double di, double dx, bool seed,
    int nneighs, const std::vector<double> &neigh_x,
    const std::vector<double> &neigh_vals,
    std::shared_ptr<EdgeLoss> L, double gamma) {

  double ri_new = 0;
  for (int k = 0; k < nneighs; ++k) {
    ri_new -= neigh_vals[k] * L->loss_gradient(xi + dx - neigh_x[k]) / gamma;
  }
  if (seed) {
    ri_new -= di * L->loss_gradient(xi + dx - 1);
  } else {
    ri_new -= di * L->loss_gradient(xi + dx);
  }
    return ri_new;
}

template<typename W>
size_t _max_nz_degree(const Graph<W> &G) {
  size_t ret = 0;
  for (size_t u = 1; u <= G.size(); ++u) {
    ret = std::max(ret, G.deg(u));
  }
  return ret;
}

double dxi_solver(const Graph<double> &G, const std::vector<double> &x,
    double kappa, double epsilon, double gamma, const std::vector<double> &r,
    const std::unordered_set<int> &seedset, double rho, int i,
    std::shared_ptr<EdgeLoss> L, std::vector<double> &buf_x,
    std::vector<double> &buf_vals, double thd1, double thd2) {

  std::cout << "dxi_solver " << i << std::endl;
  int di = G.deg(i);
  bool found_dxi = false;
  int nneighs = _buffer_neighbours(x, G, i, buf_x, buf_vals);

  int nbisect = 0;

  double ri_new = r[i - 1];
  double dx_min = 0;
  double thd_min = std::min(thd1, thd2);
  double thd_max = std::max(thd1, thd2);
  double thd = thd_max;
  double dx = thd;

  ri_new = _eval_residual_i(x[i - 1], di, dx, seedset.count(i) > 0, nneighs,
      buf_x, buf_vals, L, gamma);
  if (ri_new < 0) {
    ri_new = r[i - 1];
    thd = thd_min;
  }
  double last_dx = 0;

  int ratio = 10;
  while (ri_new > rho * kappa * di) {
    std::cout << "ri_new = " << ri_new << std::endl;
    dx = thd;
    ri_new = _eval_residual_i(x[i - 1], di, dx, seedset.count(i) > 0,
        nneighs, buf_x, buf_vals, L, gamma);

    if (nbisect >= 40) {
      std::cout << i << " " << dx << " " << di << " " << ri_new << " " << rho * kappa * di << std::endl;
    }

    last_dx = dx_min;
    dx_min = thd;
    thd *= ratio;
    nbisect += 1;
  }
  dx_min = last_dx;
  double dx_max = thd/ratio;

  double dx_mid = 0;
  while ((found_dxi == false && dx_max - dx_min > epsilon) || ri_new < 0) {
    std::cout << "searching [" << dx_min << ", " << dx_max << "], ri_new = " << ri_new << std::endl;
    double dx_mid = dx_max / 2 + dx_min / 2;
    ri_new = _eval_residual_i(x[i - 1], di, dx_mid, seedset.count(i) > 0,
        nneighs, buf_x, buf_vals, L, gamma);

    if (ri_new < rho * kappa * di) {
      dx_max = dx_mid;
    } else if (ri_new > rho * kappa * di) {
      dx_min = dx_mid;
    } else {
      found_dxi = true;
    }
  }
  return (dx_mid == 0) ? dx_max : dx_mid;
}

void residual_update(const Graph<double> &G, const std::vector<double> &x,
    double dxi, int i, const std::unordered_set<int> &seedset,
    std::vector<double> &r, double gamma,
    CircularDeque<int> &Q, double kappa, std::shared_ptr<EdgeLoss> L) {

  r[i - 1] = 0;
  for (const auto &a : G.neighbours(i)) {
    int j = a.first;
    double dri = L->loss_gradient(x[j - 1] - x[i - 1] - dxi);
    double drij = a.second * (L->loss_gradient(x[j - 1] - x[i - 1]) - dri);
    drij /= gamma;
    double rj_old = r[j - 1];
    r[j - 1] += drij;
    r[i - 1] += a.second * dri / gamma;
    if (rj_old <= kappa * G.deg(j) && r[j - 1] > kappa * G.deg(j)) {
      Q.push(j);
    }
  }

  if (seedset.count(i)) {
    r[i - 1] -= G.deg(i) * L->loss_gradient(x[i - 1] + dxi - 1);
  } else {
    r[i - 1] -= G.deg(i) * L->loss_gradient(x[i - 1] + dxi);
  }
  if (r[i - 1] > kappa * G.deg(i)) {
    Q.push(i);
  }
}

std::tuple<std::vector<double>, std::vector<double>, int>
slq_diffusion(const Graph<double> &G, const std::vector<size_t> &S, double gamma,
    double kappa, double rho, std::shared_ptr<EdgeLoss> L, int max_iters,
    double epsilon) {

  int n = G.size();
  std::vector<double> x(n, 0), r(n, 0);

  int max_deg = _max_nz_degree(G);

  std::vector<double> buf_x(max_deg, 0), buf_vals(max_deg, 0);
  CircularDeque<int> Q(n);

  for (size_t i : S) {
    r[i - 1] = G.deg(i);
    Q.push(i);
  }
  std::unordered_set<int> seedset(S.begin(), S.end());

  int iter = 0;

  int Svol = 0, Gvol = 0;
  for (int u : S) {
    Svol += G.deg(u);
  }
  for (int u = 1; u <= n; ++u) {
    Gvol += G.deg(u);
  }

  double thd1 = (*L)(1.0 * Svol / Gvol);
  double thd2 = thd1;

  while (Q.length() > 0 && iter < max_iters) {
    int i = Q.popfirst();
    double dxi = dxi_solver(G, x, kappa, epsilon, gamma, r, seedset, rho, i, L,
        buf_x, buf_vals, thd1, thd2);
    thd2 = dxi;
    residual_update(G, x, dxi, i, seedset, r, gamma, Q, kappa, L);
    x[i - 1] += dxi;

    iter += 1;
    std::cout << "finished iteration " << iter <<std::endl;
  }

  if (iter == max_iters && Q.length() > 0) {
    std::cerr << "reached maximum iterations" << std::endl;
  }

  return std::make_tuple(std::move(x), std::move(r), iter);
}
