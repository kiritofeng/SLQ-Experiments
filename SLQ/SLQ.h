#ifndef SLQ_H
#define SLQ_H
#include <memory>
#include <tuple>
#include <unordered_set>

#include "circulardeque.h"
#include "graph.h"

class EdgeLoss {
public:
  virtual double operator()(double d) const = 0;
  virtual double loss_gradient(double x) const = 0;
  virtual double loss_function(double x) const = 0;
};

class QHuberLoss : public EdgeLoss {
  double q, delta;
public:
  QHuberLoss(double q, double delta);

  virtual double operator()(double d) const;
  virtual double loss_gradient(double x) const;
  virtual double loss_function(double x) const;
};

class TwoNormLoss : public EdgeLoss {
public:
  TwoNormLoss() = default;

  virtual double operator()(double d) const;
  virtual double loss_gradient(double x) const;
  virtual double loss_function(double x) const;
};

std::shared_ptr<EdgeLoss> lossType(double q, double delta);

double dxi_solver(const Graph<double> &G, const std::vector<double> &x, double kappa,
    double epsilon, double gamma, const std::vector<double> &r,
    const std::unordered_set<int> &seedset, double rho, int i,
    std::shared_ptr<EdgeLoss> L, std::vector<double> &buf_x,
    std::vector<double> &buf_vals, double thd1, double thd2);

void residual_update(const Graph<double> &G, const std::vector<double> &x,
    double dxi, int i, const std::unordered_set<int> &seedset,
    std::vector<double> &r, double gamma,
    CircularDeque<int> &Q, double kappa, std::shared_ptr<EdgeLoss> L);

/**
 * L is either TwoNormLoss or QHuberLoss, where we have
 * - `q`: the value of q in the q-norm
 * - `delta`: the value of delta in the q-Huber function
 *
 * `gamma` is for regularization,
 *   std::numeric_limits<double>::infinity returns seed set, 0 is hard/ill-posed
 * `kappa` is the sparsity regularization term.
 * `rho` is the slack term in the KKT conditions to get faster convergence
 *   (rho=1 is slow, rho=0)
 * `eps` is the value of epsilon n the local binary search
 */
std::tuple<std::vector<double>, std::vector<double>, int>
slq_diffusion(const Graph<double> &G, const std::vector<size_t> &S, double gamma,
  double kappa, double rho, std::shared_ptr<EdgeLoss> L, int max_iters=1000,
  double epsilon=1e-8);

#endif // SLQ_H
