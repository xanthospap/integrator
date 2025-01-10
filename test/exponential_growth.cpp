#include "dop853.hpp"
#include "integration_parameters.hpp"
#include <cstdio>
#include <cmath>
#include <cassert>

/* ODE equation: dy/dx = y with initial values y(0) = 1
 * Exact solution: y(x) = exp(x)
 * Integrate ove intervals:
 * a. [0,1]
 * b. [0,2]
 */
int deriv([[maybe_unused]]double tsec, Eigen::Ref<const Eigen::VectorXd> y0,
           [[maybe_unused]]dso::IntegrationParameters *params,
           Eigen::Ref<Eigen::VectorXd> y) noexcept {
  y = y0;
  return 0;
}

int main() {
  dso::Dop853 dop853(deriv, 1, nullptr, 1e-8, 1e-8);
  dop853.set_stiffness_check(10);

  /* initial conditions */
  Eigen::VectorXd y0(1);
  y0 << 1e0;
  Eigen::VectorXd y(1);
  
  /* Test case 1: integrate [0,1] */
  if (dop853.integrate(0e0, 1e0, y0, y, 1e-1)) {
    fprintf(stderr, "ERROR. Integration failed!\n");
    return 1;
  }

  /* */
  printf("Exact solution at x=%.3f: y = %.6e\n", 1e0, std::exp(1e0));
  printf("DOP853         at x=%.3f: y = %.6e dy=%.1e\n", 1e0, y(0), std::abs(y(0)-std::exp(1e0)));
  assert(std::abs(y(0)-std::exp(1e0)) < 1e-15);
  
  /* integrate [0,2] */
  if (dop853.integrate(0e0, 2e0, y0, y, 1e-1)) {
    fprintf(stderr, "ERROR. Integration failed!\n");
    return 1;
  }

  /* */
  printf("Exact solution at x=%.3f: y = %.6e\n", 2e0, std::exp(2e0));
  printf("DOP853         at x=%.3f: y = %.6e dy=%.1e\n", 2e0, y(0), std::abs(y(0)-std::exp(2e0)));
  assert(std::abs(y(0)-std::exp(2e0)) < 1e-13);

  return 0;
}
