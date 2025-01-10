#include "dop853.hpp"
#include "integration_parameters.hpp"
#include <cstdio>
#include <cmath>
#include <cassert>

/* ODE equation(s): 
 * dy1/dx = -2*y1 +   y2
 * dy2/dx =    y1 - 2*y2
 *
 * Initial Conditions:
 * y1(0) = 1
 * y2(0) = 0
 *
 * Exact solution:
 * y1(x) = 1/2 e^(-x) + 1/2 e^(-3x)
 * y2(x) = 1/2 e^(-x) - 1/2 e^(-3x)
 *
 * Integrate ove intervals:
 * a. [0,5]
 */
int deriv([[maybe_unused]]double tsec, Eigen::Ref<const Eigen::VectorXd> y0,
           [[maybe_unused]]dso::IntegrationParameters *params,
           Eigen::Ref<Eigen::VectorXd> y) noexcept {
  y(0) = -2. * y0(0) + y0(1);
  y(1) = y0(0) - 2. * y0(1);
  return 0;
}

double exact_y1(double x) { return .5 * std::exp(-x) + .5 * std::exp(-3.*x); }
double exact_y2(double x) { return .5 * std::exp(-x) - .5 * std::exp(-3.*x); }

int main() {
  dso::Dop853 dop853(deriv, 2, nullptr, 1e-9, 1e-12);
  dop853.set_stiffness_check(10);

  /* initial conditions */
  Eigen::VectorXd y0(2);
  y0 << 1e0, 0e0;
  Eigen::VectorXd y(2);
  
  /* Test case 1: integrate [0,5] */
  const double xend = 5;
  if (dop853.integrate(0e0, xend, y0, y, 1e-1)) {
    fprintf(stderr, "ERROR. Integration failed!\n");
    return 1;
  }

  /* */
  printf("Exact solution at x=5: y = (%.6e, %.6e)\n", exact_y1(xend),
         exact_y2(xend));
  printf("DOP853         at x=5: y = (%.6e, %.6e) dy=(%.1e, %.1e)\n", y(0),
         y(1), std::abs(y(0) - exact_y1(xend)),
         std::abs(y(1) - exact_y2(xend)));
  assert(std::abs(y(0) - exact_y1(xend)) < 1e-15);
  assert(std::abs(y(1) - exact_y2(xend)) < 1e-15);

  return 0;
}
