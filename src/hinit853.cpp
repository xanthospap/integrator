#include "dop853.hpp"
#include <cmath>

int dso::Dop853::hinit853(double t, double hmax, double posneg, double &h,
                          const Eigen::Ref<const Eigen::VectorXd> y,
                          const Eigen::Ref<const Eigen::VectorXd> f0) noexcept {
  /* first guess for initial step */
  Eigen::ArrayXd sk = atol.array() + rtol.array() * y.array().abs();
  // printf("\thinit853: y= [");
  // for (int i = 0; i < 6; i++)
  //   printf("%.6e ", y(i));
  // printf("]\n\thinit853: dy=[");
  // for (int i = 0; i < 6; i++)
  //   printf("%.6e ", f0(i));
  // printf("]\n\thinit853: rtol=[");
  // for (int i = 0; i < 6; i++)
  //   printf("%.3e ", rtol(i));
  // printf("]\n\thinit853: atol=[");
  // for (int i = 0; i < 6; i++)
  //   printf("%.3e ", atol(i));
  // printf("]\n\thinit853: sk=[");
  // for (int i = 0; i < 6; i++)
  //   printf("%.3e ", sk(i));
  const double dnf = ((f0.array() / sk).square()).sum();
  const double dny = ((y.array() / sk).square()).sum();
  // printf("]\n");

  if ((dnf <= 1e-10) || (dny <= 1e-10)) {
    h = 1.0e-6;
  } else {
    h = std::sqrt(dny / dnf) * 1e-2;
  }

  h = std::min(h, hmax) * posneg;
  // printf("\thinit853: h(guess)=%.9f, dnf=%.3e dny=%.3e\n", h, dnf, dny);

  /* Perform an explicit Euler step */
  wp.col(1) = y + h * f0;
  if (fcn(t + h, wp.col(1), params, wp.col(2))) {
    fprintf(stderr, "[ERROR] Failed to compute partials (traceback: %s)\n",
            __func__);
    return 1;
  }

  /* Estimate the second derivative of the solution */
  double der2 = (((wp.col(2) - f0).array() / sk).square()).sum();
  der2 = std::sqrt(der2) / h;

  /* final, initial step */
  const double der12 = std::max(std::abs(der2), std::sqrt(dnf));

  const double h1 = (der12 <= 1e-15) ? std::max(1e-6, std::abs(h) * 1e-3)
                                     : std::pow(1e-2 / der12, 1. / 8.);
  h = std::min(1e2 * std::abs(h), std::min(h1, hmax)) * posneg;

  return 0;
}
