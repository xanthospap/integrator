#ifndef __DSO_ORBIT_INTEGRATION_DOP853_HPP__
#define __DSO_ORBIT_INTEGRATION_DOP853_HPP__

#include "integration_parameters.hpp"
#include "eigen3/Eigen/Eigen"
#include <cstdint>

namespace dso {

// Function prototype for first order differential equations
// i.e. y' = f(t,y)
// where y' is output-ed in the yp vector
typedef int (*ODEfun)(double t,                       // Independent variable
                      const Eigen::Ref<const Eigen::VectorXd>& y,  // State (function values)
                      const IntegrationParameters *params,
                      Eigen::Ref<Eigen::VectorXd> yp /* Partials/Derivative */) noexcept;

class Dop853 {
private:
  ODEfun fcn;
  const IntegrationParameters *params{nullptr};
  /* number of equations in the system */
  int neqn;
  /* relative error tollerance */
  Eigen::VectorXd rtol;
  /* absolute error tolerance */
  Eigen::VectorXd atol;

  /* parameters for stiffness detection: */
  /* parameter: stifflag
   * if 0: then no stiffness testing is performed. 
   * If > 0, then we are performing stiffness testing whenever the following 
   * condition holds:
   * naccpt % stiffflag == 0 or iasti > 0
   * where naccpt is the (current) number of accepted steps.
   *
   * See functions test_for_stiffness and nstiff.
   */
  uint16_t stifflag{100};
  /* consecutive steps that no sitffness was detected. If more than 6 are 
   * encountered (in a row, assuming stifflag>0) then iasti will be set
   * to false.
   */
  int8_t nonsti{0};
  /* 'stiffness encountered' flag */
  int8_t iasti{0};

  /* max allowed steps (rejected+accepted) per integration interval */
  int max_allowed_steps{1000};

  /* step size control and stabilization parameters
   * 1) 0 <= beta < 0.4
   * fac1 and fac2 define bounds for the allowed step size change ratio:
   * fac1: is the minimum factor by which the step size can decrease, default
   *       value = .333
   * fac2: is the maximum factor by which the step size can increase, default 
   *       value is 6.0
   * safe: scale factor to reduce the aggressiveness of step size changes
   * 1e-4 <= safe < 1
   */
  double beta{0e0};
  double fac1{333e-3};
  double fac2{6e1};
  double safe{9e-1};

  /* scratch space: A matrix of neqn rows and 12 columns (used e.g. for 
   * storing the k's vectors. do as you please with them.
   */
  Eigen::Matrix<double, Eigen::Dynamic, 12> wp;

  /** @brief Are we checking for stiffness */
  bool test_for_stiffness() const noexcept { return stifflag; }

  /** @brief Return the number N, which controlls stiffness testing in the 
   * sense that: we are checking for stiffness (at least) every N accepted 
   * steps during integration .
   */
  int nstiff() const noexcept { return (int)stifflag; }

  /** @brief Core integrator for the DOP853 method.
   * 
   * This function performs the integration of a system of ordinary 
   * differential equations (ODEs) using the DOP853 explicit Runge–Kutta 
   * method. The method features adaptive step size control, error estimation, 
   * and optional stiffness detection. It operates on an already initialized 
   * system, stepping through the integration until the specified endpoint 
   * is reached.
   *
   * @param[in] t       Initial value of the independent variable.
   * @param[in] tend    Final/target value of the independent variable.
   * @param[inout] y    At input, initial values of the dependent variables.
   *                    At (successefult) output, holds the computed values of 
   *                    the dependent variable, i.e. the state vector 
   *                    corresponding to tend.
   * @param[in] hmax    Maximum value of step size.
   * @param[in] hinit   Initial step size; if set to 0, then the hinit853
   *                    will be used to compute it.
   * @return            Anything other than 0 denotes an error.
   *
   */
   int dp86co(double t, double tend, double hmax, double hinit,
             Eigen::VectorXd &y) noexcept;

  /** @brief Compute an initial step size for the DOP853 integration method.
   *
   * This function calculates a suitable initial step size `h` for use in the 
   * DOP853 integration method. The step size is based on the norms of the 
   * initial state vector `y` and its derivatives `f0`.
   *
   * Note that the function will overwrite the seconds and third columns of 
   * the scratch matrix wp (i.e. wp.col(1) and wp.col(2)), so nothing of 
   * (further) used should be placed there before a call to the function.
   *
   * @param[in] t       Initial value of the independent variable.
   * @param[in] hmax    Maximum value of step size.
   * @param[in] posneg  Indicates the direction of integration:
   *                    +1 for forward integration, 
   *                    -1 for backward integration.
   * @param[in] y0      Initial values of the dependent variables.
   * @param[in] f0      Derivatives at the initial point (dy/dx = f(x, y)).
   * @param[out] h      The computed initial step size `h`.
   * @return            Anything other than 0 denotes an error.
   *
   * @note If the computed `h` exceeds the maximum step size (`hmax`), it will 
   *       be capped at `hmax`. The sign of `h` is adjusted to match the 
   *       direction of integration (`posneg`).
   */
  int hinit853(double t, double hmax, double posneg,
                          double &h, const Eigen::VectorXd &y0,
                          const Eigen::VectorXd &f0) noexcept;

public:
  Dop853(ODEfun f, int neq, const IntegrationParameters *p, double rt,
         double at) noexcept
      : fcn{f}, params(p), neqn(neq) {
    rtol = Eigen::VectorXd::Constant(neqn, rt);
    atol = Eigen::VectorXd::Constant(neqn, at);
    wp = Eigen::Matrix<double, Eigen::Dynamic, 12>::Zero(neqn, 12);
  }

  /** @brief Set stiffness testing option.
   *
   * If set to 0, then no stiffness testing is performed. If set to a value
   * > 0, then stiffness testing will be performed in either of two cases:
   *  1. Current step is a multiple of stifflag (i.e. test every stifflag
   *      steps), or
   *  2. Stiffness already detected in previous step, and no 6-consecutive
   *     steps where taken, clear of stiffness.
   */
  void set_stiffness_check(int s) noexcept { stifflag = s; }

}; /* Dop853 */

} /* namespace dso */

#endif
