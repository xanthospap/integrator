#include "dop853.hpp"
#include <cmath>

namespace {
constexpr const double a21 = 0.0526001519587677318785587544488;
constexpr const double a31 = 0.0197250569845378994544595329183;
constexpr const double a32 = 0.0591751709536136983633785987549;
constexpr const double a41 = 0.0295875854768068491816892993775;
constexpr const double a43 = 0.0887627564304205475450678981324;
constexpr const double a51 = 0.241365134159266685502369798665;
constexpr const double a53 = -0.884549479328286085344864962717;
constexpr const double a54 = 0.924834003261792003115737966543;
constexpr const double a61 = 0.0370370370370370370370370370370;
constexpr const double a64 = 0.170828608729473871279604482173;
constexpr const double a65 = 0.125467687566822425016691814123;
constexpr const double a71 = 0.0371093750000000000000000000000;
constexpr const double a74 = 0.170252211019544039314978060272;
constexpr const double a75 = 0.0602165389804559606850219397278;
constexpr const double a76 = -0.0175781250000000000000000000000;
constexpr const double a81 = 0.0370920001185047927108779319836;
constexpr const double a84 = 0.170383925712239993810214054705;
constexpr const double a85 = 0.107262030446373284651809199168;
constexpr const double a86 = -0.0153194377486244017527936158236;
constexpr const double a87 = 0.00827378916381402288758473766002;
constexpr const double a91 = 0.624110958716075717114429577812;
constexpr const double a94 = -3.36089262944694129406857109825;
constexpr const double a95 = -0.868219346841726006818189891453;
constexpr const double a96 = 27.5920996994467083049415600797;
constexpr const double a97 = 20.1540675504778934086186788979;
constexpr const double a98 = -43.4898841810699588477366255144;
constexpr const double a101 = 0.477662536438264365890433908527;
constexpr const double a104 = -2.48811461997166764192642586468;
constexpr const double a105 = -0.590290826836842996371446475743;
constexpr const double a106 = 21.2300514481811942347288949897;
constexpr const double a107 = 15.2792336328824235832596922938;
constexpr const double a108 = -33.2882109689848629194453265587;
constexpr const double a109 = -0.0203312017085086261358222928593;
constexpr const double a111 = -0.93714243008598732571704021658;
constexpr const double a114 = 5.18637242884406370830023853209;
constexpr const double a115 = 1.09143734899672957818500254654;
constexpr const double a116 = -8.14978701074692612513997267357;
constexpr const double a117 = -18.5200656599969598641566180701;
constexpr const double a118 = 22.7394870993505042818970056734;
constexpr const double a119 = 2.49360555267965238987089396762;
constexpr const double a1110 = -3.04676447189821950038236690220;
constexpr const double a121 = 2.27331014751653820792359768449;
constexpr const double a124 = -10.5344954667372501984066689879;
constexpr const double a125 = -2.00087205822486249909675718444;
constexpr const double a126 = -17.9589318631187989172765950534;
constexpr const double a127 = 27.9488845294199600508499808837;
constexpr const double a128 = -2.85899827713502369474065508674;
constexpr const double a129 = -8.87285693353062954433549289258;
constexpr const double a1210 = 12.3605671757943030647266201528;
constexpr const double a1211 = 0.643392746015763530355970484046;
constexpr const double c2 = 0.526001519587677318785587544488e-01;
constexpr const double c3 = 0.789002279381515978178381316732e-01;
constexpr const double c4 = 0.118350341907227396726757197510e+00;
constexpr const double c5 = 0.281649658092772603273242802490e+00;
constexpr const double c6 = 0.333333333333333333333333333333e+00;
constexpr const double c7 = 0.25e+00;
constexpr const double c8 = 0.307692307692307692307692307692e+00;
constexpr const double c9 = 0.651282051282051282051282051282e+00;
constexpr const double c10 = 0.6e+00;
constexpr const double c11 = 0.857142857142857142857142857142e+00;
constexpr const double c14 = 0.1e+00;
constexpr const double c15 = 0.2e+00;
constexpr const double c16 = 0.777777777777777777777777777778e+00;
constexpr const double b1 = 5.42937341165687622380535766363e-2;
constexpr const double b6 = 4.45031289275240888144113950566e0;
constexpr const double b7 = 1.89151789931450038304281599044e0;
constexpr const double b8 = -5.8012039600105847814672114227e0;
constexpr const double b9 = 3.1116436695781989440891606237e-1;
constexpr const double b10 = -1.52160949662516078556178806805e-1;
constexpr const double b11 = 2.01365400804030348374776537501e-1;
constexpr const double b12 = 4.47106157277725905176885569043e-2;
constexpr const double bhh1 = 0.244094488188976377952755905512e+00;
constexpr const double bhh2 = 0.733846688281611857341361741547e+00;
constexpr const double bhh3 = 0.220588235294117647058823529412e-01;
constexpr const double er1 = 0.1312004499419488073250102996e-01;
constexpr const double er6 = -0.1225156446376204440720569753e+01;
constexpr const double er7 = -0.4957589496572501915214079952e+00;
constexpr const double er8 = 0.1664377182454986536961530415e+01;
constexpr const double er9 = -0.3503288487499736816886487290e+00;
constexpr const double er10 = 0.3341791187130174790297318841e+00;
constexpr const double er11 = 0.8192320648511571246570742613e-01;
constexpr const double er12 = -0.2235530786388629525884427845e-01;
} /* namespace */

int dso::Dop853::dp86co(double t, double tend, double hmax, double hinit,
                        Eigen::Ref<Eigen::VectorXd> y) noexcept {
  // printf("Called dp86co from %.9f to %.9f with initial state=(%.3f, %.3f
  // %.3f, "
  //        "%.6f %.6f %.6f)\n",
  //        t, tend, y(0), y(1), y(2), y(3), y(4), y(5));

  /* Initial factor for step size adjustment */
  const double expo1 = 1e0 / 8e0 - beta * 0.2;
  const double facc1 = 1e0 / fac1;
  const double facc2 = 1e0 / fac2;
  double facold = 1e-4;

  /* Determine the direction of integration */
  const double posneg = ((tend - t) >= 0) ? 1 : -1;

  /* number of steps taken (total, rejected, accepted) */
  int nstep{0};
  int naccpt{0};
  int nrejct{0};

  /* number of function calls */
  unsigned nfcn{0};

  /* evaluate the function at the initial point, store at wp.col(0) */
  if (fcn(t, y, params, wp.col(0))) {
    fprintf(stderr,
            "[ERROR] Integrator failed! Failed computing derivative at t0! "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  /* if initial step is zero, compute a valid value. note: hinit853 will
   * overwrite wp.col(1) and wp.col(2)
   */
  hmax = std::abs(hmax);
  double h = hinit;
  if (h == 0e0) {
    // printf("\tComputing initial step size, ...\n");
    if (hinit853(t, hmax, posneg, h, y, wp.col(0))) {
      fprintf(stderr,
              "[ERROR] Failed computing initial guess for step size! "
              "(traceback: %s)\n",
              __func__);
      return 1;
    }
    // printf("\tinitial step size, h=%.12f\n", h);
    ++nfcn;
  }

  /* basic integration step: keep takin steps (either accepted or rejected),
   * until a) we reached (or nearly reached) tend, b) reached the maximum
   * steps allowed, or c) an error occured, signalled by `error`!=0.
   */
  int last_step = 0;
  int last_step_rejected = false;
  int error = 0;
  while (((!error) && (nstep++ < max_allowed_steps)) && (!last_step)) {

    /* check if step size is too small */
    if (1e-1 * std::abs(h) <=
        std::abs(t) * std::numeric_limits<double>::epsilon()) {
      fprintf(stderr,
              "[ERROR] Integrator failed! Step size underflow! (traceback: "
              "%s)\n",
              __func__);
      return 1;
    }

    /* adjust step size for the last step */
    if ((t + 1.01 * h - tend) * posneg > 0e0) {
      h = tend - t;
      last_step = 1;
    }

    /* increment step counter (regrdless if accepted or not) */
    nstep += 1;

    /* the twelve stages (stage 1 already computed in wp.col(0)) */
    // Stage 2
    Eigen::VectorXd y1 = y + h * a21 * wp.col(0);
    error += fcn(t + c2 * h, y1, params, wp.col(1));

    // Stage 3
    y1 = y + h * (a31 * wp.col(0) + a32 * wp.col(1));
    error += fcn(t + c3 * h, y1, params, wp.col(2));

    // Stage 4
    y1 = y + h * (a41 * wp.col(0) + a43 * wp.col(2));
    error += fcn(t + c4 * h, y1, params, wp.col(3));

    // Stage 5
    y1 = y + h * (a51 * wp.col(0) + a53 * wp.col(2) + a54 * wp.col(3));
    error += fcn(t + c5 * h, y1, params, wp.col(4));

    // Stage 6
    y1 = y + h * (a61 * wp.col(0) + a64 * wp.col(3) + a65 * wp.col(4));
    error += fcn(t + c6 * h, y1, params, wp.col(5));

    // Stage 7
    y1 = y + h * (a71 * wp.col(0) + a74 * wp.col(3) + a75 * wp.col(4) +
                  a76 * wp.col(5));
    error += fcn(t + c7 * h, y1, params, wp.col(6));

    // Stage 8
    y1 = y + h * (a81 * wp.col(0) + a84 * wp.col(3) + a85 * wp.col(4) +
                  a86 * wp.col(5) + a87 * wp.col(6));
    error += fcn(t + c8 * h, y1, params, wp.col(7));

    // Stage 9
    y1 = y + h * (a91 * wp.col(0) + a94 * wp.col(3) + a95 * wp.col(4) +
                  a96 * wp.col(5) + a97 * wp.col(6) + a98 * wp.col(7));
    error += fcn(t + c9 * h, y1, params, wp.col(8));

    // Stage 10
    y1 = y + h * (a101 * wp.col(0) + a104 * wp.col(3) + a105 * wp.col(4) +
                  a106 * wp.col(5) + a107 * wp.col(6) + a108 * wp.col(7) +
                  a109 * wp.col(8));
    error += fcn(t + c10 * h, y1, params, wp.col(9));

    // Stage 11 (note that we are storing result in col 1)
    y1 = y + h * (a111 * wp.col(0) + a114 * wp.col(3) + a115 * wp.col(4) +
                  a116 * wp.col(5) + a117 * wp.col(6) + a118 * wp.col(7) +
                  a119 * wp.col(8) + a1110 * wp.col(9));
    error += fcn(t + c11 * h, y1, params, wp.col(1));

    // Stage 12 (note that we are storing result in col 2)
    y1 = y + h * (a121 * wp.col(0) + +a1211 * wp.col(1) + a124 * wp.col(3) +
                  a125 * wp.col(4) + a126 * wp.col(5) + a127 * wp.col(6) +
                  a128 * wp.col(7) + a129 * wp.col(8) + a1210 * wp.col(9));
    error += fcn(/*t + c16 * h*/ t + h, y1, params, wp.col(2));
    nfcn += 11;

    wp.col(3) = b1 * wp.col(0) + b11 * wp.col(1) + b12 * wp.col(2) +
                b6 * wp.col(5) + b7 * wp.col(6) + b8 * wp.col(7) +
                b9 * wp.col(8) + b10 * wp.col(9);
    wp.col(4) = y + h * wp.col(3);

    /* vector tolerances of size n=neqn */
    Eigen::ArrayXd sk =
        atol.array() +
        rtol.array() * (y.array().abs().max(wp.col(4).array().abs()));

    /* Compute the first error component */
    Eigen::VectorXd verr =
        wp.col(3) - bhh1 * wp.col(0) - bhh2 * wp.col(8) - bhh3 * wp.col(2);
    const double err2 = ((verr.array() / sk).square()).sum();

    /* Compute the second error component */
    verr = er1 * wp.col(0) + er11 * wp.col(1) + er12 * wp.col(2) +
           er6 * wp.col(5) + er7 * wp.col(6) + er8 * wp.col(7) +
           er9 * wp.col(8) + er10 * wp.col(9);
    const double err1 = ((verr.array() / sk).square()).sum();

    /* Combine the error components into a denominator */
    double deno = err1 + 1e-2 * err2;
    if (deno <= 0e0)
      deno = 1e0;

    /* Final error estimate */
    const double err = std::abs(h) * err1 * std::sqrt(1e0 / (neqn * deno));

    /* Compute factor for step size adjustment */
    const double fac11 = std::pow(err, expo1);
    /* Lund stabilization */
    const double _fac = fac11 / std::pow(facold, beta);
    /* bound the scaling factor */
    const double fac = std::max(facc2, std::min(facc1, _fac / safe));
    /* new step size */
    double hnew = h / fac;

    if (err <= 1e0) {
      /* Step is accepted */
      facold = std::max(err, 1e-4);
      /* increment count of accepted steps */
      naccpt += 1;

      /* recompute k4 using the updated step */
      error += fcn(t + h, wp.col(4), params, wp.col(3));
      /* increment function call counter */
      ++nfcn;

      /* stiffness detection */
      if (test_for_stiffness() && ((naccpt % nstiff() == 0) || (iasti > 0))) {
        const double stnum = (wp.col(3) - wp.col(2)).array().square().sum();
        const double stden = (wp.col(4) - y1).array().square().sum();

        if (stden > 0e0) {
          /* stiffness indicator */
          const double hlamb = std::abs(h) * std::sqrt(stnum / stden);
          if (hlamb > 6.1) {
            /* set non-stiff flag/counter */
            nonsti = 0;
            /* increment stiffness counter */
            iasti += 1;
            if (iasti == 15) {
              /* trigger stiffness warning and give up */
              fprintf(stderr,
                      "[WRNNG] Integrator possibly encoutered stiff problem! "
                      "(traceback: %s)\n",
                      __func__);
              error += -100;
            }
          } else {
            /* increment non-stiff counter */
            nonsti += 1;
            if (nonsti == 6)
              /* reset stiffness counter */
              iasti = 0;
          }
        }
      } /* end stiffness detection */

      /* step accepted: update for the next iteration */
      t += h;

      /* step accepted: update the state vector */
      wp.col(0) = wp.col(3);
      y = wp.col(4);

      /* next step size */
      if (std::abs(hnew) > hmax)
        hnew = posneg * hmax;
      if (last_step_rejected)
        hnew = posneg * std::min(std::abs(hnew), std::abs(h));
      last_step_rejected = 0;
      // printf("step accepted, new step size = %.9f\n", hnew);

    } else {
      /* Step rejected */
      nrejct += 1;
      hnew = h / std::min(facc1, fac11 / safe);
      last_step_rejected = 1;
      last_step = 0;
      // printf("step rejected, new step size = %.9f\n", hnew);
    }

    /* update step size for next iteration. note that if the step was
     * rejected, t will have not changed, hence we are going forward, from
     * the same starting t, with a new step size (h).
     */
    h = hnew;

  } /* main loop (while) */

  /* Check if maximum steps are exceeded */
  if (nstep >= max_allowed_steps) {
    fprintf(stderr,
            "[ERROR] Integrator failed! Too many number of steps! "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  /* error checking: problem is stiff! */
  if (error == -100) {
    fprintf(stderr,
            "[ERROR] Integrator possibly encoutered stiff problem! "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  /* error in the calling function */
  if (error) {
    fprintf(stderr,
            "[ERROR] Integrator failed calling the derivative function "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  /* if everything is ok, we should have reached the last step; y now contains
   * the state at tend
   */
  printf("\tNumber of deriv calls in dp86 is %u\n", nfcn);
  return (!last_step);
}
