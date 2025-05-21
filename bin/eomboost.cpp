#include "datetime/tpdate.hpp"
#include "geodesy/units.hpp"
#include "iers/earth_rotation.hpp"
#include "iers/fundarg.hpp"
#include "iers/gravity.hpp"
#include "iers/iau.hpp"
#include "iers/icgemio.hpp"
#include "iers/ocean_tide.hpp"
#include "iers/relativity.hpp"
#include "integration_parameters.hpp"
#include "sp3.hpp"
#include "yaml-cpp/yaml.h"
#include <boost/numeric/odeint.hpp>
#include <cstdio>
#include <iostream> // DEBUG

constexpr const double GM_Moon = /*0.49028010560e13;*/ 4902.800076e9;
constexpr const double GM_Sun = /*1.32712442076e20;*/ 132712440040.944e9;

/* Define the State Type for ODEINT */
typedef std::array<double, 6> state_type;
typedef boost::numeric::odeint::runge_kutta_dopri5<state_type> dopri5_type;

/* Compute relevant quaternions for the ITRS/GCRS transformation
 */
int prep_c2i(const dso::MjdEpoch &tai, dso::EopSeries &eops,
             Eigen::Quaterniond &q_c2tirs, Eigen::Quaterniond &q_tirs2i,
             double *fargs, dso::EopRecord &eopr) noexcept {

  /* epoch of request in TT */
  const auto tt = tai.tai2tt();
  // printf("\tPreparing rotations for TAI=%.15f\n", tai.as_mjd());

  /* compute (X,Y) CIP and fundamental arguments (we are doing this here
   * to compute fargs).
   */
  double Xcip, Ycip;
  dso::xycip06a(tt, Xcip, Ycip, fargs);

  /* interpolate EOPs */
  if (dso::EopSeries::out_of_bounds(eops.interpolate(tt, eopr))) {
    fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
    return 1;
  }

  /* compute gmst using an approximate value for UT1 (linear interpolation) */
  double dut1_approx;
  eops.approx_dut1(tt, dut1_approx);
  const double gmst = dso::gmst(tt, tt.tt2ut1(dut1_approx));

  /* add libration effect [micro as] */
  {
    double dxp, dyp, dut1, dlod;
    dso::deop_libration(fargs, gmst, dxp, dyp, dut1, dlod);
    eopr.xp() += dxp * 1e-6;
    eopr.yp() += dyp * 1e-6;
    eopr.dut() += dut1 * 1e-6;
    eopr.lod() += dlod * 1e-6;
  }

  /* add ocean tidal effect [micro as] */
  {
    double dxp, dyp, dut1, dlod;
    dso::deop_ocean_tide(fargs, gmst, dxp, dyp, dut1, dlod);
    eopr.xp() += dxp * 1e-6;
    eopr.yp() += dyp * 1e-6;
    eopr.dut() += dut1 * 1e-6;
    eopr.lod() += dlod * 1e-6;
  }

  /* de-regularize */
  {
    double ut1_cor;
    double lod_cor;
    double omega_cor;
    dso::deop_zonal_tide(fargs, ut1_cor, lod_cor, omega_cor);
    /* apply (note: microseconds to seconds) */
    eopr.dut() += (ut1_cor * 1e-6);
    eopr.lod() += (lod_cor * 1e-6);
  }

  /* use fundamental arguments to compute s */
  const double s = dso::s06(tt, Xcip, Ycip, fargs);

  /* apply CIP corrections */
  Xcip += dso::sec2rad(eopr.dX());
  Ycip += dso::sec2rad(eopr.dY());

  /* spherical crd for CIP (E, d) */
  double d, e;
  dso::detail::xycip2spherical(Xcip, Ycip, d, e);

  /* Earth rotation angle */
  const double era = dso::era00(tt.tt2ut1(eopr.dut()));

  /* compute rotation quaternions */
  q_c2tirs = dso::detail::c2tirs(era, s, d, e);
  q_tirs2i = dso::detail::tirs2i(dso::sec2rad(eopr.xp()),
                                 dso::sec2rad(eopr.yp()), dso::sp00(tt));

  return 0;
}

struct EomSystem {
  dso::IntegrationParameters *params;
  EomSystem(dso::IntegrationParameters *p) noexcept : params(p) {}

  void operator()(const state_type &x, state_type &dxdt, double t) {
    Eigen::Matrix<double, 6, 1> y0;
    for (int i = 0; i < 6; i++)
      y0(i) = x[i];

    /* epoch of request in TT */
    const auto tt =
        (params->t0().add_seconds(dso::FractionalSeconds(t))).tai2tt();

    /* accumulated acceleration and gradient in ITRS */
    Eigen::Vector3d ai = Eigen::Vector3d::Zero();
    [[maybe_unused]] Eigen::Matrix<double, 3, 3> gi;
    /* accumulated acceleration in GCRS */
    Eigen::Vector3d ac = Eigen::Vector3d::Zero();

    /* GCRS/ITRS at t */
    double fargs[14];
    dso::EopRecord eopr;
    Eigen::Quaterniond q_c2tirs, q_tirs2i;
    prep_c2i(tt.tt2tai(), params->eops(), q_c2tirs, q_tirs2i, fargs, eopr);
    Eigen::Vector3d omega;
    omega << 0e0, 0e0, dso::earth_rotation_rate(eopr.lod());

    /* state in ITRS (from GCRS) */
    Eigen::Matrix<double, 6, 1> itrs = Eigen::Matrix<double, 6, 1>::Zero();
    itrs.segment<3>(0) = q_tirs2i * (q_c2tirs * y0.segment<3>(0));
    itrs.segment<3>(3) = q_tirs2i * (q_c2tirs * y0.segment<3>(3) -
                                     omega.cross(q_c2tirs * y0.segment<3>(0)));

    /* gravity acceleration */
    if (dso::sh2gradient_cunningham(params->earth_gravity(), itrs.segment<3>(0),
                                    ai, gi,
                                    params->earth_gravity().max_degree(),
                                    params->earth_gravity().max_order(), -1, -1,
                                    &(params->tw()), &(params->tm()))) {
      fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
      assert(1 == 2);
    }

    /* Third Body perturbations */
    {
      /* get Sun position & velocity in ICRF */
      Eigen::Matrix<double, 6, 1> rtb_sun;
      if (dso::planet_state(dso::Planet::SUN, tt, rtb_sun)) {
        fprintf(stderr, "ERROR Failed to compute Sun position!\n");
        assert(1 == 2);
      }
      ac += dso::point_mass_acceleration(y0.segment<3>(0),
                                         rtb_sun.segment<3>(0), GM_Sun);

      /* get Moon position in ICRF */
      Eigen::Matrix<double, 3, 1> rtb_moon;
      if (dso::planet_pos(dso::Planet::MOON, tt, rtb_moon)) {
        fprintf(stderr, "ERROR Failed to compute Moon position!\n");
        assert(1 == 2);
      }
      ac += dso::point_mass_acceleration(y0.segment<3>(0), rtb_moon, GM_Moon);
    }

    /* form the derivative vector */
    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];
    // yp.segment<3>(0) = y0.segment<3>(3);
    Eigen::Vector3d a = ac + q_c2tirs.conjugate() * (q_tirs2i.conjugate() * ai);
    dxdt[3] = a(0);
    dxdt[4] = a(1);
    dxdt[5] = a(2);

    return;
  }
}; /* struct EomSystem */

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s [CONFIG] [SP3] \n", argv[0]);
    return 1;
  }

  /* create an Sp3c instance */
  dso::Sp3c sp3(argv[2]);

  /* choose satellite */
  dso::sp3::SatelliteId sv = sp3.sattellite_vector()[0];
  if (argc > 3) {
    /* check if the satellite is included in the Sp3 */
    if (!sp3.has_sv(dso::sp3::SatelliteId(argv[3]))) {
      fprintf(stderr, "Error. Satellite [%s] not included in sp3 file!\n",
              argv[3]);
      return 1;
    }
    sv = dso::sp3::SatelliteId{argv[3]};
  }

  /* get starting epoch in TT */
  auto start_t = sp3.start_epoch();
  if (!std::strcmp(sp3.time_sys(), "GPS")) {
    start_t = start_t.gps2tai();
  }

  /* load parameters from YAML */
  auto t2 = start_t.add_seconds(dso::seconds(5 * 86400));
  dso::IntegrationParameters params = dso::IntegrationParameters::from_config(
      argv[1], dso::MjdEpoch(start_t), dso::MjdEpoch(t2));

  /* dummy */
  double fargs[14];
  dso::EopRecord eopr;

  /* integrator */
  auto stepper = make_controlled(1e-8, 1e-8, dopri5_type());

  Eigen::VectorXd y0 = Eigen::Matrix<double, 6, 1>::Zero();
  Eigen::VectorXd y = Eigen::Matrix<double, 6, 1>::Zero();
  Eigen::Quaterniond q_c2tirs, q_tirs2i;
  Eigen::Vector3d omega;
  dso::MjdEpoch t0;
  std::size_t it = 0;
  dso::Sp3DataBlock block;
  int sp3err = 0;

  while (!sp3err) {
    /* get next record from sp3 */
    sp3err = sp3.get_next_data_block(sv, block);
    if (sp3err > 0) {
      printf("Something went wrong ....status = %3d\n", sp3err);
      return 1;
    }

    /* time of current block in TAI */
    auto block_tai = block.t;
    if (!std::strcmp(sp3.time_sys(), "GPS")) {
      block_tai = block_tai.gps2tai();
    }

    /* current time, TAI */
    const auto tai = dso::MjdEpoch(block_tai);
    // printf("-> Current time is %.12f (TAI)\n", tai.as_mjd());

    /* GCRS/ITRS at t */
    prep_c2i(tai, params.eops(), q_c2tirs, q_tirs2i, fargs, eopr);
    omega << 0e0, 0e0, dso::earth_rotation_rate(eopr.lod());

    /* parse state for current epoch, ITRS */
    y << block.state[0] * 1e3, block.state[1] * 1e3, block.state[2] * 1e3,
        block.state[4] * 1e-1, block.state[5] * 1e-1, block.state[6] * 1e-1;

    if (!it) {
      /* transform state to GCRS (from ITRS) */
      Eigen::VectorXd yc = Eigen::Matrix<double, 6, 1>::Zero();
      yc.segment<3>(0) =
          q_c2tirs.conjugate() * (q_tirs2i.conjugate() * y.segment<3>(0));
      yc.segment<3>(3) = q_c2tirs.conjugate() *
                         (q_tirs2i.conjugate() * y.segment<3>(3) +
                          omega.cross(q_tirs2i.conjugate() * y.segment<3>(0)));
      y0 = yc;
      t0 = tai;
      params.t0() = tai;
    } else {
      /* seconds since initial epoch */
      dso::FractionalSeconds sec =
          tai.diff<dso::DateTimeDifferenceType::FractionalSeconds>(t0);

      /* integrate from first SP3 record to here - setup initial conditions */
      params.t0() = t0;
      state_type x = {y0(0), y0(1), y0(2), y0(3), y0(4), y0(5)};
      boost::numeric::odeint::integrate_adaptive(stepper, EomSystem(&params), x,
                                                 0e0, sec.seconds(), 1e-3);
      Eigen::VectorXd yc = Eigen::Matrix<double, 6, 1>::Zero();
      yc << x[0], x[1], x[2], x[3], x[4], x[5];

      /* compare in ITRS */
      Eigen::VectorXd yi = Eigen::Matrix<double, 6, 1>::Zero();
      yi.segment<3>(0) = q_tirs2i * (q_c2tirs * yc.segment<3>(0));
      yi.segment<3>(3) = q_tirs2i * (q_c2tirs * yc.segment<3>(3) -
                                     omega.cross(q_c2tirs * yc.segment<3>(0)));

      printf("%.9f %.3f %.3f %.3f %.6f %.6f %.6f\n", sec.seconds(),
             yi(0) - y(0), yi(1) - y(1), yi(2) - y(2), yi(3) - y(3),
             yi(4) - y(4), yi(5) - y(5));
    }

    ++it;
    if (it > 5)
      break;
  }

  printf("# Num of epochs parsed/used: %ld\n", it);
}