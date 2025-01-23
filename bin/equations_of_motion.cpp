#include "dop853.hpp"
#include "iers/earth_rotation.hpp"
#include "iers/fundarg.hpp"
#include "iers/gravity.hpp"
#include "iers/iau.hpp"
#include "iers/icgemio.hpp"
#include "integration_parameters.hpp"
#include "sp3.hpp"
#include "geodesy/units.hpp"
#include "yaml-cpp/yaml.h"
#include <cstdio>

constexpr const double GM_Moon = /*0.49028010560e13;*/4902.800076e9;
constexpr const double GM_Sun = /*1.32712442076e20;*/132712440040.944e9;

int gcrf2ecef(const dso::MjdEpoch &tai, dso::EopSeries &eops,
              Eigen::Matrix<double, 3, 3> &R,
              Eigen::Matrix<double, 3, 3> &dRdt,
              double *fargs,
              dso::EopRecord &eopr) noexcept {
  const auto tt = tai.tai2tt();
  double X, Y;

  /* compute (X, Y)_{cip} and (14) fundamental arguments */
  dso::xycip06a(tt, X, Y, fargs);

  /* get (interpolate EOPs) */
  if (dso::EopSeries::out_of_bounds(eops.interpolate(tt, eopr))) {
    fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
    return 1;
  }

  /* compute gmst using an approximate value for UT1 (linear interpolation) */
  double dut1_approx;
  eops.approx_dut1(tt, dut1_approx);
  [[maybe_unused]]const double gmst = dso::gmst(tt, tt.tt2ut1(dut1_approx));

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

  R = dso::detail::gcrs2itrs(dso::era00(tt.tt2ut1(eopr.dut())),
                             dso::s06(tt, X, Y), dso::sp00(tt), X, Y,
                             dso::sec2rad(eopr.xp()), dso::sec2rad(eopr.yp()),
                             eopr.lod(), dRdt);
  return 0;
}

int deriv(double tsec, Eigen::Ref<const Eigen::VectorXd> y0,
          dso::IntegrationParameters *params,
          Eigen::Ref<Eigen::VectorXd> y) noexcept {

  /* current time in TAI */
  dso::MjdEpoch t = params->t0();
  t.add_seconds(dso::FractionalSeconds(tsec));

  /* Dealunay args (14) */
  double fargs[14];

  /* Celestial to Terrestrial Matrix */
  dso::EopRecord eopr;
  Eigen::Matrix<double, 3, 3> R, dRdt;
  if (gcrf2ecef(dso::MjdEpoch(t), *(params->eops()), R, dRdt, fargs, eopr)) {
    return 8;
  }

  /* Gravity field stokes coefficients */
  auto stokes = *(params->grav());

  /* compute acceleration for given epoch/position (ECEF) */
  [[maybe_unused]] Eigen::Matrix<double, 3, 3> g;
  const Eigen::VectorXd r_ecef = R * y0.segment<3>(0);
  Eigen::Matrix<double, 3, 1> acc;
  if (dso::sh2gradient_cunningham(stokes, r_ecef, acc, g,
                                  stokes.max_degree(), stokes.max_order(), -1,
                                  -1, &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }

  // TODO adding sun worsens results !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Eigen::Matrix<double, 3, 1> acc_moon;
  Eigen::Matrix<double, 3, 1> acc_sun;
  /* get Sun position in ICRF */
  Eigen::Matrix<double, 3, 1> rtb_sun;
  if (dso::planet_pos(dso::Planet::SUN, t.tai2tt(), rtb_sun)) {
    fprintf(stderr, "ERROR Failed to compute Sun position!\n");
    return 2;
  }
  acc_sun = dso::point_mass_acceleration(y0.segment<3>(0), rtb_sun, GM_Sun);

  /* get Moon position in ICRF */
  Eigen::Matrix<double, 3, 1> rtb_moon;
  if (dso::planet_pos(dso::Planet::MOON, t.tai2tt(), rtb_moon)) {
    fprintf(stderr, "ERROR Failed to compute Moon position!\n");
    return 2;
  }
  acc_moon = dso::point_mass_acceleration(y0.segment<3>(0), rtb_moon, GM_Moon);

  /* Solid Earth Tide (ITRF) */
  Eigen::Matrix<double, 3, 1> acc_set;
  params->mse_tide->stokes_coeffs(t.tai2tt(), t.tai2ut1(eopr.dut()),
                                  R * rtb_moon, R * rtb_sun, fargs);
  if (dso::sh2gradient_cunningham(params->mse_tide->stokes_coeffs(), r_ecef,
                                  acc_set, g, -1, -1, -1, -1, &(params->tw()),
                                  &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }

  {
    // print for debugging
    const auto tgps = t.tai2gps();
    printf("%.15f %.15f %.15f %.15f\n", tgps.as_mjd(), acc_moon(0), acc_moon(1), acc_moon(2));
    //printf("%.15f %.15f %.15f %.15f\n", tgps.as_mjd(), acc_sun(0), acc_sun(1), acc_sun(2));
    //printf("%.15f %.15f %.15f %.15f\n", tgps.as_mjd(), acc_set(0), acc_set(1), acc_set(2));
  }

  /* set velocity vector (ICRF) */
  y.segment<3>(0) = y0.segment<3>(3);

  /* ECEF to ICRF note that y = (v, a) and y0 = (r, v) */
  //y.segment<3>(3) =
  //    R.transpose() * acc + acc_tb;
  y.segment<3>(3) = (acc_moon + acc_sun);
  y.segment<3>(3) += R.transpose() * acc_set;
  y.segment<3>(3) += R.transpose() * acc;

  return 0;
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s CONFIG SP3_file [SAT_ID]\n", argv[0]);
    return 1;
  }

  /* parse the yaml configuration file */
  YAML::Node config = YAML::LoadFile(argv[1]);

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

  /* get starting epoch in TAI */
  auto start_t = sp3.start_epoch();
  if (!std::strcmp(sp3.time_sys(), "GPS")) {
    start_t = start_t.gps2tai();
  }

  /* EOPs */
  dso::EopSeries eop;
  {
    std::string tmp = config["eop"].as<std::string>();
    auto t1 = start_t;
    t1.add_seconds(dso::seconds(-86400));
    auto t2 = start_t;
    t2.add_seconds(dso::seconds(5 * 86400));
    if (dso::parse_iers_C04(tmp.c_str(), dso::MjdEpoch(t1), dso::MjdEpoch(t2),
                            eop)) {
      fprintf(stderr, "ERROR Failed parsing eop file\n");
      return 1;
    }
    eop.regularize();
  }

  /* load planetary ephemeris kernels */
  std::string tmp = config["planetary-ephemeris"]["bsp"].as<std::string>();
  dso::load_spice_kernel(tmp.c_str());
  tmp = config["planetary-ephemeris"]["tls"].as<std::string>();
  dso::load_spice_kernel(tmp.c_str());

  /* Earth's gravity field */
  dso::StokesCoeffs stokes;
  {
    tmp = config["gravity"]["model"].as<std::string>();
    int DEGREE = config["gravity"]["degree"].as<int>();
    int ORDER = config["gravity"]["order"].as<int>();
    dso::Icgem icgem(tmp.c_str());
    if (icgem.parse_data(DEGREE, ORDER, start_t, stokes)) {
      fprintf(stderr, "ERROR Failed reading gravity model!\n");
      return 1;
    }
    /* checks */
    assert(stokes.max_degree() == DEGREE);
    assert(stokes.max_order() == ORDER);
  }

  /* Solid Earth Tides */
  dso::SolidEarthTide setide(iers2010::GMe, iers2010::Re, GM_Sun, GM_Moon);

  /* setup integration parameters */
  dso::IntegrationParameters params;
  params.meops = &eop;
  params.mgrav = &stokes;
  params.mse_tide = &setide;
  params.mtai0 = dso::MjdEpoch(start_t);

  /* setup the integrator */
  dso::Dop853 dop853(deriv, 6, &params, 1e-9, 1e-12);
  dop853.set_stiffness_check(10);

  /* dummy */
  double fargs[14];
  dso::EopRecord eopr;

  /* Just for testing Vs costg */
  Eigen::VectorXd state = Eigen::Matrix<double, 6, 1>::Zero();
  Eigen::VectorXd y = Eigen::Matrix<double, 6, 1>::Zero();
  Eigen::Matrix<double, 3, 3> R, dRdt;
  std::size_t it = 0;
  dso::Sp3DataBlock block;
  int sp3err = 0;
  while (!sp3err) {
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
    /* GCRF to ITRF rotation matrix */
    if (gcrf2ecef(dso::MjdEpoch(block_tai), eop, R, dRdt, fargs, eopr)) {
      return 8;
    }
    /* get state for current epoch ITRF */
    state << block.state[0] * 1e3, block.state[1] * 1e3, block.state[2] * 1e3,
        block.state[4] * 1e-1, block.state[5] * 1e-1, block.state[6] * 1e-1;
    /* transform state to GCRF */
    y.segment<3>(0) = R.transpose() * state.segment<3>(0);
    y.segment<3>(3) = R.transpose() * state.segment<3>(3) +
                      dRdt.transpose() * state.segment<3>(0);
    state = y;
    /* seconds since initial epoch */
    dso::FractionalSeconds sec =
        block_tai.diff<dso::DateTimeDifferenceType::FractionalSeconds>(start_t);
    /* compute derivative at this epoch */
    if (deriv(sec.seconds(), state, &params, y)) {
      fprintf(stderr, "ERROR. Failed computing derivative!\n");
      return 1;
    }
    ++it;
    if (it >= 2000)
      break;
  }

  //Eigen::VectorXd state = Eigen::Matrix<double, 6, 1>::Zero();
  //Eigen::VectorXd y = Eigen::Matrix<double, 6, 1>::Zero();
  //Eigen::Matrix<double, 3, 3> R, dRdt;
  //std::size_t it = 0;
  //dso::Sp3DataBlock block;
  //int sp3err = 0;
  //while (!sp3err) {
  //  sp3err = sp3.get_next_data_block(sv, block);
  //  if (sp3err > 0) {
  //    printf("Something went wrong ....status = %3d\n", sp3err);
  //    return 1;
  //  }
  //  bool position_ok = !block.flag.is_set(dso::Sp3Event::bad_abscent_position);
  //  if (position_ok && (!sp3err)) {
  //    auto block_tai = block.t;
  //    if (!std::strcmp(sp3.time_sys(), "GPS")) {
  //      block_tai = block_tai.gps2tai();
  //    }
  //    if (gcrf2ecef(dso::MjdEpoch(block_tai), eop, R, dRdt, fargs, eopr)) {
  //      return 8;
  //    }
  //    if (!it) {
  //      /* first state of satellite in file; transform to celestial and store 
  //       * as state and start_t. This is where we start integrating from.
  //       */
  //      state << block.state[0] * 1e3, block.state[1] * 1e3,
  //          block.state[2] * 1e3, block.state[4] * 1e-1, block.state[5] * 1e-1,
  //          block.state[6] * 1e-1;
  //      y.segment<3>(0) = R.transpose() * state.segment<3>(0);
  //      y.segment<3>(3) = R.transpose() * state.segment<3>(3) +
  //                        dRdt.transpose() * state.segment<3>(0);
  //      start_t = block_tai;
  //      state = y;
  //    } else {
  //      /* new entry; seconds since intial epoch */
  //      dso::FractionalSeconds sec =
  //          block_tai.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
  //              start_t);
  //      /* integrate from initial conditions to this epoch */
  //      if (dop853.integrate(dso::MjdEpoch(start_t), sec.seconds(), state, y)) {
  //        fprintf(stderr, "ERROR. Integration failed! sec is %.3f\n",
  //                sec.seconds());
  //        return 1;
  //      }
  //      /* transform integration results to ECEF */
  //      Eigen::VectorXd yt = y;
  //      y.segment<3>(0) = R * yt.segment<3>(0);
  //      y.segment<3>(3) = R * yt.segment<3>(3) + dRdt * yt.segment<3>(0);
  //      printf("%.12f %.6f %.6f %.6f %.9f %.9f %.9f %.6f %.6f %.6f %.9f %.9f "
  //             "%.9f\n",
  //             sec.seconds(),
  //             block.state[0] * 1e3, block.state[1] * 1e3, block.state[2] * 1e3,
  //             block.state[4] * 1e-1, block.state[5] * 1e-1,
  //             block.state[6] * 1e-1, y(0), y(1), y(2), y(3), y(4), y(5));
  //    }
  //  }
  //  ++it;
  //  if (it >= 100)
  //    break;
  //}

  return sp3err;
}
