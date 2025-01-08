#include "dop853.hpp"
#include "iers/earth_rotation.hpp"
#include "iers/fundarg.hpp"
#include "iers/gravity.hpp"
#include "iers/iau.hpp"
#include "iers/icgemio.hpp"
#include "integration_parameters.hpp"
#include "sp3.hpp"
#include "yaml-cpp/yaml.h"
#include <cstdio>

int gcrf2ecef(const dso::MjdEpoch &tai, dso::EopSeries &eops,
              Eigen::Matrix<double, 3, 3> &R,
              Eigen::Matrix<double, 3, 3> &dRdt) noexcept {
  double fargs[14];
  dso::EopRecord eopr;
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
  const double gmst = dso::gmst(tt, tt.tt2ut1(dut1_approx));

  /* compute fundamental arguments at given epoch */
  // dso::fundarg(tt, fargs);

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
                             dso::s06(tt, X, Y), dso::sp00(tt), X, Y, eopr.xp(),
                             eopr.yp(), eopr.lod(), dRdt);
  // Eigen::Matrix<double, 3, 3> R = R1.transpose();
  // Eigen::Matrix<double, 3, 3> dRdt = R1.transpose();
  // yi.segment<3>(0) = R1 * ye.segment<3>(0);
  // yi.segment<3>(3) = R1 * ye.segment<3>(3) + dRdt * ye.segment<3>(0);

  return 0;
}

int deriv(double tsec, Eigen::Ref<const Eigen::VectorXd> y0,
          dso::IntegrationParameters *params,
          Eigen::Ref<Eigen::VectorXd> y) noexcept {

  dso::MjdEpoch t = params->t0();
  t.add_seconds(dso::FractionalSeconds(tsec));

  auto stokes = *(params->grav());

  /* compute acceleration for given epoch/position */
  [[maybe_unused]] Eigen::Matrix<double, 3, 3> g;
  Eigen::Matrix<double, 3, 1> acc = y.segment<3>(3);
  if (dso::sh2gradient_cunningham(stokes, y0.segment<3>(0), acc, g,
                                  stokes.max_degree(), stokes.max_order(), -1,
                                  -1, &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }
  printf("[deriv] At r=(%.3f, %.3f, %.3f) acc=(%.9f, %.9f, %.9f) [m/sec^2]\n", y0(0), y0(1), y0(2), acc(0), acc(1), acc(2));

  y.segment<3>(0) = y0.segment<3>(3);
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

  /* get states */
  auto start_t = sp3.start_epoch();

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

  ///* solid earth tide */
  // dso::SolidEarthTide se_tide;
  //

  /* setup integration parameters */
  dso::IntegrationParameters params;
  params.meops = &eop;
  params.mgrav = &stokes;

  /* setup the integrator */
  dso::Dop853 dop853(deriv, 6, &params, 1e-6, 1e-6);
  dop853.set_stiffness_check(10);

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
    bool position_ok = !block.flag.is_set(dso::Sp3Event::bad_abscent_position);
    if (position_ok && (!sp3err)) {
      gcrf2ecef(dso::MjdEpoch(block.t), eop, R, dRdt);
      if (!it) {
        state << block.state[0] * 1e3, block.state[1] * 1e3,
            block.state[2] * 1e3, block.state[4] * 1e-1, block.state[5] * 1e-1,
            block.state[6] * 1e-1;
        y.segment<3>(0) = R.transpose() * state.segment<3>(0);
        y.segment<3>(3) = R.transpose() * state.segment<3>(3) +
                          dRdt.transpose() * state.segment<3>(0);
        start_t = block.t;
        state = y;
        printf("Read first state off from sp3 file: %.12f %.6f %.6f %.6f %.9f "
               "%.9f %.9f\n",
               block.t.imjd().as_underlying_type() +
                   block.t.fractional_days().days(),
               block.state[0] * 1e3, block.state[1] * 1e3, block.state[2] * 1e3,
               block.state[4] * 1e-1, block.state[5] * 1e-1,
               block.state[6] * 1e-1);
        printf("Transformed to inertial           : %.12f %.6f %.6f %.6f %.9f "
               "%.9f %.9f\n",
               block.t.imjd().as_underlying_type() +
                   block.t.fractional_days().days(),
               state(0), state(1), state(2), state(3), state(4), state(5));
      } else {
        dso::FractionalSeconds sec =
            block.t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
                start_t);
        printf("Integrating %.3f sec from %.12f; IC=(%.6f %.6f %.6f %.9f %.9f "
               "%.9f)\n",
               sec.seconds(), dso::MjdEpoch(start_t).as_mjd(), state(0),
               state(1), state(2), state(3), state(4), state(5));
        if (dop853.integrate(dso::MjdEpoch(start_t), sec.seconds(), state, y)) {
          fprintf(stderr, "ERROR. Integration failed! sec is %.3f\n",
                  sec.seconds());
          return 1;
        }
        Eigen::VectorXd yt = y;
        y.segment<3>(0) = R * yt.segment<3>(0);
        y.segment<3>(3) = R * yt.segment<3>(3) + dRdt * yt.segment<3>(0);
        printf("%.12f %.6f %.6f %.6f %.9f %.9f %.9f %.6f %.6f %.6f %.9f %.9f "
               "%.9f\n",
               block.t.imjd().as_underlying_type() +
                   block.t.fractional_days().days(),
               block.state[0] * 1e3, block.state[1] * 1e3, block.state[2] * 1e3,
               block.state[4] * 1e-1, block.state[5] * 1e-1,
               block.state[6] * 1e-1, y(0), y(1), y(2), y(3), y(4), y(5));
      }
    }
    ++it;
    if (it >= 100)
      break;
  }

  return sp3err;
}
