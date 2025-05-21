#include "datetime/datetime_write.hpp"
#include "dop853.hpp"
#include "geodesy/units.hpp"
#include "iers/earth_rotation.hpp"
#include "iers/fundarg.hpp"
#include "iers/gravity.hpp"
#include "iers/iau.hpp"
#include "iers/icgemio.hpp"
#include "integration_parameters.hpp"
#include "sp3.hpp"
#include "yaml-cpp/yaml.h"
#include <cassert>
#include <cstdio>
#include <datetime/core/datetime_io_core.hpp>

constexpr const double GM_Moon = /*0.49028010560e13;*/ 4902.800076e9;
constexpr const double GM_Sun = /*1.32712442076e20;*/ 132712440040.944e9;

int prep_c2i(const dso::MjdEpoch &tai, dso::EopSeries &eops,
             Eigen::Quaterniond &q_c2tirs, Eigen::Quaterniond &q_tirs2i,
             double *fargs, dso::EopRecord &eopr) noexcept {

  const auto tt = tai.tai2tt();

  /* compute (X,Y) CIP and fundamental arguments (we are doing this here
   * to compute fargs).
   */
  double Xcip, Ycip;
  dso::xycip06a(tt, Xcip, Ycip, fargs);

  /* get (interpolate EOPs) */
  if (dso::EopSeries::out_of_bounds(eops.interpolate(tt, eopr))) {
    fprintf(stderr, "Failed to interpolate: Epoch is out of bounds!\n");
    return 1;
  }

  /* compute gmst using an approximate value for UT1 (linear interpolation) */
  double dut1_approx;
  eops.approx_dut1(tt, dut1_approx);
  [[maybe_unused]] const double gmst = dso::gmst(tt, tt.tt2ut1(dut1_approx));

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
  /* spherical crd for CIP */
  double d, e;
  dso::detail::xycip2spherical(Xcip, Ycip, d, e);
  const double era = dso::era00(tt.tt2ut1(eopr.dut()));

  q_c2tirs = dso::detail::c2tirs(era, s, d, e);
  q_tirs2i = dso::detail::tirs2i(dso::sec2rad(eopr.xp()),
                                 dso::sec2rad(eopr.yp()), dso::sp00(tt));

  return 0;
}

int deriv(double tsec, Eigen::Ref<const Eigen::VectorXd> y0,
          dso::IntegrationParameters *params,
          Eigen::Ref<Eigen::VectorXd> y) noexcept {

  /* current time in TAI */
  dso::MjdEpoch t = params->t0().add_seconds(dso::FractionalSeconds(tsec));

  /* Dealunay args (14) */
  double fargs[14];

  /* Celestial to Terrestrial Matrix */
  dso::EopRecord eopr;
  Eigen::Quaterniond q_c2tirs, q_tirs2i;
  prep_c2i(t, params->eops(), q_c2tirs, q_tirs2i, fargs, eopr);

  /* satellite position in ECEF */
  const Eigen::VectorXd r_ecef = q_tirs2i * (q_c2tirs * y0.segment<3>(0));

  /* Gravity field stokes coefficients */
  dso::StokesCoeffs stokes(params->earth_gravity());
  /* compute acceleration for given epoch/position (ECEF) */
  [[maybe_unused]] Eigen::Matrix<double, 3, 3> g;
  Eigen::Matrix<double, 3, 1> acc_grav;
  if (dso::sh2gradient_cunningham(stokes, r_ecef, acc_grav, g,
                                  stokes.max_degree(), stokes.max_order(), -1,
                                  -1, &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }

  /* Third body perturbation */
  Eigen::Matrix<double, 3, 1> acc_moon;
  Eigen::Matrix<double, 3, 1> acc_sun;
  /* get Sun position & velocity in ICRF */
  Eigen::Matrix<double, 6, 1> rtb_sun;
  if (dso::planet_state(dso::Planet::SUN, t.tai2tt(), rtb_sun)) {
    fprintf(stderr, "ERROR Failed to compute Sun position!\n");
    return 2;
  }
  acc_sun = dso::point_mass_acceleration(y0.segment<3>(0),
                                         rtb_sun.segment<3>(0), GM_Sun);

  /* get Moon position in ICRF */
  Eigen::Matrix<double, 3, 1> rtb_moon;
  if (dso::planet_pos(dso::Planet::MOON, t.tai2tt(), rtb_moon)) {
    fprintf(stderr, "ERROR Failed to compute Moon position!\n");
    return 2;
  }
  acc_moon = dso::point_mass_acceleration(y0.segment<3>(0), rtb_moon, GM_Moon);

  /* ECEF to ICRF note that y = (v, a) and y0 = (r, v) */
  y.segment<3>(0) = y0.segment<3>(3);
  y.segment<3>(3) = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * acc_grav);

  // printf("Requested ODE to %.12f, i.e. %.9f sec away,\n", t.as_mjd(), tsec);
  // printf("\t(r, v) = (%.3f, %.3f, %.3f, %.6f, %.6f, %.6f) ->\n\t(v, a) =
  // (%.3f, %.3f, %.3f, %.6f, %.6f, %.6f)\n", y0(0), y0(1), y0(2), y0(3), y0(4),
  // y0(5), y(0), y(1), y(2), y(3), y(4), y(5));
  //{
  //  /* DEBUG */
  //  Eigen::Matrix<double, 3, 1> gravn;
  //  if (dso::sh2gradient_cunningham(stokes, r_ecef, gravn, g,
  //                                  stokes.max_degree(), stokes.max_order(),
  //                                  -1, -1, &(params->tw()), &(params->tm())))
  //                                  {
  //    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
  //    return 1;
  //  }
  //  Eigen::Matrix<double, 3, 1> grav1;
  //  if (dso::sh2gradient_cunningham(stokes, r_ecef, grav1, g,
  //                                  1, 1, -1,
  //                                  -1, &(params->tw()), &(params->tm()))) {
  //    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
  //    return 1;
  //  }
  //  const auto an1 = gravn - grav1;
  //  printf("%.15f %.15f %.15f %.15f\n", t.tai2gps().as_mjd(), an1(0), an1(1),
  //  an1(2));
  //}

  return 0;
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s CONFIG SP3_file [SAT_ID]\n", argv[0]);
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

  /* get starting epoch in TAI */
  auto start_t = sp3.start_epoch();
  if (!std::strcmp(sp3.time_sys(), "GPS")) {
    start_t = start_t.gps2tai();
  }

  {
    char buf[64];
    printf("Note: starting epoch in Sp3 is %s\n",
           dso::to_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSSF>(
               start_t, buf));
  }

  auto t2 = start_t.add_seconds(dso::seconds(10 * 86400));
  dso::IntegrationParameters params = dso::IntegrationParameters::from_config(
      argv[1], dso::MjdEpoch(start_t), dso::MjdEpoch(t2));
  params.mtai0 = dso::MjdEpoch(start_t);

  /* setup the integrator */
  dso::Dop853 dop853(deriv, 6, &params, 1e-9, 1e-12);
  dop853.set_stiffness_check(10);

  /* dummy */
  double fargs[14];
  dso::EopRecord eopr;

  Eigen::VectorXd state = Eigen::Matrix<double, 6, 1>::Zero();
  Eigen::VectorXd y = Eigen::Matrix<double, 6, 1>::Zero();
  std::size_t it = 0;
  dso::Sp3DataBlock block;
  int sp3err = 0;
  while (!sp3err) {
    /* get next redord from sp3 */
    sp3err = sp3.get_next_data_block(sv, block);
    if (sp3err > 0) {
      printf("Something went wrong ....status = %3d\n", sp3err);
      return 1;
    }

    bool position_ok = !block.flag.is_set(dso::Sp3Event::bad_abscent_position);
    if (position_ok && (!sp3err)) {
      /* got new SP3 record, flag is ok */
      auto block_tai = block.t;
      if (!std::strcmp(sp3.time_sys(), "GPS")) {
        block_tai = block_tai.gps2tai();
      }

      /* GRCF/ITRF transformation quaternions (at t of new sp3 record) */
      Eigen::Quaterniond q_c2tirs, q_tirs2i;
      prep_c2i(dso::MjdEpoch(block_tai), params.eops(), q_c2tirs, q_tirs2i,
               fargs, eopr);

      if (!it) {
        /*
         * first state of satellite in file; transform to celestial and store
         * as state and start_t. This is where we start integrating from.
         */
        state << block.state[0] * 1e3, block.state[1] * 1e3,
            block.state[2] * 1e3, block.state[4] * 1e-1, block.state[5] * 1e-1,
            block.state[6] * 1e-1;

        /* transform state to GCRF (from ITRF) */
        Eigen::Vector3d omega;
        omega << 0e0, 0e0, dso::earth_rotation_rate(eopr.lod());
        y.segment<3>(0) =
            q_c2tirs.conjugate() * (q_tirs2i.conjugate() * state.segment<3>(0));
        y.segment<3>(3) =
            q_c2tirs.conjugate() *
            (q_tirs2i.conjugate() * state.segment<3>(3) +
             omega.cross(q_tirs2i.conjugate() * state.segment<3>(0)));

        start_t = block_tai;
        state = y;

      } else {
        /* new entry; seconds since intial epoch */
        dso::FractionalSeconds sec =
            block_tai.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
                start_t);

        /* integrate from initial conditions to this epoch */
        if (dop853.integrate(dso::MjdEpoch(start_t), sec.seconds(), state, y)) {
          fprintf(stderr, "ERROR. Integration failed! sec is %.3f\n",
                  sec.seconds());
          return 1;
        }

        /* transform integration results to ECEF */
        Eigen::VectorXd yt = y;
        Eigen::Vector3d omega;
        omega << 0e0, 0e0, dso::earth_rotation_rate(eopr.lod());
        yt.segment<3>(0) = q_tirs2i * (q_c2tirs * y.segment<3>(0));
        yt.segment<3>(3) = q_tirs2i * (q_c2tirs * y.segment<3>(3) -
                                       omega.cross(q_c2tirs * y.segment<3>(0)));

        printf("%.12f %.6f %.6f %.6f %.9f %.9f %.9f %.6f %.6f %.6f %.9f %.9f "
               "%.9f\n",
               sec.seconds(), block.state[0] * 1e3, block.state[1] * 1e3,
               block.state[2] * 1e3, block.state[4] * 1e-1,
               block.state[5] * 1e-1, block.state[6] * 1e-1, y(0), y(1), y(2),
               y(3), y(4), y(5));
      }
    } /* new sp3 record, flag ok */

    ++it;
    if (it >= 100)
      break;
  }

  return sp3err;
}
