#include "dop853.hpp"
#include "geodesy/units.hpp"
#include "iers/earth_rotation.hpp"
#include "iers/fundarg.hpp"
#include "iers/gravity.hpp"
#include "iers/iau.hpp"
#include "iers/icgemio.hpp"
#include "iers/relativity.hpp"
#include "integration_parameters.hpp"
#include "sp3.hpp"
#include "sysnsats/occultation.hpp"
#include "sysnsats/srp.hpp"
#include "yaml-cpp/yaml.h"
#include <cassert>
#include <cstdio>
#ifdef USE_BOOST
#include <boost/numeric/odeint.hpp>
#include <stdexcept>
#endif

constexpr const double GM_Moon = /*0.49028010560e13;*/ 4902.800076e9;
constexpr const double GM_Sun = /*1.32712442076e20;*/ 132712440040.944e9;
constexpr const double EVERY_SEC = 180e0;

#ifdef USE_BOOST
/* Define the State Type for ODEINT */
typedef std::array<double, 6> state_type;
typedef boost::numeric::odeint::runge_kutta_dopri5<state_type> dopri5_type;
#endif

/* Compute relevant quaternions for the ITRS/GCRS transformation
 */
int prep_c2i(const dso::MjdEpoch &tai, dso::EopSeries &eops,
             Eigen::Quaterniond &q_c2tirs, Eigen::Quaterniond &q_tirs2i,
             double *fargs, dso::EopRecord &eopr) noexcept {

  /* epoch of request in TT */
  const auto tt = tai.tai2tt();

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

#ifdef USE_BOOST
struct EomSystem {
  dso::IntegrationParameters *params;
  EomSystem(dso::IntegrationParameters *p) noexcept : params(p) {}

  void operator()(const state_type &x, state_type &dxdt, double tsec) {
    Eigen::Matrix<double, 6, 1> y0;
    for (int i = 0; i < 6; i++)
      y0(i) = x[i];
#else
int deriv(double tsec, Eigen::Ref<const Eigen::VectorXd> y0,
          dso::IntegrationParameters *params,
          Eigen::Ref<Eigen::VectorXd> yp) noexcept {
#endif
    /* epoch of request in TT */
    const auto tt =
        (params->t0().add_seconds(dso::FractionalSeconds(tsec))).tai2tt();

    /* GCRS/ITRS at t */
    double fargs[14];
    dso::EopRecord eopr;
    Eigen::Quaterniond q_c2tirs, q_tirs2i;
    prep_c2i(tt.tt2tai(), params->eops(), q_c2tirs, q_tirs2i, fargs, eopr);
    Eigen::Vector3d omega;
    omega << 0e0, 0e0, dso::earth_rotation_rate(eopr.lod());

    /* get Sun position & velocity in ICRF */
    Eigen::Matrix<double, 6, 1> rsun;
    if (dso::planet_state(dso::Planet::SUN, tt, rsun)) {
      fprintf(stderr, "ERROR Failed to compute Sun position!\n");
#ifdef USE_BOOST
      throw std::runtime_error(
          "ERROR. Failed computing derivative [Sun position]\n");
#else
    return 100;
#endif
    }

    /* get Moon position in ICRF */
    Eigen::Matrix<double, 3, 1> rmoon;
    if (dso::planet_pos(dso::Planet::MOON, tt, rmoon)) {
      fprintf(stderr, "ERROR Failed to compute Moon position!\n");
#ifdef USE_BOOST
      throw std::runtime_error(
          "ERROR. Failed computing derivative [Moon position]\n");
#else
    return 101;
#endif
    }

    /* state in ITRS (from GCRS) */
    Eigen::Matrix<double, 6, 1> itrs = Eigen::Matrix<double, 6, 1>::Zero();
    itrs.segment<3>(0) = q_tirs2i * (q_c2tirs * y0.segment<3>(0));
    itrs.segment<3>(3) = q_tirs2i * (q_c2tirs * y0.segment<3>(3) -
                                     omega.cross(q_c2tirs * y0.segment<3>(0)));

    /* accumulated acceleration and gradient in ITRS */
    Eigen::Vector3d ai = Eigen::Vector3d::Zero();
    [[maybe_unused]] Eigen::Matrix<double, 3, 3> gi;
    /* accumulated acceleration in GCRS */
    Eigen::Vector3d ac = Eigen::Vector3d::Zero();

    /* accumulated SH coeffs
    TODO!! WARNING!! What if some other SH coeffs (e.g. ocean tide) have (n,m)>
    gravity(n,m)? write a function as member of IntegrationParameters that
    return a StokesCoeffs of some degree and order
    */
    auto acstokes{params->earth_gravity()};

    /* add Solid Earth Tides to SH coeffs */
    if (params->solid_earth_tide()) {
      params->solid_earth_tide()->stokes_coeffs(
          tt, tt.tt2ut1(eopr.dut()), q_tirs2i * (q_c2tirs * rmoon),
          q_tirs2i * (q_c2tirs * rsun.segment<3>(0)), fargs);
      /* add SET effect */
      acstokes += params->solid_earth_tide()->stokes_coeffs();
    }

    /* add Ocean Tides to SH coeffs */
    if (params->ocean_tide()) {
      params->ocean_tide()->stokes_coeffs(tt, tt.tt2ut1(eopr.dut()), fargs);
      acstokes += params->ocean_tide()->stokes_coeffs();
    }

    /* add Pole Tide to SH coeffs */
    if (params->pole_tide()) {
      double dC21, dS21;
      params->pole_tide()->stokes_coeffs(tt, eopr.xp(), eopr.yp(), dC21, dS21);
      acstokes.C(2, 1) += dC21;
      acstokes.S(2, 1) += dS21;
    }

    /* add Ocean Pole Tide to SH coeffs */
    if (params->ocean_pole_tide()) {
      if (params->ocean_pole_tide()->stokes_coeffs(tt, eopr.xp(), eopr.yp())) {
        fprintf(stderr, "ERROR Failed computing Stokes Coefficients\n");
#ifdef USE_BOOST
        throw std::runtime_error(
            "ERROR. Failed computing derivative [Ocean Pole Tide]\n");
#else
      return 102;
#endif
      }
      acstokes += params->ocean_pole_tide()->stokes_coeffs();
    }

    /* add deAliasing to SH coeffs */
    if (params->dealias()) {
      /*
      TODO!! WARNING!! The dealias instance should have a function that appends
      the coefficients at a StokesCoeffs instance!
      */
      auto tempstokes{params->earth_gravity()};
      if (params->dealias()->coefficients_at(
              dso::from_mjdepoch<dso::nanoseconds>(tt), tempstokes)) {
        fprintf(stderr, "Failed interpolating dealiasing coefficients\n");
#ifdef USE_BOOST
        throw std::runtime_error(
            "ERROR. Failed computing derivative [deAliasing]\n");
#else
      return 103;
#endif
      }
      acstokes += tempstokes;
    }

    /* add atmospheric tides to SH coeffs */
    if (params->atmospheric_tide()) {
      params->atmospheric_tide()->stokes_coeffs(tt, tt.tt2ut1(eopr.dut()),
                                                fargs);
      acstokes += params->atmospheric_tide()->stokes_coeffs();
    }

    /* acceleration from accumulated SH expansion */
    if (dso::sh2gradient_cunningham(acstokes,
                                    // params->earth_gravity(),
                                    itrs.segment<3>(0), ai, gi,
                                    params->earth_gravity().max_degree(),
                                    params->earth_gravity().max_order(), -1, -1,
                                    &(params->tw()), &(params->tm()))) {
      fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
#ifdef USE_BOOST
      throw std::runtime_error(
          "ERROR. Failed computing derivative [sh2gradient]\n");
#else
    return 104;
#endif
    }

    /* Third Body perturbations and Relativity (IERS 2010) */
    {
      ac += dso::point_mass_acceleration(y0.segment<3>(0), rsun.segment<3>(0),
                                         GM_Sun);
      ac += dso::point_mass_acceleration(y0.segment<3>(0), rmoon, GM_Moon);
      /* Relativistic Correction */
      ac += dso::iers2010_relativistic_acceleration(y0, rsun);
    }

    if (params->mCr) {
      const double of =
          /* dso::conical_occultation(y0.segment<3>(0), rsun.segment<3>(0)); */
          dso::conical_occultation(y0.segment<3>(0), rsun.segment<3>(0), Eigen::Matrix<double,3,1>::Zero(), ::iers2010::Re) * 
          dso::conical_occultation(y0.segment<3>(0), rsun.segment<3>(0), rmoon, 1737e3);
      
      if (of > 0e0) {
        /* Solar Radiation Pressure */
        if (params->matt) {
          /* get attitude */
          if (params->matt->attitude_at(tt, *(params->mattdata))) {
            fprintf(stderr, "[ERROR] Failed getting attitude!\n");
#ifdef USE_BOOST
            throw std::runtime_error(
                "ERROR. Failed computing derivative [attitude]\n");
#else
          return 200;
#endif
          }
          /* we may need (depending on satellite) the satellite-to-sun vector */
          Eigen::Vector3d sat2sun =  rsun.segment<3>(0) - y0.segment<3>(0);
          /* compute SRP acceleration */
          const Eigen::Vector3d tmp = 
          /*ac +=*/ (params->mCr * of) *
                dso::solar_radiation_pressure(
                    params->msatmm->rotate_macromodel(
                        params->mattdata->quaternions(),
                        params->mattdata->angles(), &sat2sun),
                    y0.segment<3>(0), rsun.segment<3>(0), params->msatmm->satellite_mass());
          // printf("%.12f %.15f %.15f %.15f\n", tt.as_mjd(), tmp(0), tmp(1), tmp(2));
          ac += tmp;
        } else {
          ac += (params->mCr * of) *
                dso::solar_radiation_pressure(
                    params->msatmm->srp_cannoball_area(), y0.segment<3>(0),
                    rsun.segment<3>(0), params->msatmm->satellite_mass());
        }
      }
    }

#ifdef USE_BOOST
    Eigen::Matrix<double, 6, 1> yp = Eigen::Matrix<double, 6, 1>::Zero();
#endif

    /* form the derivative vector */
    yp.segment<3>(0) = y0.segment<3>(3);
    yp.segment<3>(3) = ac + q_c2tirs.conjugate() * (q_tirs2i.conjugate() * ai);

#ifdef USE_BOOST
    /* form the derivative vector */
    dxdt[0] = yp(3);
    dxdt[1] = yp(4);
    dxdt[2] = yp(5);
    dxdt[3] = yp(0);
    dxdt[4] = yp(1);
    dxdt[5] = yp(2);

    return;
  }
}; /* struct EomSystem */
#else
  return 0;
}
#endif

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
#ifdef USE_BOOST
  auto stepper = make_controlled(1e-8, 1e-8, dopri5_type());
#else
  dso::Dop853 dop853(deriv, 6, &params, 1e-9, 1e-12);
  dop853.set_stiffness_check(10);
#endif

  Eigen::VectorXd y0 = Eigen::Matrix<double, 6, 1>::Zero();
  Eigen::VectorXd y = Eigen::Matrix<double, 6, 1>::Zero();
  Eigen::Quaterniond q_c2tirs, q_tirs2i;
  Eigen::Vector3d omega;
  dso::MjdEpoch t0;
  dso::MjdEpoch tp;
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
      tp = t0;
      params.t0() = tai;
    } else if (tai.diff<dso::DateTimeDifferenceType::FractionalSeconds>(tp)
                   .seconds() >= EVERY_SEC) {
      /* seconds since initial epoch */
      dso::FractionalSeconds sec =
          tai.diff<dso::DateTimeDifferenceType::FractionalSeconds>(t0);

      /* integrate from first SP3 record to here - setup initial conditions */
      params.t0() = t0;
      Eigen::VectorXd yc = Eigen::Matrix<double, 6, 1>::Zero();

      /* note: we are starting from the top here, so reload the attitude
       * stream
       */
      if (params.matt) {
        if (params.matt->reload()) {
          fprintf(stderr, "ERROR. Failed reloading attitude stream!\n");
          return -100;
        }
      }

#ifdef USE_BOOST
      state_type x = {y0(0), y0(1), y0(2), y0(3), y0(4), y0(5)};
      boost::numeric::odeint::integrate_adaptive(stepper, EomSystem(&params), x,
                                                 0e0, sec.seconds(), 1e-3);
      yc << x[0], x[1], x[2], x[3], x[4], x[5];
#else
      if (dop853.integrate(0e0, sec.seconds(), y0, yc)) {
        fprintf(stderr, "ERROR. Integration failed! sec is %.3f\n",
                sec.seconds());
        return 1;
      }
#endif

      /* compare in ITRS */
      Eigen::VectorXd yi = Eigen::Matrix<double, 6, 1>::Zero();
      yi.segment<3>(0) = q_tirs2i * (q_c2tirs * yc.segment<3>(0));
      yi.segment<3>(3) = q_tirs2i * (q_c2tirs * yc.segment<3>(3) -
                                     omega.cross(q_c2tirs * yc.segment<3>(0)));

      printf("%.9f %.3f %.3f %.3f %.6f %.6f %.6f %.3f %.3f %.3f %.6f %.6f "
             "%.6f\n",
             sec.seconds(), y(0), y(1), y(2), y(3), y(4), y(5), yi(0), yi(1),
             yi(2), yi(3), yi(4), yi(5));

      tp = tai;
    }

    ++it;
    if (it > 600)
      break;
  }

  printf("# Num of epochs parsed/used: %ld\n", it);
}
