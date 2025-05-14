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
#include <cstdio>
#include <iostream> // DEBUG

constexpr const double GM_Moon = /*0.49028010560e13;*/ 4902.800076e9;
constexpr const double GM_Sun = /*1.32712442076e20;*/ 132712440040.944e9;

/* Compute and return the rotation matrix and its derivative between GCRF and
 * ECEF frames.
 */
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

  /* current time in GPST (debugging) */
  const auto tgps = t.tai2gps();

  /* Dealunay args (14) .. to be computed */
  double fargs[14];

  /* Celestial to Terrestrial Matrix */
  dso::EopRecord eopr;
  Eigen::Quaterniond q_c2tirs, q_tirs2i;
  prep_c2i(dso::MjdEpoch(t), params->eops(), q_c2tirs, q_tirs2i, fargs, eopr);

  /* satellite position in ECEF */
  const Eigen::VectorXd r_ecef = q_tirs2i * (q_c2tirs * y0.segment<3>(0));

  /* Gravity field stokes coefficients */
  dso::StokesCoeffs stokes(params->earth_gravity());
  stokes.C(0, 0) = 0e0;
  stokes.C(1, 0) = stokes.C(1, 1) = 0e0;
  stokes.S(1, 1) = 0e0;

  /* compute acceleration for given epoch/position (ECEF) */
  [[maybe_unused]] Eigen::Matrix<double, 3, 3> g;
  Eigen::Matrix<double, 3, 1> acc_grav;
  if (dso::sh2gradient_cunningham(stokes, r_ecef, acc_grav, g,
                                  stokes.max_degree(), stokes.max_order(), -1,
                                  -1, &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }
  {
    /* print acceleration in ECI (P1)*/
    const auto acc = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * acc_grav);
    printf("%.15f %.15f %.15f %.15f", tgps.as_mjd(), acc(0), acc(1), acc(2));
    y.segment<3>(0) += acc;
  }

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
  /* print acceleration in ECI (P2)*/
  printf(" %.15f %.15f %.15f", acc_sun(0), acc_sun(1), acc_sun(2));
  y.segment<3>(0) += acc_sun;

  /* get Moon position in ICRF */
  Eigen::Matrix<double, 3, 1> rtb_moon;
  if (dso::planet_pos(dso::Planet::MOON, t.tai2tt(), rtb_moon)) {
    fprintf(stderr, "ERROR Failed to compute Moon position!\n");
    return 2;
  }
  acc_moon = dso::point_mass_acceleration(y0.segment<3>(0), rtb_moon, GM_Moon);
  /* print acceleration in ECI (P3)*/
  printf(" %.15f %.15f %.15f", acc_moon(0), acc_moon(1), acc_moon(2));
  y.segment<3>(0) += acc_moon;

  /* Relativistic Correction */
  Eigen::Matrix<double, 3, 1> acc_rel =
      dso::iers2010_relativistic_acceleration(y0, rtb_sun);
  /* print acceleration in ECI (P4)*/
  printf(" %.15f %.15f %.15f", acc_rel(0), acc_rel(1), acc_rel(2));
  y.segment<3>(0) += acc_rel;

  /* Solid Earth Tide (ITRF) */
  Eigen::Matrix<double, 3, 1> acc_set;
  params->solid_earth_tide()->stokes_coeffs(
      t.tai2tt(), t.tai2ut1(eopr.dut()), q_tirs2i * (q_c2tirs * rtb_moon),
      q_tirs2i * (q_c2tirs * rtb_sun.segment<3>(0)), fargs);
  if (dso::sh2gradient_cunningham(params->solid_earth_tide()->stokes_coeffs(),
                                  r_ecef, acc_set, g, -1, -1, -1, -1,
                                  &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }
  {
    /* print acceleration in ECI (P5)*/
    // const auto acc = R.transpose() * acc_set;
    const auto acc = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * acc_set);
    printf(" %.15f %.15f %.15f", acc(0), acc(1), acc(2));
    y.segment<3>(0) += acc;
  }

  /* Solid Earth Pole Tide */
  dso::StokesCoeffs sept_stokes(3, 3);
  sept_stokes.clear();
  Eigen::Matrix<double, 3, 1> acc_setp;
  if (dso::PoleTide::stokes_coeffs(t.tai2tt(), eopr.xp(), eopr.yp(),
                                   sept_stokes.C(2, 1), sept_stokes.S(2, 1))) {
    return 1;
  }
  if (dso::sh2gradient_cunningham(sept_stokes, r_ecef, acc_setp, g, -1, -1, -1,
                                  -1, &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }
  {
    /* print acceleration in ECI (P6)*/
    // const auto acc = R.transpose() * acc_setp;
    const auto acc = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * acc_setp);
    printf(" %.15f %.15f %.15f", acc(0), acc(1), acc(2));
    y.segment<3>(0) += acc;
  }

  /* Ocean Pole Tide */
  Eigen::Matrix<double, 3, 1> acc_otp;
  if (params->ocean_pole_tide()->stokes_coeffs(
          t.tai2tt(), eopr.xp(), eopr.yp(),
          params->ocean_pole_tide()->max_degree(),
          params->ocean_pole_tide()->max_order())) {
    fprintf(stderr, "ERROR Failed computing Stokes Coefficients\n");
    return 1;
  }
  params->ocean_pole_tide()->stokes_coeffs().C(0, 0) =
      params->ocean_pole_tide()->stokes_coeffs().C(1, 0) =
          params->ocean_pole_tide()->stokes_coeffs().C(1, 1) = 0e0;
  params->ocean_pole_tide()->stokes_coeffs().S(1, 1) = 0e0;
  if (dso::sh2gradient_cunningham(params->ocean_pole_tide()->stokes_coeffs(),
                                  r_ecef, acc_otp, g,
                                  params->ocean_pole_tide()->max_degree(),
                                  params->ocean_pole_tide()->max_order(), -1,
                                  -1, &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }
  {
    // const auto acc = R.transpose() * acc_otp;
    const auto acc = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * acc_otp);
    /* print acceleration in ECI (P7)*/
    printf(" %.15f %.15f %.15f", acc(0), acc(1), acc(2));
    y.segment<3>(0) += acc;
  }

  /* Dealiasing
   * WARNING !!
   * Here, we are storing the Stokes coefficients in stokes (created earlier).
   */
  Eigen::Matrix<double, 3, 1> acc_das;
  if (params->dealias()->coefficients_at(
          dso::from_mjdepoch<dso::nanoseconds>(t.tai2tt()), stokes)) {
    fprintf(stderr, "Failed interpolating dealiasing coefficients\n");
    return 1;
  }
  stokes.C(0, 0) = stokes.C(1, 0) = stokes.C(1, 1) = 0e0;
  stokes.S(1, 1) = 0e0;
  if (dso::sh2gradient_cunningham(stokes, r_ecef, acc_das, g,
                                  stokes.max_degree(), stokes.max_order(), -1,
                                  -1, &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }
  {
    /* print acceleration in ECI (P8)*/
    // const auto acc = R.transpose() * acc_das;
    const auto acc = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * acc_das);
    printf(" %.15f %.15f %.15f", acc(0), acc(1), acc(2));
    y.segment<3>(0) += acc;
  }

  /* Atmospheric Tide */
  Eigen::Matrix<double, 3, 1> acc_atm;
  /* compute Stokes coeffs (for atm. tides) */
  params->atmospheric_tide()->stokes_coeffs(t.tai2tt(), t.tai2ut1(eopr.dut()),
                                            fargs);
  /* for the test, degree one coefficients are not taken into account */
  params->atmospheric_tide()->stokes_coeffs().C(0, 0) =
      params->atmospheric_tide()->stokes_coeffs().C(1, 0) =
          params->atmospheric_tide()->stokes_coeffs().C(1, 1) = 0e0;
  params->atmospheric_tide()->stokes_coeffs().S(1, 1) = 0e0;
  /* compute acceleration for given epoch/position (ITRF) */
  if (dso::sh2gradient_cunningham(params->atmospheric_tide()->stokes_coeffs(),
                                  r_ecef, acc_atm, g,
                                  params->atmospheric_tide()->max_degree(),
                                  params->atmospheric_tide()->max_order(), -1,
                                  -1, &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }
  {
    /* print acceleration in ECI (P9)*/
    // const auto acc = R.transpose() * acc_atm;
    const auto acc = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * acc_atm);
    printf(" %.15f %.15f %.15f", acc(0), acc(1), acc(2));
    y.segment<3>(0) += acc;
  }

  /* Ocean Tide */
  Eigen::Matrix<double, 3, 1> acc_ot;
  /* compute Stokes coeffs (for ocean tides) */
  params->ocean_tide()->stokes_coeffs(t.tai2tt(), t.tai2ut1(eopr.dut()), fargs);
  /* for the test, degree one coefficients are not taken into account */
  params->ocean_tide()->stokes_coeffs().C(0, 0) =
      params->ocean_tide()->stokes_coeffs().C(1, 0) =
          params->ocean_tide()->stokes_coeffs().C(1, 1) = 0e0;
  params->ocean_tide()->stokes_coeffs().S(1, 1) = 0e0;
  /* compute acceleration for given epoch/position (ITRF) */
  if (dso::sh2gradient_cunningham(params->ocean_tide()->stokes_coeffs(), r_ecef,
                                  acc_ot, g, params->ocean_tide()->max_degree(),
                                  params->ocean_tide()->max_order(), -1, -1,
                                  &(params->tw()), &(params->tm()))) {
    fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
    return 1;
  }
  {
    /* print acceleration in ECI (P10)*/
    // const auto acc = R.transpose() * acc_ot;
    const auto acc = q_c2tirs.conjugate() * (q_tirs2i.conjugate() * acc_ot);
    printf(" %.15f %.15f %.15f", acc(0), acc(1), acc(2));
    y.segment<3>(0) += acc;
  }

  /* set velocity vector (ICRF) */
  y.segment<3>(3) = y.segment<3>(0);
  y.segment<3>(0) = y0.segment<3>(3);

  /* ECEF to ICRF note that y = (v, a) and y0 = (r, v) */
  // y.segment<3>(3) =
  //     R.transpose() * acc + acc_tb;
  // y.segment<3>(3) = (acc_moon + acc_sun + acc_rel);
  // y.segment<3>(3) += R.transpose() * acc_set;
  // y.segment<3>(3) += R.transpose() * acc;

  printf("\n");
  return 0;
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s [CONFIG] [00orbit_itrf.sp3] \n", argv[0]);
    return 1;
  }

  /* parse the yaml configuration file */
  // YAML::Node config = YAML::LoadFile(argv[1]);

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

  auto t2 = start_t.add_seconds(dso::seconds(5 * 86400));
  dso::IntegrationParameters params = dso::IntegrationParameters::from_config(
      argv[1], dso::MjdEpoch(start_t), dso::MjdEpoch(t2));

  ///* EOPs */
  // dso::EopSeries eop;
  //{
  //   std::string tmp = config["eop"].as<std::string>();
  //   auto t1 = start_t;
  //   t1.add_seconds(dso::seconds(-86400));
  //   auto t2 = start_t;
  //   t2.add_seconds(dso::seconds(5 * 86400));
  //   if (dso::parse_iers_C04(tmp.c_str(), dso::MjdEpoch(t1),
  //   dso::MjdEpoch(t2),
  //                           eop)) {
  //     fprintf(stderr, "ERROR Failed parsing eop file\n");
  //     return 1;
  //   }
  //   eop.regularize();
  // }

  ///* load planetary ephemeris kernels */
  // std::string tmp = config["planetary-ephemeris"]["bsp"].as<std::string>();
  // dso::load_spice_kernel(tmp.c_str());
  // tmp = config["planetary-ephemeris"]["tls"].as<std::string>();
  // dso::load_spice_kernel(tmp.c_str());

  ///* Earth's gravity field */
  // dso::StokesCoeffs stokes;
  //{
  //   tmp = config["gravity"]["model"].as<std::string>();
  //   int DEGREE = config["gravity"]["degree"].as<int>();
  //   int ORDER = config["gravity"]["order"].as<int>();
  //   dso::Icgem icgem(tmp.c_str());
  //   if (icgem.parse_data(DEGREE, ORDER, start_t, stokes)) {
  //     fprintf(stderr, "ERROR Failed reading gravity model!\n");
  //     return 1;
  //   }
  //   /* checks */
  //   assert(stokes.max_degree() == DEGREE);
  //   assert(stokes.max_order() == ORDER);
  // }

  ///* Solid Earth Tides */
  // dso::SolidEarthTide setide(iers2010::GMe, iers2010::Re, GM_Sun, GM_Moon);

  ///* Solid Earth Pole Tide */
  //// no need to initialize, will use the static function:
  //// PoleTide::stokes_coeffs

  ///* Ocean Pole Tide */
  // tmp = config["ocean-pole-tide"]["model"].as<std::string>();
  // dso::OceanPoleTide *opt = nullptr;
  // if (tmp == "Desai02") {
  //   int DEGREE = config["ocean-pole-tide"]["degree"].as<int>();
  //   int ORDER = config["ocean-pole-tide"]["order"].as<int>();
  //   tmp = config["ocean-pole-tide"]["coeffs"].as<std::string>();
  //   opt = new dso::OceanPoleTide(DEGREE, ORDER, tmp.c_str());
  // }

  ///* Dealiasing */
  // tmp = config["dealiasing"]["model"].as<std::string>();
  // dso::Aod1bDataStream<dso::AOD1BCoefficientType::GLO> *dap = nullptr;
  // if (tmp != "") {
  //   const auto f = config["dealiasing"]["data-file"].as<std::string>();
  //   const auto d = config["dealiasing"]["data-dir"].as<std::string>();
  //   dap = new dso::Aod1bDataStream<dso::AOD1BCoefficientType::GLO>(f.c_str(),
  //                                                                  d.c_str());
  //   dap->initialize();
  // }

  ///* Atmospheric Tide */
  // tmp = config["atmospheric-tide"]["model"].as<std::string>();
  // dso::AtmosphericTide atm;
  // if (tmp == "AOD1B RL06") {
  //   int DEGREE = config["atmospheric-tide"]["degree"].as<int>();
  //   int ORDER = config["atmospheric-tide"]["order"].as<int>();
  //   std::string data_dir =
  //       config["atmospheric-tide"]["data_dir"].as<std::string>();
  //   YAML::Node tide_atlas =
  //   config["atmospheric-tide"]["tide_atlas_from_aod1b"]; for
  //   (YAML::const_iterator it = tide_atlas.begin(); it != tide_atlas.end();
  //        ++it) {
  //     std::string filename = it->second.as<std::string>();
  //     std::string full_path = data_dir + "/" + filename;
  //     atm.append_wave(full_path.c_str(), DEGREE, ORDER);
  //   }
  // }

  ///* Ocean Tide */
  // tmp = config["ocean-tide"]["model"].as<std::string>();
  // dso::OceanTide *ot = nullptr;
  //{
  //   int DEGREE = config["ocean-tide"]["degree"].as<int>();
  //   int ORDER = config["ocean-tide"]["order"].as<int>();
  //   std::string data_dir =
  //   config["ocean-tide"]["data_dir"].as<std::string>(); std::string
  //   groops_file_list =
  //       config["ocean-tide"]["groops_file_list"].as<std::string>();
  //   ot = new dso::OceanTide(dso::groops_ocean_atlas(
  //       groops_file_list.c_str(), data_dir.c_str(), DEGREE, ORDER));
  // }

  ///* setup integration parameters */
  // dso::IntegrationParameters params;
  // params.meops = &eop;
  // params.mgrav = &stokes;
  // params.mse_tide = &setide;
  params.mtai0 = dso::MjdEpoch(start_t);
  // params.mop_tide = opt;
  // params.mdealias = dap;
  // if (params.mdealias) {
  //   params.mdealias_maxdegree = config["dealiasing"]["degree"].as<int>();
  //   params.mdealias_maxorder = config["dealiasing"]["order"].as<int>();
  // }
  // params.mat_tide = &atm;
  // if (params.mat_tide) {
  //   params.matm_maxdegree = config["atmospheric-tide"]["degree"].as<int>();
  //   params.matm_maxorder = config["atmospheric-tide"]["order"].as<int>();
  // }
  // params.moc_tide = ot;
  // if (params.moc_tide) {
  //   params.moc_maxdegree = config["ocean-tide"]["degree"].as<int>();
  //   params.moc_maxorder = config["ocean-tide"]["order"].as<int>();
  // }

  /* setup the integrator */
  // dso::Dop853 dop853(deriv, 6, &params, 1e-9, 1e-12);
  // dop853.set_stiffness_check(10);

  /* dummy */
  double fargs[14];
  dso::EopRecord eopr;

  /* Just for testing Vs costg */
  Eigen::VectorXd state = Eigen::Matrix<double, 6, 1>::Zero();
  Eigen::VectorXd y = Eigen::Matrix<double, 6, 1>::Zero();
  // Eigen::Matrix<double, 3, 3> R, dRdt;
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

    /* time of current block in TAI */
    auto block_tai = block.t;
    if (!std::strcmp(sp3.time_sys(), "GPS")) {
      block_tai = block_tai.gps2tai();
    }

    /* GCRF to ITRF rotation matrix */
    Eigen::Quaterniond q_c2tirs, q_tirs2i;
    prep_c2i(dso::MjdEpoch(block_tai), params.eops(), q_c2tirs, q_tirs2i,
             fargs, eopr);

    /* get state for current epoch ITRF */
    state << block.state[0] * 1e3, block.state[1] * 1e3, block.state[2] * 1e3,
        block.state[4] * 1e-1, block.state[5] * 1e-1, block.state[6] * 1e-1;

    /* transform state to GCRF (from ITRF) */
    Eigen::Vector3d omega;
    omega << 0e0, 0e0, dso::earth_rotation_rate(eopr.lod());
    y.segment<3>(0) =
        q_c2tirs.conjugate() * (q_tirs2i.conjugate() * state.segment<3>(0));
    y.segment<3>(3) = q_c2tirs.conjugate() *
                      (q_tirs2i.conjugate() * state.segment<3>(3) +
                       omega.cross(q_tirs2i.conjugate() * state.segment<3>(0)));
    state = y;

    /* seconds since initial epoch */
    dso::FractionalSeconds sec =
        block_tai.diff<dso::DateTimeDifferenceType::FractionalSeconds>(start_t);

    /* compute accelerations at this epoch */
    if (deriv(sec.seconds(), state, &params, y)) {
      fprintf(stderr, "ERROR. Failed computing derivative!\n");
      return 1;
    }

    ++it;
  }

  printf("# Num of epochs parsed/used: %ld\n", it);

  return sp3err;
}
