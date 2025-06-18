#include "iers/icgemio.hpp"
#include "integration_parameters.hpp"
#include "yaml-cpp/yaml.h"

bool is_nonempty_string(const YAML::Node &node) {
  return node.IsDefined() && node.IsScalar() && !node.as<std::string>().empty();
}

dso::IntegrationParameters
dso::IntegrationParameters::from_config(const char *fn,
                                        const dso::MjdEpoch &tt_start,
                                        const dso::MjdEpoch &tt_stop) {
  /* parse the yaml configuration file */
  YAML::Node config = YAML::LoadFile(fn);

  IntegrationParameters params;

  /* check dates */
  if (tt_stop < tt_start) {
    std::string err_msg =
        "[ERROR] Invalid start/stop dates in config file (traceback:" +
        std::string(__func__) + ")\n";
    throw std::runtime_error(err_msg);
  }

  /* find mid-epoch (TT) */
  auto dsec =
      tt_stop.diff<dso::DateTimeDifferenceType::FractionalSeconds>(tt_start)
          .seconds();
  const auto tt_midepoch =
      tt_start.add_seconds(dso::FractionalSeconds(dsec / 2.));

  /* fill EOP Series of the IntegrationParameters instance */
  {
    std::string tmp = config["eop"].as<std::string>();
    const auto t1 = tt_start.add_seconds(dso::FractionalSeconds(-86400));
    const auto t2 = tt_stop.add_seconds(dso::FractionalSeconds(86400));
    if (dso::parse_iers_C04(tmp.c_str(), dso::MjdEpoch(t1), dso::MjdEpoch(t2),
                            params.meops)) {
      std::string err_msg = "[ERROR] Failed parsing EOP file " + tmp +
                            " (traceback:" + std::string(__func__) + ")\n";
      throw std::runtime_error(err_msg);
    }
    params.meops.regularize();
  }

  /* load planetary ephemeris kernels */
  std::string tmp = config["planetary-ephemeris"]["bsp"].as<std::string>();
  dso::load_spice_kernel(tmp.c_str());
  tmp = config["planetary-ephemeris"]["tls"].as<std::string>();
  dso::load_spice_kernel(tmp.c_str());

  /* read Earth's gravity field into IntegrationParameters instance */
  {
    tmp = config["gravity"]["model"].as<std::string>();
    int DEGREE = config["gravity"]["degree"].as<int>();
    int ORDER = config["gravity"]["order"].as<int>();
    dso::Icgem icgem(tmp.c_str());
    if (icgem.parse_data(DEGREE, ORDER,
                         dso::from_mjdepoch<dso::nanoseconds>(tt_midepoch),
                         params.mgrav)) {
      std::string err_msg = "[ERROR] Failed parsing gravity field model " +
                            tmp + " (traceback:" + std::string(__func__) +
                            ")\n";
      throw std::runtime_error(err_msg);
    }
    /* checks */
    assert(params.mgrav.max_degree() == DEGREE);
    assert(params.mgrav.max_order() == ORDER);
  }

  /* Solid Earth Tides (all parameters are default initialized) */
  if (is_nonempty_string(config["solid-earth-tide"]["model"])) {
    tmp = config["solid-earth-tide"]["model"].as<std::string>();
    if (tmp == "IERS2010") {
      params.mse_tide = new dso::SolidEarthTide();
    } else {
      std::string err_msg = "[ERROR] Unknown Solid Earth Tide model given " +
                            tmp + " (traceback:" + std::string(__func__) +
                            ")\n";
      throw std::runtime_error(err_msg);
    }
  }

  /* Ocean Tide (if model name is non-empty) */
  dso::OceanTide *ot = nullptr;
  if (is_nonempty_string(config["ocean-tide"]["model"])) {
    /* basic parameters, must exist */
    tmp = config["ocean-tide"]["model"].as<std::string>();
    int DEGREE = config["ocean-tide"]["degree"].as<int>();
    int ORDER = config["ocean-tide"]["order"].as<int>();
    std::string data_dir = config["ocean-tide"]["data_dir"].as<std::string>();
    std::string groops_file_list =
        config["ocean-tide"]["groops_file_list"].as<std::string>();
    /* OceanTide model with admittance */
    if (is_nonempty_string(config["ocean-tide"]["groops_doodson02_file"]) &&
        is_nonempty_string(config["ocean-tide"]["groops_admittance03_file"])) {
      std::string file02 =
          config["ocean-tide"]["groops_doodson02_file"].as<std::string>();
      std::string file03 =
          config["ocean-tide"]["groops_admittance03_file"].as<std::string>();

      try {
        ot = new dso::OceanTide(dso::groops_ocean_atlas(
            groops_file_list.c_str(), file02.c_str(), file03.c_str(),
            data_dir.c_str(), DEGREE, ORDER));
      } catch (std::exception &e) {
        fprintf(stderr, e.what());
        std::string err_msg = "[ERROR] Failed to construct Ocean Tide model "
                              "with admittance from "
                              "config parameters " +
                              data_dir + ", " + groops_file_list + ", " +
                              file02 + ", " + file03 +
                              " (traceback:" + std::string(__func__) + ")\n";
        throw std::runtime_error(err_msg);
      }
    } else {
      /* OceanTide model without admittance */
      try {
        ot = new dso::OceanTide(dso::groops_ocean_atlas(
            groops_file_list.c_str(), data_dir.c_str(), DEGREE, ORDER));
      } catch (std::exception &e) {
        fprintf(stderr, e.what());
        std::string err_msg =
            "[ERROR] Failed to construct Ocean Tide model from "
            "config parameters " +
            data_dir + ", " + groops_file_list +
            " (traceback:" + std::string(__func__) + ")\n";
        throw std::runtime_error(err_msg);
      }
    }
  }
  params.moc_tide = ot;

  /* Solid Earth Pole Tide */
  if (is_nonempty_string(config["solid-earth-pole-tide"]["model"])) {
    tmp = config["solid-earth-pole-tide"]["model"].as<std::string>();
    if (tmp == "IERS2010") {
      params.mep_tide = new dso::PoleTide();
    } else {
      std::string err_msg =
          "[ERROR] Unknown (Solid Earth) Pole Tide model given " + tmp +
          " (traceback:" + std::string(__func__) + ")\n";
      throw std::runtime_error(err_msg);
    }
  }

  /* Ocean Pole Tide (if model name is non-empty) */
  dso::OceanPoleTide *opt = nullptr;
  {
    if (is_nonempty_string(config["ocean-pole-tide"]["model"])) {
      tmp = config["ocean-pole-tide"]["model"].as<std::string>();
      if (tmp == "Desai02") {
        /* Desai 2002 model, need degree, order and coeffs */
        int DEGREE = config["ocean-pole-tide"]["degree"].as<int>();
        int ORDER = config["ocean-pole-tide"]["order"].as<int>();
        std::string tfn = config["ocean-pole-tide"]["coeffs"].as<std::string>();
        opt = new dso::OceanPoleTide(DEGREE, ORDER, tfn.c_str());
      } else {
        std::string err_msg = "[ERROR] Unknown Ocean Pole Tide model given " +
                              tmp + " (traceback:" + std::string(__func__) +
                              ")\n";
        throw std::runtime_error(err_msg);
      }
    }
  }
  params.mop_tide = opt;

  /* Atmospheric Tide (if model name is non-empty) */
  dso::AtmosphericTide *atm = nullptr;
  if (is_nonempty_string(config["atmospheric-tide"]["model"])) {
    int DEGREE = config["atmospheric-tide"]["degree"].as<int>();
    int ORDER = config["atmospheric-tide"]["order"].as<int>();
    std::string data_dir =
        config["atmospheric-tide"]["data_dir"].as<std::string>();
    if (is_nonempty_string(config["atmospheric-tide"]["groops_file_list"])) {
      /* generating an AtmosphericTide instance from GROOPS files */
      std::string fn1 =
          config["atmospheric-tide"]["groops_file_list"].as<std::string>();
      try {
        atm = new dso::AtmosphericTide(dso::groops_atmospheric_atlas(
            fn1.c_str(), data_dir.c_str(), DEGREE, ORDER));
      } catch (std::exception &e) {
        fprintf(stderr, e.what());
        std::string err_msg =
            "[ERROR] Failed to construct Atmospheric Tide model from "
            "config parameters " +
            data_dir + ", " + fn1 + " (traceback:" + std::string(__func__) +
            ")\n";
        throw std::runtime_error(err_msg);
      }
    } else {
      /* generating an AtmosphericTide instance from AOD1B files */
      atm = new dso::AtmosphericTide();
      YAML::Node tide_atlas =
          config["atmospheric-tide"]["tide_atlas_from_aod1b"];
      for (YAML::const_iterator it = tide_atlas.begin(); it != tide_atlas.end();
           ++it) {
        std::string filename = it->second.as<std::string>();
        std::string full_path = data_dir + "/" + filename;
        if (atm->append_wave(full_path.c_str(), DEGREE, ORDER)) {
          std::string err_msg =
              "[ERROR] Failed appending tide (file) " + full_path +
              " to Atmospheric Tide model (traceback:" + std::string(__func__) +
              ")\n";
          throw std::runtime_error(err_msg);
        }
      }
    }
  }
  params.mat_tide = atm;

  /* Dealiasing (if model name is non-empty) */
  dso::Aod1bDataStream<dso::AOD1BCoefficientType::GLO> *dap = nullptr;
  if (is_nonempty_string(config["dealiasing"]["model"])) {
    const auto f = config["dealiasing"]["data-file"].as<std::string>();
    const auto d = config["dealiasing"]["data-dir"].as<std::string>();
    try {
      dap = new dso::Aod1bDataStream<dso::AOD1BCoefficientType::GLO>(f.c_str(),
                                                                     d.c_str());
      dap->initialize();
    } catch (std::exception &e) {
      fprintf(stderr, e.what());
      std::string err_msg =
          "[ERROR] Failed to construct an Aod1bDataStream<GLO> stream from "
          "config parameters " +
          f + ", " + d + " (traceback:" + std::string(__func__) + ")\n";
      throw std::runtime_error(err_msg);
    }
  }
  params.mdealias = dap;

  /* get the satellite (mandatory) */
  dso::SATELLITE sat(dso::SATELLITE::UNKNOWN);
  if (is_nonempty_string(config["satellite-attitude"]["satellite"])) {
    const auto s = config["satellite-attitude"]["satellite"].as<std::string>();
    try {
      sat = dso::translate_satid(s.c_str());
    } catch (std::exception &e) {
      fprintf(stderr, e.what());
      std::string err_msg = "[ERROR] Failed to translate satellite name "
                            "from "
                            "config parameters " +
                            s + " (traceback:" + std::string(__func__) + ")\n";
      throw std::runtime_error(err_msg);
    }
  }

  /* Satellite attitude (if model satellite-attitude::data_file is non-empty) */
  dso::SatelliteAttitude *satat = nullptr;
  dso::attitude_details::MeasuredAttitudeData *mad = nullptr;
  if (is_nonempty_string(config["satellite-attitude"]["data_file"])) {
    const auto d = config["satellite-attitude"]["data_file"].as<std::string>();
    try {
      satat = new dso::MeasuredAttitude(sat, d.c_str());
      mad = new dso::attitude_details::MeasuredAttitudeData(
          dso::measured_attitude_data_factory(sat));
    } catch (std::exception &e) {
      fprintf(stderr, e.what());
      std::string err_msg = "[ERROR] Failed to construct a SatelliteAttitude "
                            "and/or MeasuredAttitudeData instance from "
                            "config parameters " +
                            d + " (traceback:" + std::string(__func__) + ")\n";
      throw std::runtime_error(err_msg);
    }
  }
  params.matt = satat;   // attitude
  params.mattdata = mad; // attitude data (for retrieving attitude from stream)

  /* create the macromodel (independent of attitude) */
  dso::SatelliteMacromodel *satmm = nullptr;
  {
    try {
      satmm = new dso::SatelliteMacromodel(
          dso::SatelliteMacromodel::createSatelliteMacromodel(sat));
    } catch (std::exception &e) {
      fprintf(stderr, e.what());
      std::string err_msg = "[ERROR] Failed to construct a SatelliteMacromodel"
                            " instance from config parameters (traceback:" +
                            std::string(__func__) + ")\n";
      throw std::runtime_error(err_msg);
    }
    if (is_nonempty_string(config["satellite-attitude"]["cnes_sat_file"])) {
      const auto c =
          config["satellite-attitude"]["cnes_sat_file"].as<std::string>();
      if (satmm->load_satellite_mass_correction(c.c_str(), tt_midepoch)) {
        std::string err_msg = "[ERROR] Failed loading CNES satellite file " +
                              c + " (traceback:" + std::string(__func__) +
                              ")\n";
        throw std::runtime_error(err_msg);
      }
    }
  }
  params.msatmm = satmm; // macromodel (maybe used or not depending on attitude)

  /* atmospheric density */
  {
    if (is_nonempty_string(config["atmospheric-density-model"]["model"])) {
      if (is_nonempty_string(config["space-weather-data"]["celestrak_csv"])) {
        const auto swd =
            config["space-weather-data"]["celestrak_csv"].as<std::string>();
        const auto t1 =
            tt_start.add_seconds(dso::FractionalSeconds(-5 * 86400));
        const auto t2 = tt_stop.add_seconds(dso::FractionalSeconds(5 * 86400));
        std::vector<dso::SpaceWeatherData> *tswdata =
            new std::vector<dso::SpaceWeatherData>();
        *tswdata = dso::load_celestrak_sw(swd.c_str(), t1, t2);
        params.mspwdata = tswdata;

        params.matmdens = new dso::Nrlmsise00;
        params.matmdens->setup_data_hunter(*(params.mspwdata));
      } else {
        std::string err_msg = "[ERROR] Need Space Weather data to compute "
                              "atmospheric density; field is empty!"
                              " (traceback:" +
                              std::string(__func__) + ")\n";
        throw std::runtime_error(err_msg);
      }
    }
  }

  /* Dynamic parameters */
  {
    params.mCr = config["dynamic-parameters"]["Cr"].as<double>();
    params.mCd = config["dynamic-parameters"]["Cd"].as<double>();
  }

  /* return */
  return params;
}
