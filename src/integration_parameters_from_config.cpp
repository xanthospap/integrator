#include "integration_parameters.hpp"
#include "yaml-cpp/yaml.h"

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
  auto tt_midepoch = tt_start;
  tt_midepoch.add_seconds(dso::FractionalSeconds(dsec / 2.));

  /* fill EOP Series of the IntegrationParameters instance */
  {
    std::string tmp = config["eop"].as<std::string>();
    auto t1 = tt_start;
    t1.add_seconds(dso::seconds(-86400));
    auto t2 = tt_stop;
    t2.add_seconds(dso::seconds(86400));
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
    if (icgem.parse_data(DEGREE, ORDER, tt_midepoch, params.mgrav)) {
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
  if (config["solid-earth-tide"]["model"]) {
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
  if (config["ocean-tide"]["model"]) {
    /* basic parameters, must exist */
    tmp = config["ocean-tide"]["model"].as<std::string>();
    int DEGREE = config["ocean-tide"]["degree"].as<int>();
    int ORDER = config["ocean-tide"]["order"].as<int>();
    std::string data_dir = config["ocean-tide"]["data_dir"].as<std::string>();
    std::string groops_file_list =
        config["ocean-tide"]["groops_file_list"].as<std::string>();
    /* OceanTide model with admittance */
    if (config["ocean-tide"]["groops_doodson02_file"] &&
        config["ocean-tide"]["groops_admittance03_file"]) {
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
  if (config["solid-earth-pole-tide"]["model"]) {
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
  if (config["ocean-pole-tide"]["model"]) {
    tmp = config["ocean-pole-tide"]["model"].as<std::string>();
    if (tmp == "Desai02") {
      /* Desai 2002 model, need degree, order and coeffs */
      int DEGREE = config["ocean-pole-tide"]["degree"].as<int>();
      int ORDER = config["ocean-pole-tide"]["order"].as<int>();
      std::string fn = config["ocean-pole-tide"]["coeffs"].as<std::string>();
      opt = new dso::OceanPoleTide(DEGREE, ORDER, fn.c_str());
    } else {
      std::string err_msg = "[ERROR] Unknown Ocean Pole Tide model given " +
                            tmp + " (traceback:" + std::string(__func__) +
                            ")\n";
      throw std::runtime_error(err_msg);
    }
  }
  params.mop_tide = opt;

  /* Atmospheric Tide (if model name is non-empty) */
  dso::AtmosphericTide *atm = nullptr;
  if (config["atmospheric-tide"]["model"]) {
    int DEGREE = config["atmospheric-tide"]["degree"].as<int>();
    int ORDER = config["atmospheric-tide"]["order"].as<int>();
    std::string data_dir =
        config["atmospheric-tide"]["data_dir"].as<std::string>();
    if (config["atmospheric-tide"]["groops_file_list"]) {
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
  if (config["dealiasing"]["model"]) {
    // tmp = config["dealiasing"]["model"].as<std::string>();
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

  /* return */
  return params;
}