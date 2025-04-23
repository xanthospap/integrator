#ifndef __DSO_POD_INTEGRATION_PARAMETERS_HPP__
#define __DSO_POD_INTEGRATION_PARAMETERS_HPP__

#include "datetime/calendar.hpp"
#include "iers/aod1b_data_stream.hpp"
#include "iers/atmospheric_tides.hpp"
#include "iers/ocean_tide.hpp"
#include "iers/planets.hpp"
#include "iers/pole_tide.hpp"
#include "iers/solid_earth_tide.hpp"

namespace dso {

class IntegrationParameters {
  static constexpr const int PI_MAX_DEGREE = 200;
  static constexpr const int PI_MAX_ORDER = 200;

public:
  /* TAI epoch */
  MjdEpoch mtai0;
  /* EOPs */
  EopSeries meops;
  /* Earth's gravity field (Stokes) */
  StokesCoeffs mgrav;
  /* tidal phenomena ... */
  SolidEarthTide *mse_tide{nullptr};
  OceanTide *moc_tide{nullptr};
  PoleTide *mep_tide{nullptr};
  OceanPoleTide *mop_tide{nullptr};
  AtmosphericTide *mat_tide{nullptr};
  /* dealiasing */
  Aod1bDataStream<AOD1BCoefficientType::GLO> *mdealias{nullptr};

  /* third body gravity
   *                       | Bit Nr. | Default state
   * Bit sequence: Moon    | (0)     | 1
   *               Sun     | (1)     | 1
   *               Mercury | (2)     | 0
   *               Venus   | (3)     | 0
   *               Mars    | (4)     | 0
   *               Jupiter | (5)     | 0
   *
   * If you want a planet other than Moon/Sun to be included in the
   * computation, toggle the relevant bit. E.g. to set Mars to `on`, use:
   * tbg |= (1 << 4); (in general: tbg |= (1 << Bit_Nr)
   * To check if a planet will be included:
   * bool included = tbg & (1 << Bit_Nr);
   * To exclude a planet:
   * tbg &= ~(1 << Bit_Nr);
   */
  uint8_t tbg = 0x03;

  /* allocate scratch space for computations */
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mW{PI_MAX_DEGREE + 3,
                                                           PI_MAX_DEGREE + 3};
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> mM{PI_MAX_DEGREE + 3,
                                                           PI_MAX_DEGREE + 3};

public:
  /** @brief Constructor */
  IntegrationParameters() noexcept {};

  /** @brief No copy constructor */
  IntegrationParameters(const IntegrationParameters &) = delete;

  /** @brief No assignment operator */
  IntegrationParameters &operator=(const IntegrationParameters &) = delete;

  /** @brief Move constructor */
  IntegrationParameters(IntegrationParameters &&other) noexcept
      : mtai0(other.mtai0), meops(std::move(other.meops)),
        mgrav(std::move(other.mgrav)), mse_tide(other.mse_tide),
        moc_tide(other.moc_tide), mep_tide(other.mep_tide),
        mop_tide(other.mop_tide), mat_tide(other.mat_tide), tbg(other.tbg),
        mW(std::move(other.mW)), mM(std::move(other.mM)) {

    other.mse_tide = nullptr;
    other.moc_tide = nullptr;
    other.mep_tide = nullptr;
    other.mop_tide = nullptr;
    other.mat_tide = nullptr;
  }

  /** @brief Move assignment operator */
  IntegrationParameters &operator=(IntegrationParameters &&other) noexcept {
    if (this != &other) {
      mtai0 = other.mtai0;
      meops = std::move(other.meops);
      mgrav = std::move(other.mgrav);
      mse_tide = other.mse_tide;
      moc_tide = other.moc_tide;
      mep_tide = other.mep_tide;
      mop_tide = other.mop_tide;
      mat_tide = other.mat_tide;
      tbg = other.tbg;
      mW = std::move(other.mW);
      mM = std::move(other.mM);
      other.mse_tide = nullptr;
      other.moc_tide = nullptr;
      other.mep_tide = nullptr;
      other.mop_tide = nullptr;
      other.mat_tide = nullptr;
    }
    return *this;
  }

  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &tw() noexcept {
    return mW;
  }
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &tm() noexcept {
    return mM;
  }

  const MjdEpoch &t0() const noexcept { return mtai0; }
  MjdEpoch &t0() noexcept { return mtai0; }

  /** @brief Construct instance from a YAML configuration file. */
  static IntegrationParameters from_config(const char *config_fn,
                                           const MjdEpoch &tt_start,
                                           const MjdEpoch &tt_stop);
}; /* class IntegrationParameters */

} /* namespace dso */

#endif
