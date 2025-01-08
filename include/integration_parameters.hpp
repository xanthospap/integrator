#ifndef __DSO_POD_INTEGRATION_PARAMETERS_HPP__
#define __DSO_POD_INTEGRATION_PARAMETERS_HPP__

#include "datetime/calendar.hpp"
#include "iers/atmospheric_tides.hpp"
#include "iers/aod1b_data_stream.hpp"
#include "iers/ocean_tide.hpp"
#include "iers/pole_tide.hpp"
#include "iers/solid_earth_tide.hpp"
#include "iers/planets.hpp"

namespace dso {

class IntegrationParameters {
  static constexpr const int PI_MAX_DEGREE = 200;
  static constexpr const int PI_MAX_ORDER  = 200;

// private:
public:
  /* TAI epoch */
  MjdEpoch mtai0;
  /* EOPs */
  EopSeries *meops;
  /* Earth's gravity field */
  StokesCoeffs *mgrav;
  /* tidal phenomena ... */
  SolidEarthTide  *mse_tide{nullptr};
  OceanTide       *moc_tide{nullptr};
  PoleTide        *mep_tide{nullptr};
  OceanPoleTide   *mop_tide{nullptr};
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
  EopSeries *eops() noexcept { return meops; }
  StokesCoeffs *grav() noexcept { return mgrav; }
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &tw() noexcept {
    return mW;
  }
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> &tm() noexcept {
    return mM;
  }
  MjdEpoch t0() const noexcept { return mtai0; }
}; /* class IntegrationParameters */

} /* namespace dso */

#endif
