#ifndef __DSO_POD_INTEGRATION_PARAMETERS_HPP__
#define __DSO_POD_INTEGRATION_PARAMETERS_HPP__

#include "datetime/calendar.hpp"
#include "iers/atmospheric_tides.hpp"
#include "iers/aod1b_data_stream.hpp"
#include "iers/ocean_tide.hpp"
#include "iers/pole_tide.hpp"
#include "iers/solid_earth_tide.hpp"

namespace dso {

class IntegrationParameters {
  static constexpr const int PI_MAX_DEGREE = 200;
  static constexpr const int PI_MAX_ORDER  = 200;

private:
  /* TAI epoch */
  MjdEpoch tai0;
  /* Earth's gravity field */
  StokesCoeffs stokes;
  /* tidal phenomena ... */
  SolidEarthTide  *se_tide{nullptr};
  OceanTide       *oc_tide{nullptr};
  PoleTide        *ep_tide{nullptr};
  OceanPoleTide   *op_tide{nullptr};
  AtmosphericTide *at_tide{nullptr};
  /* dealiasing */
  Aod1bDataStream<AOD1BCoefficientType::GLO> *dealias{nullptr};

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
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> W{PI_MAX_DEGREE + 3,
                                                          PI_MAX_DEGREE + 3};
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> M{PI_MAX_DEGREE + 3,
                                                          PI_MAX_DEGREE + 3};

public:
}; /* class IntegrationParameters */

} /* namespace dso */

#endif
