#include "dop853.hpp"
#include "sp3.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s SP3_file [SAT_ID]\n", argv[0]);
    return 1;
  }

  /* create an Sp3c instance */
  dso::Sp3c sp3(argv[1]);

  /* choose satellite */
  dso::sp3::SatelliteId sv = sp3.sattellite_vector()[0];
  if (argc == 3) {
    /* check if the satellite is included in the Sp3 */
    if (!sp3.has_sv(dso::sp3::SatelliteId(argv[2]))) {
      fprintf(stderr, "Error. Satellite [%s] not included in sp3 file!\n", argv[2]);
      return 1;
    }
    sv = dso::sp3::SatelliteId{argv[2]};
  }

  /* get states */
  // auto start_t = sp3.start_epoch();
  dso::Sp3DataBlock block;

  int sp3err=0;
  while (!sp3err) {
    sp3err = sp3.get_next_data_block(sv, block);
    if (sp3err > 0) {
      printf("Something went wrong ....status = %3d\n", sp3err);
      return 1;
    } 
    bool position_ok = !block.flag.is_set(dso::Sp3Event::bad_abscent_position);
    if (position_ok && (!sp3err))
      printf("%15.6f %15.7f %15.7f %15.7f\n",
             block.t.imjd().as_underlying_type() +
                 block.t.fractional_days().days(),
             block.state[0], block.state[1], block.state[2]);
  }

  /* load planetary ephemeris kernels */
  dso::load_spice_kernel(argv[?]);
  dso::load_spice_kernel(argv[?]);

  /* Earth's gravity field */
  dso::StokesCoeffs stokes;
  {
  dso::Icgem icgem(argv[?]);
  if (icgem.parse_data(DEGREE, ORDER, start_t, stokes)) {
    fprintf(stderr, "ERROR Failed reading gravity model!\n");
    return 1;
  }
  /* checks */
  assert(stokes.max_degree() == DEGREE);
  assert(stokes.max_order() == ORDER);
  }

  /* solid earth tide */
  dso::SolidEarthTide se_tide;

  return sp3err;
}
