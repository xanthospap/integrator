# integrator
Integrator for POD

## Installation

You will need:
 * the `datetime` library by [DSO](https://github.com/DSOlab/ggdatetime),
 * the `iers` library by [DSO](https://github.com/xanthospap/iers2010),
 * the `geodesy` library dy [DSO](https://github.com/DSOlab/ggeodesy),
 * the `sp3` library by [DSO](https://github.com/xanthospap/sp3),
 * the [yaml-cpp](https://github.com/jbeder/yaml-cpp) library for parsing `yaml` (configuration) files
 * the [cspice](https://naif.jpl.nasa.gov/naif/toolkit.html) library (if not already installed with `iers`). For instructions, see the [iers repository](https://github.com/xanthospap/iers2010)

Installation is performed via [cmake](https://cmake.org/):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/usr/local/lib
cmake --build build --target all --config=Release -- -j4
cd build && sudo make install
```
