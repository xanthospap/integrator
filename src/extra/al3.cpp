#include <cmath>

/** @brief State vector to classical orbital elements.
 *
 * Classical orbit elements is the set of {α, e, i, Ω, ω, θ} (semimajor axis, 
 * eccentricity, inclination, longitude of the ascending node, argument of 
 * pericenter and true anomaly).
 *
 * This algorithm follows the approach described in Flores and Fantino (2024).
 * In the publication, the algorithm is named AL3.
 *
 * The approach here is not entirely general; it will not work for zero 
 * orbital angular momentum. However, this case is of minimal relevance in 
 * practical applications (the trajectory is a straight line passing through 
 * the primary) and can usually be ignored.
 *
 * @param[in] r  Position vector of satellite, in an inertial Cartesian 
 *               frame (ECI) in [m].
 * @param[in] v  Velocity vector of satellite, in an inertial Cartesian frame 
 *               (ECI) in [m/sec].
 * @param[in] mu Gravitational parameter, μ=GM in [m**3 s**−2].
 *
 * R. Flores, E. Fantino, COMPUTING CLASSICAL ORBITAL ELEMENTS WITH IMPROVED 
 * EFFICIENCY AND ACCURACY, Advances in Space Research, 2024, ISSN 0273-1177,
 * https://doi.org/10.1016/j.asr.2024.12.048.
 */
int sv2coe(const Eigen::Matrix<double, 3, 1> &r,
           const Eigen::Matrix<double, 3, 1> &v, double mu) noexcept {

  /* step 1 */
  const Eigen::Vector3d h = r.cross(v);
  const double i = std::atan2(std::sqrt(h(0) * h(0) + h(1) * h(1)), h(2));
  const double Omega = std::atan2(h(0), -h(1));

  /* step 2 (note that b and n are unit vectors) */
  Eigen::Vector3d n;
  n << std::cos(Omega), std::sin(Omega), 0e0;
  const Eigen::Vector3d b = (h.normalized()).cross(n);

  /* step 3 */
  const double a = 1e0 / (2e0 / r.norm() - v.dot(v) / mu);
  const Eigen::Vector3d e = v.cross(h) / mu - r.normalized();

  /* step 4 (note that theta is true anomaly) */
  const double omega = std::atan2(e.dot(b), e.dot(n));
  const double theta = std::atan2(r.dot(b), r.dot(n)) - omega;

  /* all done */
  return;
}
