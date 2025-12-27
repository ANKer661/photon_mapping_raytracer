#ifndef Sampling_h
#define Sampling_h

#include "glm/glm.hpp"
#include <cmath>

/** Build an orthonormal basis (T, B, N) given a normal vector (Pixar)
 * @param normal The input normal vector
 * @param T Output tangent vector
 * @param B Output bitangent vector
 */
inline void build_orthonormal_basis(const glm::vec3 &normal, glm::vec3 &T,
                                    glm::vec3 &B) {
  float sign = copysign(1.0f, normal.x);
  const float a = -1.0f / (sign + normal.z);
  const float b = normal.x * normal.y * a;
  T = glm::vec3(1.0f + sign * normal.x * normal.x * a, sign * b,
                -sign * normal.x);
  B = glm::vec3(b, sign + normal.y * normal.y * a, -normal.y);
}

/** Sample a random direction in hemisphere with cosine weighting
 * @param normal The normal vector defining the hemisphere orientation
 * @return A sampled direction vector
 */
inline glm::vec3 sample_hemisphere(const glm::vec3 &normal) {
  glm::vec3 T;
  glm::vec3 B;
  build_orthonormal_basis(normal, T, B);

  float u = static_cast<float>(rand()) / RAND_MAX;
  float v = static_cast<float>(rand()) / RAND_MAX;
  float z = sqrt(u);
  float r = sqrt(1.0f - u);
  float phi = 2.0f * M_PI * v;

  return glm::normalize(r * cos(phi) * T + r * sin(phi) * B + z * normal);
}

/** Overloaded version of sample_hemisphere with given T and B
 * @param normal The normal vector defining the hemisphere orientation
 * @param T Tangent vector
 * @param B Bitangent vector
 * @return A sampled direction vector
 */
inline glm::vec3 sample_hemisphere(const glm::vec3 &normal, const glm::vec3 &T,
                                   const glm::vec3 &B) {
  float u = static_cast<float>(rand()) / RAND_MAX;
  float v = static_cast<float>(rand()) / RAND_MAX;
  float z = sqrt(u);
  float r = sqrt(1.0f - u);
  float phi = 2.0f * M_PI * v;

  return glm::normalize(r * cos(phi) * T + r * sin(phi) * B + z * normal);
}

/** Sample a random direction uniformly on a sphere
 * used for sampling a point light direction
 * @return A sampled direction vector
 */
inline glm::vec3 sample_uniform_sphere() {
  float u = static_cast<float>(rand()) / RAND_MAX;
  float v = static_cast<float>(rand()) / RAND_MAX;
  float z = 1.0f - 2.0f * v;
  float theta = 2.0f * M_PI * u;
  float r = sqrt(1.0f - z * z);
  return glm::vec3(r * cos(theta), r * sin(theta), z);
}

#endif /* Sampling_h */