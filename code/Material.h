//
//  Material.h
//  Raytracer
//
//  Created by Piotr Didyk on 14.07.21.
//

#ifndef Material_h
#define Material_h

#include "glm/glm.hpp"

/**
 Structure describing a material of an object
 */
struct Material {
  glm::vec3 ambient = glm::vec3(0.0);
  glm::vec3 diffuse = glm::vec3(1.0);
  glm::vec3 specular = glm::vec3(0.0);
  float shininess = 0.0;
  float reflect_portion = 0.0; // portion of perfect reflection
  bool refractive = false;
  float refractive_index = 1.0;
  bool fresnel = false; // if true, determine reflect portion dynmically
  glm::vec3 emissive = glm::vec3(0.0); // for self-emission
};


/**
 * Function to compute the perfect reflection portion using Fresnel equations.
 * @param incident The incident direction
 * @param refract_normal The normal at the refraction point
 * @param refract_direction The refracted direction
 * @param material The material of the object
 * @param outside Boolean indicating whether the ray is outside the object
 */
float fresnel_perfect_reflect(const glm::vec3 &incident,
                              const glm::vec3 &refract_normal,
                              const glm::vec3 &refract_direction,
                              const Material &material, bool outside) {
  float idx1, idx2;
  float cos1 = -glm::dot(incident, refract_normal);
  float cos2 = -glm::dot(refract_direction, refract_normal);
  if (outside) {
    idx1 = 1.0f;
    idx2 = material.refractive_index;
  } else {
    idx1 = material.refractive_index;
    idx2 = 1.0f;
  }
  float i1c1 = idx1 * cos1;
  float i2c2 = idx2 * cos2;
  float i1c2 = idx1 * cos2;
  float i2c1 = idx2 * cos1;

  return 0.5 * (pow((i1c1 - i2c2) / (i1c1 + i2c2), 2) +
                pow((i1c2 - i2c1) / (i1c2 + i2c1), 2));
}

#endif /* Material_h */
