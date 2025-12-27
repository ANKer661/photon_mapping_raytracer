#ifndef Light_h
#define Light_h

#include "glm/glm.hpp"

/**
 Light class
 */
class Light {
public:
  glm::vec3 position; ///< Position of the light source
  glm::vec3 color;    ///< Color/intentisty of the light source
  bool area_light = false;
  Light(glm::vec3 position) : position(position) { color = glm::vec3(1.0); }
  Light(glm::vec3 position, glm::vec3 color)
      : position(position), color(color) {}
};

class AreaLight : public Light {
public:
  glm::vec3 u_vec;  
  glm::vec3 v_vec;
  glm::vec3 normal; // normal vector (given by u x v)
  glm::vec3 corner; // lower left corner: center - (u+v)/2
  /**
   * @param position Center position of the area light
   * @param color Color/intensity of the area light
   * @param u_vec Vector along width edge
   * @param v_vec Vector along height edge
   */
  AreaLight(glm::vec3 position, glm::vec3 u_vec, glm::vec3 v_vec,
            glm::vec3 color)
      : Light(position, color), u_vec(u_vec), v_vec(v_vec) {
    normal = -glm::normalize(glm::cross(u_vec, v_vec));  // we are left-handed
    corner = position - 0.5f * (u_vec + v_vec);
    this->area_light = true;
    float area = glm::length(glm::cross(u_vec, v_vec));
    this->color = color / area;
  }

  // not used currently
  void sample_light(int samples_per_axis,
                    std::vector<glm::vec3> &sample_positions) {
    float du = 1.0f / samples_per_axis;
    float dv = 1.0f / samples_per_axis;
    for (int i = 0; i < samples_per_axis; i++) {
      for (int j = 0; j < samples_per_axis; j++) {
        float su = (i + static_cast<float>(rand()) / RAND_MAX) * du;
        float sv = (j + static_cast<float>(rand()) / RAND_MAX) * dv;
        glm::vec3 sample_pos = corner + su * u_vec + sv * v_vec;
        sample_positions.push_back(sample_pos);
      }
    }
  }
};

#endif /* Light_h */