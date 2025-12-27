/**
@file main.cpp
*/

#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include "Image.h"
#include "Light.h"
#include "Material.h"
#include "Object.h"
#include "Photon.h"
#include "Ray.h"
#include "Sampling.h"

using namespace std;

vector<Object *> objects; ///< A list of all objects in the scene

Hit find_closest_hit(const Ray &ray) {
  Hit closest_hit;
  closest_hit.hit = false;
  closest_hit.distance = INFINITY;

  for (Object *obj : objects) {
    Hit hit = obj->intersect(ray);
    if (hit.hit && hit.distance < closest_hit.distance) {
      closest_hit = hit;
    }
  }
  return closest_hit;
}

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(0.001);

bool visualize_photons = false;
bool visualize_photonmap = false;
int max_recursion_depth = 5;
float secondary_roulette_prob = 0.75f;
int area_light_samples_per_axis = 16;
float caustic_radius = visualize_photons ? 0.01f : 0.2f;
float global_radius = visualize_photons ? 0.01f : 0.6f;
bool super_sampling = false;
int num_photons = 8'000'000;

// photon mapping instance
PhotonMap global_map;  // for indirect photons
PhotonMap caustic_map; // for caustic photons

// to use this function in PhongModel
glm::vec3 trace_ray(Ray ray, int depth);

void trace_photon(Ray &ray, glm::vec3 power, int depth, bool caustic);

/** Function for emitting photons into the scene
 * @param num_photons Number of photons to be emitted of each light source
 * @param total_power Total power of each light source
 */
void emit_photons(int num_photons, float total_power) {
  for (int light_id = 0; light_id < lights.size(); light_id++) {
    if (lights[light_id]->area_light) {
      // area light
      AreaLight *area_light = static_cast<AreaLight *>(lights[light_id]);
      glm::vec3 normal = area_light->normal;
      glm::vec3 T;
      glm::vec3 B;
      build_orthonormal_basis(normal, T, B);
      for (int k = 0; k < num_photons; k++) {
        float r1 = static_cast<float>(rand()) / RAND_MAX;
        float r2 = static_cast<float>(rand()) / RAND_MAX;
        glm::vec3 position = area_light->corner + area_light->u_vec * r1 +
                             area_light->v_vec * r2;
        // sample random direction in hemisphere
        glm::vec3 dir = sample_hemisphere(normal, T, B);
        // emit the ray
        Ray photon_ray(position, dir);
        glm::vec3 power = (total_power / num_photons) *
                          lights[light_id]->color; // power per photon
        // trace the photon
        trace_photon(photon_ray, power, 0, false);
      }
    } else {
      // point light
      for (int i = 0; i < num_photons; i++) {
        // sample random direction
        glm::vec3 dir = sample_uniform_sphere();
        // emit the ray
        Ray photon_ray(lights[light_id]->position, dir);
        glm::vec3 power = (total_power / num_photons) *
                          lights[light_id]->color; // power per photon
        // trace the photon
        trace_photon(photon_ray, power, 0, false);
      }
    }
  }

  // building photon maps
  global_map.build();
  caustic_map.build();
}

void trace_photon(Ray &ray, glm::vec3 power, int depth, bool caustic) {
  if (depth > max_recursion_depth)
    return;

  // find nearest intersection
  Hit closest_hit;
  closest_hit.hit = false;
  closest_hit.distance = INFINITY;

  for (Object *obj : objects) {
    if (glm::length(obj->getMaterial().emissive) > 0.0f) {
      continue; // self-emissive object does not interact with photons
    }

    Hit hit = obj->intersect(ray);
    if (hit.hit && hit.distance < closest_hit.distance) {
      closest_hit = hit;
    }
  }
  if (!closest_hit.hit) {
    return;
  }

  Material material = closest_hit.object->getMaterial();

  // handle fresnel effect: dynamic reflect portion
  // handle reflect portion: russian roulette -> perfect reflection or
  // refraction/diffuse
  float reflect_portion = material.reflect_portion;
  if (material.fresnel) {
    glm::vec3 incident = ray.direction;
    bool outside = glm::dot(closest_hit.normal, incident) < 0.0f;
    glm::vec3 refract_normal =
        outside ? closest_hit.normal : -closest_hit.normal;
    float eta = outside ? (1.0f / material.refractive_index)
                        : material.refractive_index;
    glm::vec3 refract_direction =
        glm::normalize(glm::refract(incident, refract_normal, eta));
    reflect_portion = fresnel_perfect_reflect(
        incident, refract_normal, refract_direction, material, outside);
  }

  float roulette = static_cast<float>(rand()) / RAND_MAX;
  if (roulette < reflect_portion) {
    // perfect reflection
    glm::vec3 reflect_direction =
        glm::reflect(-ray.direction, closest_hit.normal);
    Ray reflect_ray(closest_hit.intersection + 0.001f * reflect_direction,
                    reflect_direction);
    trace_photon(reflect_ray, power, depth + 1, true);
  } else {
    // refraction or diffuse
    if (material.refractive) {
      // refraction
      glm::vec3 incident = ray.direction;
      bool outside = glm::dot(closest_hit.normal, incident) < 0.0f;
      glm::vec3 refract_normal =
          outside ? closest_hit.normal : -closest_hit.normal;
      float eta = outside ? (1.0f / material.refractive_index)
                          : material.refractive_index;
      glm::vec3 refract_direction =
          glm::normalize(glm::refract(incident, refract_normal, eta));
      Ray refract_ray(closest_hit.intersection + 0.001f * refract_direction,
                      refract_direction);
      trace_photon(refract_ray, power, depth + 1, true);
    } else {
      // store the photon in caustic or global map
      if (depth > 0) { // exclude direct photons
        Photon photon;
        photon.position = closest_hit.intersection;
        photon.direction = ray.direction; // pointing to intersection
        photon.color = power;
        if (caustic) {
          caustic_map.store(photon);
        } else {
          global_map.store(photon);
        }
      }

      // create secondary ray (russian roulette)
      roulette = static_cast<float>(rand()) / RAND_MAX;
      if (depth == 0 || roulette < secondary_roulette_prob) {
        glm::vec3 new_power = power * material.diffuse;
        if (depth != 0) {
          new_power = new_power * (1.0f / secondary_roulette_prob);
        }
        // sample random direction in hemisphere (cosine weighted)
        glm::vec3 new_direction = sample_hemisphere(closest_hit.normal);
        Ray new_ray(closest_hit.intersection + 0.001f * new_direction,
                    new_direction);
        trace_photon(new_ray, new_power, depth + 1, false);
      }
    }
  }
}

glm::vec3 compute_light_contribution(glm::vec3 &light_pos,
                                     glm::vec3 &light_color, glm::vec3 &point,
                                     glm::vec3 &normal,
                                     glm::vec3 &view_direction,
                                     Material &material,
                                     float perfect_reflect_portion) {
  glm::vec3 light_direction = glm::normalize(light_pos - point);
  glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

  // detect if the light source is blocked by other object
  glm::vec3 shadow_origin = point + 0.001f * light_direction;
  Ray shadow_ray(shadow_origin, light_direction);
  float d = glm::distance(light_pos, shadow_origin);
  bool blocked = false;
  bool is_transparent = true; // trick: if block by refractive, scale down but
                              // not completely block
  for (int k = 0; k < objects.size(); k++) {
    // self-emissive object does not block
    if (glm::length(objects[k]->getMaterial().emissive) > 0.0f) {
      continue;
    }
    Hit hit = objects[k]->intersect(shadow_ray);
    if (hit.hit == true && hit.distance < d) {
      blocked = true;
      Material hit_material = hit.object->getMaterial();
      is_transparent &=
          hit_material.refractive; // only if all blocking objects are
                                   // refractive, we consider it transparent
      break;
    }
  }

  // return in advance if completely blocked
  if (blocked && !is_transparent) {
    return glm::vec3(0.0);
  }

  float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);

  glm::vec3 diffuse_color = material.diffuse;
  glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
  glm::vec3 specular(0.0);
  if (material.shininess > 0.0f) {
    float VdotR =
        glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);
    specular += material.specular * glm::vec3(pow(VdotR, material.shininess));
  }

  float r = glm::distance(point, light_pos);
  r = max(r, 0.1f);
  return (light_color * (diffuse + specular) / r / r) *
         (1.0f - perfect_reflect_portion) * (blocked ? 0.75f : 1.0f)
         + material.emissive; // add emissive term
}

float cone_normalization = 3.0f;
/** Function for adding photon mapping contribution with cone filtering
 @param photons A list of photons to be considered
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @return Color contribution from the photons
*/
glm::vec3 add_photon_mapping(vector<Photon> &photons, float radius,
                             glm::vec3 &point, glm::vec3 &normal,
                             Material material, glm::vec3 &view_direction) {

  glm::vec3 photon_color(0.0);
  // lambertian brdf
  glm::vec3 brdf = material.diffuse / float(M_PI);
  for (Photon &photon : photons) {
    float distance_weight = 1 - glm::distance(point, photon.position) / radius;
    photon_color += photon.color * distance_weight;
  }

  if (!visualize_photons) { // only show photon color, no brdf
    photon_color *= brdf;
  }
  return photon_color * cone_normalization / float(M_PI * radius * radius);
}

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param view_direction A normalized direction from the point to the
 viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal,
                     glm::vec3 view_direction, Material material,
                     int depth = 0) {

  glm::vec3 color(0.0);

  float perfect_reflect_portion = material.reflect_portion;

  if (material.refractive) {
    glm::vec3 refract_normal;
    float eta; // ratio of index, incident / transmitted
    glm::vec3 incident = -view_direction;
    bool outside = glm::dot(normal, incident) < 0.0f;
    if (outside) {
      refract_normal = normal;
      eta = 1.0f / material.refractive_index;
    } else {
      // if inside, reverse normal
      refract_normal = -normal;
      eta = material.refractive_index;
    }

    glm::vec3 refract_vec = glm::refract(incident, refract_normal, eta);
    // handle total internal reflection
    if (glm::length(refract_vec) < 0.001f) {
      perfect_reflect_portion = 1.0f;
    } else {
      glm::vec3 refract_direction = glm::normalize(refract_vec);
      glm::vec3 refract_origin = point + 0.001f * refract_direction;
      Ray refract_ray(refract_origin, refract_direction);
      glm::vec3 refract_color = trace_ray(refract_ray, depth + 1);
      // handle fresnel effect
      if (material.fresnel) {
        perfect_reflect_portion = fresnel_perfect_reflect(
            incident, refract_normal, refract_direction, material, outside);
      }
      color += (1.0f - perfect_reflect_portion) * refract_color;
    }
  } else {
    for (int light_num = 0; light_num < lights.size(); light_num++) {
      if (lights[light_num]->area_light) {
        // area light
        AreaLight *area_light = static_cast<AreaLight *>(lights[light_num]);
        int total_samples =
            area_light_samples_per_axis * area_light_samples_per_axis;
        // sample at multiple sub-positions
        float du = 1.0f / area_light_samples_per_axis;
        float dv = 1.0f / area_light_samples_per_axis;
        for (int i = 0; i < area_light_samples_per_axis; i++) {
          for (int j = 0; j < area_light_samples_per_axis; j++) {
            float su = (i + static_cast<float>(rand()) / RAND_MAX) * du;
            float sv = (j + static_cast<float>(rand()) / RAND_MAX) * dv;
            glm::vec3 sample_pos = area_light->corner + su * area_light->u_vec +
                                   sv * area_light->v_vec;

            glm::vec3 light_normal = area_light->normal;
            float weight =
                glm::dot(light_normal, -glm::normalize(sample_pos - point));
            if (weight > 0.0f) { // exclude back face
              glm::vec3 light_color = area_light->color * weight;
              color += compute_light_contribution(
                           sample_pos, light_color, point, normal,
                           view_direction, material, perfect_reflect_portion) /
                       float(total_samples);
            }
          }
        }
      } else {
        // point light
        color += compute_light_contribution(
            lights[light_num]->position, lights[light_num]->color, point,
            normal, view_direction, material, perfect_reflect_portion);
      }
    }
  }

  // ambient
  color += ambient_light * material.ambient * (1.0f - perfect_reflect_portion);

  // perfect reflection part
  if (perfect_reflect_portion > 0.0f && depth < max_recursion_depth) {
    glm::vec3 pr_direction = glm::reflect(-view_direction, normal);
    glm::vec3 pr_origin = point + 0.001f * pr_direction;
    Ray pr_ray(pr_origin, pr_direction);
    glm::vec3 pr_color = trace_ray(pr_ray, depth + 1);
    color += perfect_reflect_portion * pr_color;
  }

  // gathering photon mapping
  if (visualize_photonmap || visualize_photons) {
    color = glm::vec3(0.0);
  }
  std::vector<Photon> caustic_photons;
  caustic_map.find_photons(point, caustic_radius, caustic_photons);
  if (!caustic_photons.empty()) {
    color += add_photon_mapping(caustic_photons, caustic_radius, point, normal,
                                material, view_direction);
  }

  std::vector<Photon> global_photons;
  global_map.find_photons(point, global_radius, global_photons);
  if (!global_photons.empty()) {
    color += add_photon_mapping(global_photons, global_radius, point, normal,
                                material, view_direction);
  }

  color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
  return color;
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray, int depth = 0) {
  Hit closest_hit = find_closest_hit(ray);

  glm::vec3 color(0.0);
  if (closest_hit.hit) {
    color = PhongModel(closest_hit.intersection, closest_hit.normal,
                       glm::normalize(-ray.direction),
                       closest_hit.object->getMaterial(), depth);
  } else {
    color = glm::vec3(0.0, 0.0, 0.0);
  }
  return color;
}

/**
 Function for handling area light, add a light and a plane simultaneously
*/
void add_area_light(glm::vec3 position, glm::vec3 u_vec, glm::vec3 v_vec,
                    glm::vec3 color) {
  AreaLight *area_light = new AreaLight(position, u_vec, v_vec, color);
  lights.push_back(area_light);
  // add a plane to represent the area light
  Material light_material;
  light_material.emissive = area_light->color; // visualize to the light color
  Square *light_square =
      new Square(position, u_vec, v_vec, light_material);
  objects.push_back(light_square);
}

/**
 Function defining the scene
 */
void sceneDefinition() {

  Material green_diffuse;
  green_diffuse.ambient = glm::vec3(0.7f, 0.9f, 0.7f);
  green_diffuse.diffuse = glm::vec3(0.7f, 0.9f, 0.7f);

  Material red_specular;
  red_specular.ambient = glm::vec3(1.0f, 0.3f, 0.3f);
  red_specular.diffuse = glm::vec3(1.0f, 0.3f, 0.3f);
  red_specular.specular = glm::vec3(0.5);
  red_specular.shininess = 10.0;

  Material blue_specular;
  blue_specular.ambient = glm::vec3(0.7f, 0.7f, 1.0f);
  blue_specular.diffuse = glm::vec3(0.7f, 0.7f, 1.0f);
  blue_specular.specular = glm::vec3(0.6);
  blue_specular.shininess = 100.0;

  // Material green_diffuse;
  green_diffuse.ambient = glm::vec3(0.03f, 0.1f, 0.03f);
  green_diffuse.diffuse = glm::vec3(0.3f, 1.0f, 0.3f);

  // Material red_specular;
  red_specular.diffuse = glm::vec3(1.0f, 0.2f, 0.2f);
  red_specular.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
  red_specular.specular = glm::vec3(0.5);
  red_specular.shininess = 5.0;
  // red_specular.fresnel = true;

  // Material blue_specular;
  blue_specular.ambient = glm::vec3(0.02f, 0.02f, 0.1f);
  blue_specular.diffuse = glm::vec3(0.2f, 0.2f, 1.0f);
  blue_specular.specular = glm::vec3(0.6);
  blue_specular.shininess = 100.0;
  blue_specular.reflect_portion = 1.0f;

  // Materials for refractive
  Material refractive_material;
  refractive_material.refractive = true;
  refractive_material.refractive_index = 2.5f;
  refractive_material.fresnel = true;

  objects.push_back(new Sphere(2.0, glm::vec3(3, -0.6, 8), blue_specular));
  objects.push_back(new Sphere(1.0, glm::vec3(0, -2, 6), red_specular));
  objects.push_back(
      new Sphere(2.0, glm::vec3(-3, -0.6, 8), refractive_material));

  // add a area light
  // lights.push_back(new AreaLight(glm::vec3(10, 5, 2), glm::vec3(-3, 0, -3),
  //                                glm::vec3(-3, 4, 0), glm::vec3(10)));
  // lights.push_back(new AreaLight(glm::vec3(8, -1.5, 3), glm::vec3(-1, 0, -1),
  //                                glm::vec3(0, 3, 0.3), glm::vec3(0.005)));
  // lights.push_back(new AreaLight(glm::vec3(1, -1.5, 5), glm::vec3(-1, 0, -1),
  //                                glm::vec3(0, 3, 0.2), glm::vec3(0.005)));
  // lights.push_back(new AreaLight(glm::vec3(0, 26, 5), glm::vec3(3, 0, 0),
  //                                glm::vec3(0, 0, -3), glm::vec3(9)));
  add_area_light(glm::vec3(0, 26, 5), glm::vec3(3, 0, 0), glm::vec3(0, 0, -3),
                 glm::vec3(9));
  // lights.push_back(new Light(glm::vec3(0, 9, 12), glm::vec3(0.1)));
  lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.2)));
  // lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(0.8)));

  // test light
  glm::vec3 light_direction =
      glm::normalize(glm::vec3(0, 1, 12) - glm::vec3(0, -2, 6));
  glm::vec3 u;
  glm::vec3 v;
  build_orthonormal_basis(light_direction, u, v);
  // lights.push_back(
  //     new AreaLight(glm::vec3(0, 3, 12), u * 3.0f, v * 3.0f, glm::vec3(3)));
  add_area_light(glm::vec3(0, 3, 12), u * 3.0f, v * 3.0f, glm::vec3(3));

  Material red_diffuse;
  red_diffuse.ambient = glm::vec3(0.09f, 0.06f, 0.06f);
  red_diffuse.diffuse = glm::vec3(0.9f, 0.6f, 0.6f);

  Material blue_diffuse;
  blue_diffuse.ambient = glm::vec3(0.06f, 0.06f, 0.09f);
  blue_diffuse.diffuse = glm::vec3(0.6f, 0.6f, 0.9f);
  objects.push_back(new Plane(glm::vec3(0, -3, 0), glm::vec3(0.0, 1, 0)));
  objects.push_back(
      new Plane(glm::vec3(0, 1, 30), glm::vec3(0.0, 0.0, -1.0), green_diffuse));
  objects.push_back(
      new Plane(glm::vec3(-15, 1, 0), glm::vec3(1.0, 0.0, 0.0), red_diffuse));
  objects.push_back(
      new Plane(glm::vec3(15, 1, 0), glm::vec3(-1.0, 0.0, 0.0), blue_diffuse));
  objects.push_back(new Plane(glm::vec3(0, 27, 0), glm::vec3(0.0, -1, 0)));

  Material real_blue_diffuse;
  real_blue_diffuse.ambient = glm::vec3(0.03f, 0.03f, 0.09f);
  real_blue_diffuse.diffuse = glm::vec3(0.3f, 0.3f, 0.9f);
  objects.push_back(new Plane(glm::vec3(0, 1, -8), glm::vec3(0.0, 0.0, 1.0),
                              real_blue_diffuse));

  Cone *cone2 = new Cone(green_diffuse);
  glm::mat4 translationMatrix = glm::translate(glm::vec3(6, -3, 5));
  glm::mat4 scalingMatrix = glm::scale(glm::vec3(1.0f, 3.0f, 1.0f));
  glm::mat4 rotationMatrix = glm::rotate(glm::atan(3.0f), glm::vec3(0, 0, 1));
  cone2->setTransformation(translationMatrix * rotationMatrix * scalingMatrix);
  objects.push_back(cone2);
}

glm::vec3 toneMapping(glm::vec3 intensity) {
  float gamma = 1.0 / 2.0;
  float alpha = 12.0f;
  return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)),
                    glm::vec3(0.0), glm::vec3(1.0));
}

int main(int argc, const char *argv[]) {

  clock_t t = clock(); // variable for keeping the time of the rendering

  int width = 1024; // width of the image
  int height = 768; // height of the image
  float fov = 100;  // field of view

  sceneDefinition(); // Let's define a scene

  emit_photons(num_photons, 3.f);

  Image image(width,
              height); // Create an image where we will store the result
  vector<glm::vec3> image_values(width * height);

  float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
  float X = -s * width / 2;
  float Y = s * height / 2;

  for (int i = 0; i < width; i++)
    for (int j = 0; j < height; j++) {

      float dx = X + i * s + s / 2;
      float dy = Y - j * s - s / 2;
      float dz = 1;

      // camera position
      glm::vec3 origin(0, 0, 0);
      glm::vec3 direction(dx, dy, dz);
      direction += origin;
      glm::vec3 pixel_color = glm::vec3(0.0);

      if (super_sampling) {
        float offset = s / 4.0f;
        for (int m = -1; m <= 1; m += 2) {
          for (int n = -1; n <= 1; n += 2) {
            float ddx = dx + m * offset;
            float ddy = dy + n * offset;
            glm::vec3 dir(ddx, ddy, dz);
            dir += origin;
            dir = glm::normalize(dir);
            Ray r(origin, dir);
            pixel_color += trace_ray(r) / 4.0f;
          }
        }
      } else {
        direction = glm::normalize(direction);
        Ray ray(origin, direction);
        pixel_color = trace_ray(ray);
      }
      image.setPixel(i, j, toneMapping(pixel_color));
    }

  t = clock() - t;
  cout << "It took " << ((float)t) / CLOCKS_PER_SEC
       << " seconds to render the image." << endl;
  cout << "I could render at " << (float)CLOCKS_PER_SEC / ((float)t)
       << " frames per second." << endl;

  // Writing the final results of the rendering
  if (argc == 2) {
    image.writeImage(argv[1]);
  } else {
    image.writeImage("./result.ppm");
  }

  return 0;
}
