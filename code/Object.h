#ifndef Object_h
#define Object_h

using namespace std;
#include "Material.h"
#include "Ray.h"
#include "glm/glm.hpp"

/**
 General class for the object
 */
class Object {

protected:
  glm::mat4
      transformationMatrix; ///< Matrix representing the transformation from the
                            ///< local to the global coordinate system
  glm::mat4 inverseTransformationMatrix; ///< Matrix representing the
                                         ///< transformation from the global to
                                         ///< the local coordinate system
  glm::mat4 normalMatrix; ///< Matrix for transforming normal vectors from the
                          ///< local to the global coordinate system

public:
  glm::vec3 color;   ///< Color of the object
  Material material; ///< Structure describing the material of the object
  /** A function computing an intersection, which returns the structure Hit */
  virtual Hit intersect(Ray ray) = 0;

  /** Function that returns the material struct of the object*/
  Material getMaterial() { return material; }
  /** Function that set the material
   @param material A structure describing the material of the object
  */
  void setMaterial(Material material) { this->material = material; }
  /** Functions for setting up all the transformation matrices
  @param matrix The matrix representing the transformation of the object in the
  global coordinates */
  void setTransformation(glm::mat4 matrix) {

    transformationMatrix = matrix;

    inverseTransformationMatrix = glm::inverse(matrix);
    normalMatrix = glm::transpose(inverseTransformationMatrix);
  }
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object {
private:
  float radius;     ///< Radius of the sphere
  glm::vec3 center; ///< Center of the sphere

public:
  /**
   The constructor of the sphere
   @param radius Radius of the sphere
   @param center Center of the sphere
   @param color Color of the sphere
   */
  Sphere(float radius, glm::vec3 center, glm::vec3 color)
      : radius(radius), center(center) {
    this->color = color;
  }
  Sphere(float radius, glm::vec3 center, Material material)
      : radius(radius), center(center) {
    this->material = material;
  }
  /** Implementation of the intersection function*/
  Hit intersect(Ray ray) {

    glm::vec3 c = center - ray.origin;

    float cdotc = glm::dot(c, c);
    float cdotd = glm::dot(c, ray.direction);

    Hit hit;

    float D = 0;
    if (cdotc > cdotd * cdotd) {
      D = sqrt(cdotc - cdotd * cdotd);
    }
    if (D <= radius) {
      hit.hit = true;
      float t1 = cdotd - sqrt(radius * radius - D * D);
      float t2 = cdotd + sqrt(radius * radius - D * D);

      float t = t1;
      if (t < 0)
        t = t2;
      if (t < 0) {
        hit.hit = false;
        return hit;
      }

      hit.intersection = ray.origin + t * ray.direction;
      hit.normal = glm::normalize(hit.intersection - center);
      hit.distance = glm::distance(ray.origin, hit.intersection);
      hit.object = this;
    } else {
      hit.hit = false;
    }
    return hit;
  }
};

class Plane : public Object {

private:
  glm::vec3 normal;
  glm::vec3 point;

public:
  Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal) {}
  Plane(glm::vec3 point, glm::vec3 normal, Material material)
      : point(point), normal(normal) {
    this->material = material;
  }
  Hit intersect(Ray ray) {

    Hit hit;
    hit.hit = false;

    float DdotN = glm::dot(ray.direction, normal);
    if (DdotN < 0) {

      float PdotN = glm::dot(point - ray.origin, normal);
      float t = PdotN / DdotN;

      if (t > 0) {
        hit.hit = true;
        hit.normal = normal;
        hit.distance = t;
        hit.object = this;
        hit.intersection = t * ray.direction + ray.origin;
      }
    }

    return hit;
  }
};

class Cone : public Object {
private:
  Plane *plane;

public:
  Cone(Material material) {
    this->material = material;
    plane = new Plane(glm::vec3(0, 1, 0), glm::vec3(0.0, 1, 0));
  }
  Hit intersect(Ray ray) {

    Hit hit;
    hit.hit = false;

    glm::vec3 d = inverseTransformationMatrix *
                  glm::vec4(ray.direction, 0.0); // implicit cast to vec3
    glm::vec3 o = inverseTransformationMatrix *
                  glm::vec4(ray.origin, 1.0); // implicit cast to vec3
    d = glm::normalize(d);

    float a = d.x * d.x + d.z * d.z - d.y * d.y;
    float b = 2 * (d.x * o.x + d.z * o.z - d.y * o.y);
    float c = o.x * o.x + o.z * o.z - o.y * o.y;

    float delta = b * b - 4 * a * c;

    if (delta < 0) {
      return hit;
    }

    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);

    float t = t1;
    hit.intersection = o + t * d;
    if (t < 0 || hit.intersection.y > 1 || hit.intersection.y < 0) {
      t = t2;
      hit.intersection = o + t * d;
      if (t < 0 || hit.intersection.y > 1 || hit.intersection.y < 0) {
        return hit;
      }
    };

    hit.normal =
        glm::vec3(hit.intersection.x, -hit.intersection.y, hit.intersection.z);
    hit.normal = glm::normalize(hit.normal);

    Ray new_ray(o, d);
    Hit hit_plane = plane->intersect(new_ray);
    if (hit_plane.hit && hit_plane.distance < t &&
        length(hit_plane.intersection - glm::vec3(0, 1, 0)) <= 1.0) {
      hit.intersection = hit_plane.intersection;
      hit.normal = hit_plane.normal;
    }

    hit.hit = true;
    hit.object = this;
    hit.intersection =
        transformationMatrix *
        glm::vec4(hit.intersection, 1.0); // implicit cast to vec3
    hit.normal =
        (normalMatrix * glm::vec4(hit.normal, 0.0)); // implicit cast to vec3
    hit.normal = glm::normalize(hit.normal);
    hit.distance = glm::length(hit.intersection - ray.origin);

    return hit;
  }
};

class Square : public Object {
private:
  Plane *plane;
  glm::vec3 position;
  glm::vec3 u_vec;
  glm::vec3 v_vec;
  glm::vec3 normal;
  glm::vec3 corner; // Lower-left corner
  float u_length_sq;
  float v_length_sq;

public:
  Square(glm::vec3 position, glm::vec3 u_vec, glm::vec3 v_vec, Material material)
      : position(position), u_vec(u_vec), v_vec(v_vec) {
    this->material = material;
    normal = -glm::normalize(glm::cross(u_vec, v_vec));

    corner = position - 0.5f * (u_vec + v_vec);

    plane = new Plane(corner, normal);

    u_length_sq = glm::dot(u_vec, u_vec);
    v_length_sq = glm::dot(v_vec, v_vec);
  }

  Hit intersect(Ray ray) {
    Hit hit;
    hit.hit = false;
    
    Hit plane_hit = plane->intersect(ray);
    if (!plane_hit.hit) {
      return hit;
    }

    glm::vec3 vec_to_point = plane_hit.intersection - corner;
    float alpha = glm::dot(vec_to_point, u_vec) / u_length_sq;
    float beta = glm::dot(vec_to_point, v_vec) / v_length_sq;

    if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1) {
      return hit;
    }

    hit.object = this;
    hit.hit = true;
    hit.intersection = plane_hit.intersection;
    hit.normal = plane_hit.normal;
    hit.distance = plane_hit.distance;

    return hit;
  }

  glm::vec3 getNormal() { return normal; }
};

#endif /* Object_h */