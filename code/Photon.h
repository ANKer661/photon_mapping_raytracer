#ifndef Photon_h
#define Photon_h

#include "glm/glm.hpp"
#include <vector>

struct Photon {
  glm::vec3 position;
  glm::vec3 color;
  glm::vec3 direction;
};

enum class Axis : uint8_t { X = 0, Y = 1, Z = 2 };

struct KDNode {
  Photon photon;
  int left_child = -1; // -1 for no child
  int right_child = -1;
  Axis axis; // split along which axis
};

class PhotonMap {
private:
  std::vector<Photon> photons;
  std::vector<KDNode> tree;
  int root_index = -1;

  static bool compare_photon(const Photon &a, const Photon &b, Axis axis) {
    switch (axis) {
    case Axis::X:
      return a.position.x < b.position.x;
    case Axis::Y:
      return a.position.y < b.position.y;
    case Axis::Z:
      return a.position.z < b.position.z;
    }
    return false; // should not reach here
  }

  // recursively build kd-tree
  int build_KDTree(int start, int end, int depth) {
    // use [start, end)
    if (start >= end) {
      return -1;
    }

    Axis axis = static_cast<Axis>(depth % 3);
    int mid = start + (end - start) / 2;

    // find median
    std::nth_element(photons.begin() + start, photons.begin() + mid,
                     photons.begin() + end,
                     [axis](const Photon &a, const Photon &b) {
                       return compare_photon(a, b, axis);
                     });

    KDNode node;
    node.photon = photons[mid];
    node.axis = axis;

    int current_index = tree.size();
    tree.push_back(node);

    tree[current_index].left_child = build_KDTree(start, mid, depth + 1);
    tree[current_index].right_child = build_KDTree(mid + 1, end, depth + 1);

    return current_index;
  }

  float distance_sq(const glm::vec3 &a, const glm::vec3 &b) {
    glm::vec3 diff = a - b;
    return glm::dot(diff, diff);
  }

  void search_KDTree(int node_index, const glm::vec3 &position, float radius_sq,
                     std::vector<Photon> &result) {
    if (node_index == -1) {
      return;
    }

    const KDNode &node = tree[node_index];
    float dist_sq = distance_sq(position, node.photon.position);
    // check current photon
    if (dist_sq <= radius_sq) {
      result.push_back(node.photon);
    }

    // decide which side to search
    float diff;
    switch (node.axis) {
    case Axis::X:
      diff = position.x - node.photon.position.x;
      break;
    case Axis::Y:
      diff = position.y - node.photon.position.y;
      break;
    case Axis::Z:
      diff = position.z - node.photon.position.z;
      break;
    default:
      diff = 0.0f; // this will not happen
      break;
    }
    if (diff < 0) {
      search_KDTree(node.left_child, position, radius_sq, result);
      if (diff * diff <= radius_sq) { // only search other side if intersect
        search_KDTree(node.right_child, position, radius_sq, result);
      }
    } else {
      search_KDTree(node.right_child, position, radius_sq, result);
      if (diff * diff <= radius_sq) {
        search_KDTree(node.left_child, position, radius_sq, result);
      }
    }
  }

public:
  void store(Photon &photon) { photons.push_back(photon); }

  void build() {
    if (photons.empty())
      return;

    tree.reserve(photons.size());
    root_index = build_KDTree(0, photons.size(), 0);
  }

  /** Function to find photons within a radius around a position
  @param position The center position to search around
  @param radius The search radius
  @param result Vector to store the result photons
  */
  void find_photons(const glm::vec3 &position, float radius,
                    std::vector<Photon> &result) {
    if (root_index == -1) {
      return;
    }

    float radius_sq = radius * radius;
    search_KDTree(root_index, position, radius_sq, result);
  }

  void find_knearest_photons(const glm::vec3 &position, int k,
                             std::vector<Photon> &result) {
    // todo
  }
};

#endif