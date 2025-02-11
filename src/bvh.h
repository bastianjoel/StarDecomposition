#pragma once

#include "OpenMesh/Core/Mesh/Handles.hh"
#include "vectorq.h"
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <vector>

// Axis-aligned min bounding box
struct AABB
{
    Vector3q min;
    Vector3q max;

    AABB(const std::vector<Vector3q>& positions) {
        min = max = positions[0];
        for (auto& p : positions) {
            min = min.cwiseMin(p);
            max = max.cwiseMax(p);
        }
    }

    AABB(const AABB& l, const AABB& r) {
        min = l.min.cwiseMin(r.min);
        max = l.max.cwiseMax(r.max);
    }

    bool overlaps(const AABB& other) {
        return (min[0] <= other.max[0] && max[0] >= other.min[0]) &&
               (min[1] <= other.max[1] && max[1] >= other.min[1]) &&
               (min[2] <= other.max[2] && max[2] >= other.min[2]);
    }
};

// Bounding Volume Hierarchy Node
struct BVHNode
{
    BVHNode* l;
    BVHNode* r;
    AABB aabb;
    bool is_leaf;
    std::vector<OpenMesh::FaceHandle> faces;

    BVHNode(const AABB& aabb, std::vector<OpenMesh::FaceHandle>& faces) : is_leaf(true), aabb(aabb), faces(faces), l(nullptr), r(nullptr) {}

    ~BVHNode() {
        if (!is_leaf) {
            delete l;
            delete r;
            is_leaf = true;
        }
    }
};
