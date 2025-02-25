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

    bool ray_intersects(const Vector3q& origin, const Vector3q& direction) {
        mpq_class t_entry;
        mpq_class t_exit;

        return ray_intersects(origin, direction, t_entry, t_exit);
    }

    // Ray-AABB intersection using the slab method
    bool ray_intersects(const Vector3q& origin, const Vector3q& direction, mpq_class& t_entry, mpq_class& t_exit) {
        int valid = 0;
        for (int i = 0; i < 3; i++) {
            if (direction[i] != 0) {
                valid = i;
            }
        }

        Vector3q t0, t1;
        for (int i = 0; i < 3; i++) {
            if (direction[i] != 0) {
                t0[i] = (min[i] - origin[i]) / direction[i];
                t1[i] = (max[i] - origin[i]) / direction[i];
            } else {
                t0[valid] = (min[valid] - origin[valid]) / direction[valid];
                t1[valid] = (max[valid] - origin[valid]) / direction[valid];
            }
        }

        Vector3q tClose = t0.cwiseMin(t1);
        Vector3q tFar = t0.cwiseMax(t1);

        t_entry = tClose.maxCoeff();
        t_exit = tFar.minCoeff();

        return t_exit >= t_entry;
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
