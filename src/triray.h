/* 
 * triray.h
 *
 * Triangle-ray intersection test based on Möller–Trumbore intersection algorithm
 * https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
 */

#pragma once

#include "vectorq.h"


inline bool tri_ray_intersect(const Vector3q& ray_origin, const Vector3q& ray_vector, const Vector3q& V0, const Vector3q& V1, const Vector3q& V2, mpq_class& factor)
{
    Vector3q edge1 = V1 - V0;
    Vector3q edge2 = V2 - V0;
    Vector3q ray_cross_e2 = ray_vector.cross(edge2);
    mpq_class det = edge1.dot(ray_cross_e2);

    if (det == 0)
        return false;    // This ray is parallel to this triangle.

    mpq_class inv_det = 1.0f / det;
    Vector3q s = ray_origin - V0;
    mpq_class u = inv_det * s.dot(ray_cross_e2);

    if (u < 0.0f || u > 1.0f)
        return false;

    Vector3q s_cross_e1 = s.cross(edge1);
    mpq_class v = inv_det * ray_vector.dot(s_cross_e1);

    if (v < 0.0f || u + v > 1.0f)
        return false;

    // At this stage we can compute t to find out where the intersection point is on the line.
    factor = inv_det * edge2.dot(s_cross_e1);

    if (factor > 0) // ray intersection
    {
        // return  Vector3q(ray_origin + ray_vector * t);
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}

inline bool tri_ray_intersect(const Vector3q& ray_origin, const Vector3q& ray_vector, const Vector3q& V0, const Vector3q& V1, const Vector3q& V2)
{
    mpq_class factor;
    return tri_ray_intersect(ray_origin, ray_vector, V0, V1, V2, factor);
}
