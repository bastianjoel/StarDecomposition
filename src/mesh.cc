#include "mesh.h"
#include "triray.h"
#include <algorithm>
#include <set>

Vector3q Mesh::face_center(OpenMesh::FaceHandle fh) {
    Vector3q center = Vector3q(0, 0, 0);
    for (auto v : fv_range(fh)) {
        center += data(v).point_q();
    }
    return center / 3;
}

// https://math.stackexchange.com/questions/51326/determining-if-an-arbitrary-point-lies-inside-a-triangle-defined-by-three-points
// https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/point_in_triangle.html
bool Mesh::point_on_face(OpenMesh::FaceHandle fh, Vector3q p) {
    std::vector<Vector3q> t;
    for (auto hv : fv_range(fh)) {
        t.push_back(data(hv).point_q());
    }
    t[0] -= p;
    t[1] -= p;
    t[2] -= p;

    auto u = t[1].cross(t[2]);
    auto v = t[2].cross(t[0]);
    auto w = t[0].cross(t[1]);

    if (u.dot(v) < 0 || u.dot(w) < 0) {
        return false;
    }

    return true;
}

OpenMesh::SmartVertexHandle Mesh::add_vertex_q(const Vector3q p) {
    auto fixVertex = add_vertex(p.unaryExpr([](mpq_class x) { return x.get_d(); }));
    data(fixVertex).set_point_q(p);

    return fixVertex;
}

Vector3q Mesh::update_normal_q(OpenMesh::FaceHandle fh) {
    std::vector<Vector3q> triangle;
    for (auto fv : fv_range(fh)) {
        triangle.push_back(data(fv).point_q());
    }
    Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
    set_normal_q(fh, n);

    return n;
}

void Mesh::set_normal_q(OpenMesh::FaceHandle fh, const Vector3q& n) {
    data(fh).set_normal_q(n);
    set_normal(fh, n.unaryExpr([](mpq_class x) { return x.get_d(); }));
}

void Mesh::move(OpenMesh::VertexHandle vh, const Vector3q& p) {
    data(vh).set_point_q(p);
    set_point(vh, p.unaryExpr([](mpq_class x) { return x.get_d(); }));
    for (auto f : vf_range(vh)) {
        update_normal_q(f);
    }
}

OpenMesh::FaceHandle Mesh::triangle_intersects(const std::vector<Vector3q>& t, const std::vector<OpenMesh::VertexHandle>& borderVertices) {
    Vector3q n = (t[1] - t[0]).cross(t[2] - t[0]);
    if (!_bvh) {
        for (auto face : this->faces()) {
            if (triangle_intersects(t, n, borderVertices, face)) {
                return face;
            }
        }
    } else {
        return triangle_intersects_bvh(_bvh, t, n, borderVertices);
    }

    return OpenMesh::FaceHandle();
}

OpenMesh::FaceHandle Mesh::triangle_intersects_bvh(BVHNode* node, const std::vector<Vector3q>& t, const Vector3q& n, const std::vector<OpenMesh::VertexHandle>& borderVertices) {
    AABB aabb = AABB(t);

    if (!aabb.overlaps(node->aabb)) {
        return OpenMesh::FaceHandle();
    }

    if (node->is_leaf) {
        for (auto f : node->faces) {
            if (triangle_intersects(t, n, borderVertices, f)) {
                return f;
            }
        }
    } else {
        OpenMesh::FaceHandle left = triangle_intersects_bvh(node->l, t, n, borderVertices);
        if (left.is_valid()) {
            return left;
        }

        OpenMesh::FaceHandle right = triangle_intersects_bvh(node->r, t, n, borderVertices);
        if (right.is_valid()) {
            return right;
        }
    }

    return OpenMesh::FaceHandle();
}

bool Mesh::triangle_intersects(const std::vector<Vector3q>& t, const Vector3q& n, const std::vector<OpenMesh::VertexHandle>& borderVertices, const OpenMesh::FaceHandle& face) {
    std::vector<OpenMesh::VertexHandle> fv = {
        this->to_vertex_handle(this->halfedge_handle(face)),
        this->to_vertex_handle(this->next_halfedge_handle(this->halfedge_handle(face))),
        this->to_vertex_handle(this->next_halfedge_handle(this->next_halfedge_handle(this->halfedge_handle(face))))
    };

    // Move connecting vertices slightly to avoid intersection
    std::vector<bool> shareV = { false, false, false };
    for (auto v : borderVertices) {
        for (int i = 0; i < 3; i++) {
            if (v == fv[i]) {
                shareV[i] = true;
            }
        }
    }
    short sharedVertices = shareV[0] + shareV[1] + shareV[2];
    ASSERT(sharedVertices < 3);

    Vector3q fNormal = this->data(face).normal_q();
    std::vector<Vector3q> vq = { this->data(fv[0]).point_q(), this->data(fv[1]).point_q(), this->data(fv[2]).point_q() };
    if (shareV[0] && shareV[1]) {
        vq[0] = vq[0] + (vq[2] - vq[0]) * 1e-6;
        vq[1] = vq[1] + (vq[2] - vq[1]) * 1e-6;
    } else if (shareV[0] && shareV[2]) {
        vq[0] = vq[0] + (vq[1] - vq[0]) * 1e-6;
        vq[2] = vq[2] + (vq[1] - vq[2]) * 1e-6;
    } else if (shareV[1] && shareV[2]) {
        vq[1] = vq[1] + (vq[0] - vq[1]) * 1e-6;
        vq[2] = vq[2] + (vq[0] - vq[2]) * 1e-6;
    } else if (shareV[0]) {
        vq[0] = vq[0] + (vq[2] - vq[0]) * 1e-6 + (vq[1] - vq[0]) * 1e-6;
    } else if (shareV[1]) {
        vq[1] = vq[1] + (vq[2] - vq[1]) * 1e-6 + (vq[0] - vq[1]) * 1e-6;
    } else if (shareV[2]) {
        vq[2] = vq[2] + (vq[1] - vq[2]) * 1e-6 + (vq[0] - vq[2]) * 1e-6;
    }

    return tri_tri_intersect(t[0], t[1], t[2], n, vq[0], vq[1], vq[2], fNormal);
}

OpenMesh::FaceHandle Mesh::ray_intersects(const Vector3q& o, const Vector3q& n) {
    mpq_class distResult;
    return ray_intersects(o, n, distResult);
}

OpenMesh::FaceHandle Mesh::ray_intersects(const Vector3q& o, const Vector3q& n, mpq_class& distResult) {
    distResult = -1;
    if (!_bvh) {
        OpenMesh::FaceHandle opposite;
        for (auto face : this->faces()) {
            std::vector<Vector3q> t;
            for (auto v : fv_range(face)) {
                t.push_back(data(v).point_q());
            }

            mpq_class dist;
            bool res = tri_ray_intersect(o, n, t[0], t[1], t[2], dist);
            if (!res) {
                continue;
            }

            if (distResult == -1 || dist < distResult) {
                distResult = dist;
                opposite = face;
            }
        }

        return opposite;
    } else {
        auto intersection = ray_intersects_bvh(_bvh, o, n);
        distResult = intersection.second;

        return intersection.first;
    }

    return OpenMesh::FaceHandle();
}

std::pair<OpenMesh::FaceHandle, mpq_class> Mesh::ray_intersects_bvh(BVHNode* node, const Vector3q& o, const Vector3q& n) {
    if (!node->aabb.ray_intersects(o, n)) {
        return std::make_pair(OpenMesh::FaceHandle(), -1);
    }

    if (node->is_leaf) {
        mpq_class distResult = -1;
        OpenMesh::FaceHandle opposite;
        for (auto f : node->faces) {
            std::vector<Vector3q> t;
            for (auto v : fv_range(f)) {
                t.push_back(data(v).point_q());
            }

            mpq_class dist;
            if (tri_ray_intersect(o, n, t[0], t[1], t[2], dist) && (distResult == -1 || dist < distResult)) {
                distResult = dist;
                opposite = f;
            }
        }

        return std::make_pair(opposite, distResult);
    } else {
        mpq_class leftDist, rightDist, rayExit;
        node->l->aabb.ray_intersects(o, n, leftDist, rayExit);
        node->r->aabb.ray_intersects(o, n, rightDist, rayExit);

        if (leftDist < rightDist) {
            auto left = ray_intersects_bvh(node->l, o, n);
            if (left.first.is_valid()) {
                return left;
            }

            auto right = ray_intersects_bvh(node->r, o, n);
            if (right.first.is_valid()) {
                return right;
            }
        } else {
            auto right = ray_intersects_bvh(node->r, o, n);
            if (right.first.is_valid()) {
                return right;
            }

            auto left = ray_intersects_bvh(node->l, o, n);
            if (left.first.is_valid()) {
                return left;
            }
        }
    }

    return std::make_pair(OpenMesh::FaceHandle(), -1);
}

std::pair<bool, Vector3q> Mesh::star_center() {
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> normals;
    for (auto face : this->faces()) {
        if (face.is_valid() && !face.deleted()) {
            positions.push_back(this->point(this->from_vertex_handle(this->halfedge_handle(face))));
            normals.push_back(this->normal(face).normalized());
        }
    }

    Vector3q newCenter = kernel_chebyshev_center(positions, normals).cast<mpq_class>();
    for (auto face : this->faces()) {
        std::vector<Vector3q> triangle;
        for (auto fv : this->fv_range(face)) {
            triangle.push_back(this->data(fv).point_q());
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        if ((triangle[0] - newCenter).dot(n) <= 0) {
            return std::make_pair(false, Vector3q());
        }
    }

    return std::make_pair(true, newCenter);
}

std::pair<bool, mpq_class> Mesh::intersection_factor(const Vector3q& p, const Vector3q& dir, OpenMesh::FaceHandle face) {
    auto faceNormal = this->data(face).normal_q();
    auto div = faceNormal.transpose() * dir;
    if (div[0] == 0) {
        return std::make_pair(false, mpq_class());
    }

    // Calculate intersection
    auto c = faceNormal.transpose() * this->data(this->to_vertex_handle(this->halfedge_handle(face))).point_q();
    auto nominator = faceNormal.transpose() * p - c;
    auto denominator = faceNormal.transpose() * dir;
    return std::make_pair(true, (-nominator / denominator)[0]);
}

std::vector<OpenMesh::HalfedgeHandle> Mesh::boundary_halfedges() {
    std::vector<OpenMesh::HalfedgeHandle> boundary_halfedges = {};
    for (auto h : halfedges()) {
        if (is_boundary(h)) {
            boundary_halfedges.push_back(h);
        }
    }
    /*
    if (!lastBoundaryHalfedge.is_valid() || !is_boundary(lastBoundaryHalfedge)) {
        lastBoundaryHalfedge = OpenMesh::HalfedgeHandle();
        for (auto h : halfedges()) {
            if (is_boundary(h)) {
                lastBoundaryHalfedge = h;
            }
        }

        if (!lastBoundaryHalfedge.is_valid() || !is_boundary(lastBoundaryHalfedge)) {
            return std::vector<OpenMesh::HalfedgeHandle>();
        }
    }

    std::vector<OpenMesh::HalfedgeHandle> boundary_halfedges = { lastBoundaryHalfedge };
    auto heh = next_halfedge_handle(lastBoundaryHalfedge);
    while(heh != lastBoundaryHalfedge && heh.is_valid()) {
        boundary_halfedges.push_back(heh);
        heh = next_halfedge_handle(heh);
    }
    */

    return boundary_halfedges;
}

void Mesh::generate_bvh() {
    if (_bvh != nullptr) {
        delete _bvh;
    }

    std::vector<OpenMesh::FaceHandle> faces;
    for (auto f : this->faces()) {
        faces.push_back(f);
    }
    _bvh = generate_bvh(faces);
}

BVHNode* Mesh::generate_bvh(std::vector<OpenMesh::FaceHandle>& faces, int depth) {
    std::set<OpenMesh::VertexHandle> vertices;
    for (auto f : faces) {
        for (auto v : fv_range(f)) {
            vertices.insert(v);
        }
    }

    std::vector<Vector3q> positions;
    for (auto v : vertices) {
        positions.push_back(data(v).point_q());
    }
    AABB aabb(positions);
    BVHNode* node = new BVHNode(aabb, faces);

    // Stop criterion: Maximum depth or few triangles
    if (faces.size() <= 5  || depth > 20) {
        return node;
    }

    node->is_leaf = false;

    // Choose split axis
    int split_axis = 0;
    mpq_class max_extent = node->aabb.max[0] - node->aabb.min[0];
    for (int axis = 1; axis < 3; ++axis) {
        mpq_class extent = node->aabb.max[axis] - node->aabb.min[axis];
        if (extent > max_extent) {
            max_extent = extent;
            split_axis = axis;
        }
    }

    // Split faces along median
    std::sort(faces.begin(), faces.end(), [&](const OpenMesh::FaceHandle& a, const OpenMesh::FaceHandle& b) {
        mpq_class center_a = face_center(a)[split_axis];
        mpq_class center_b = face_center(b)[split_axis];

        return center_a < center_b;
    });
    size_t median_index = faces.size() / 2;
    std::vector<OpenMesh::FaceHandle> left_faces(faces.begin(), faces.begin() + median_index);
    std::vector<OpenMesh::FaceHandle> right_faces(faces.begin() + median_index, faces.end());

    node->l = generate_bvh(left_faces, depth + 1);
    node->r = generate_bvh(right_faces, depth + 1);

    return node;
}
