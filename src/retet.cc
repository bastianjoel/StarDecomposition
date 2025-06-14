#include "retet.h"
#include "OpenMesh/Core/Mesh/Handles.hh"

bool retetrahedrize(VolumeMesh& mesh) {
    tetgenbehavior behavior;
    behavior.plc = 1;
    behavior.quality = 1;
    behavior.nobisect = 1;
    behavior.quiet = 1;
    tetgenio in;
    in.firstnumber = 0;
    std::vector<Vertex> boundary_vertices = mesh.boundary_vertices();
    auto idx = mesh.property<Vertex, int>();
    in.numberofpoints = boundary_vertices.size();
    in.pointlist = new double[in.numberofpoints * 3];
    for (int i = 0; i < in.numberofpoints; i++) {
        Vertex v = boundary_vertices[i];
        idx[v] = i;
        Eigen::Vector3d p = mesh.position(v);
        for (int j = 0; j < 3; j++) {
            in.pointlist[3 * i + j] = p(j);
        }
    }
    std::vector<Halfface> boundary_halffaces = mesh.boundary_halffaces();
    in.numberoffacets = boundary_halffaces.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    for (int i = 0; i < in.numberoffacets; i++) {
        Halfface h = boundary_halffaces[i];
        tetgenio::polygon polygon;
        polygon.numberofvertices = 3;
        polygon.vertexlist = new int[polygon.numberofvertices];
        int j = 0;
        for (auto hv : mesh.halfface_vertices(h)) {
            polygon.vertexlist[j] = idx[hv];
            j++;
        }
        tetgenio::facet facet;
        facet.numberofpolygons = 1;
        facet.polygonlist = new tetgenio::polygon[facet.numberofpolygons];
        facet.polygonlist[0] = polygon;
        facet.numberofholes = 0;
        facet.holelist = nullptr;
        in.facetlist[i] = facet;
    }
    tetgenio out;
    try {
        tetrahedralize(&behavior, &in, &out);
    } catch (int error_code) {
        std::cout << "TetGen failed, error code " << error_code << std::endl;
        return false;
    } catch (...) {
        std::cout << "TetGen failed" << std::endl;
        return false;
    }
    mesh = VolumeMesh();
    std::vector<Vertex> vertices;
    for (int i = 0; i < out.numberofpoints; i++) {
        Eigen::Vector3d p;
        for (int j = 0; j < 3; j++) {
            p(j) = out.pointlist[3 * i + j];
        }
        Vertex v = mesh.add_vertex(p);
        vertices.push_back(v);
    }
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        std::vector<Vertex> t;
        for (int j = 0; j < 4; j++) {
            t.push_back(vertices[out.tetrahedronlist[4 * i + j]]);
        }
        mesh.add_cell(t, true);
    }
    for (auto c : mesh.cells()) {
        if (mesh.degenerate_or_inverted(c)) {
            std::cout << "TetGen generated degenerate or inverted tetrahedron" << std::endl;
            return false;
        }
    }
    return true;
}

bool retetrahedrize(const Mesh& mesh, VolumeMesh& outMesh) {
    tetgenbehavior behavior;
    behavior.plc = 1;
    behavior.quality = 1;
    behavior.nobisect = 1;
    behavior.quiet = 1;
    tetgenio in;
    in.firstnumber = 0;
    in.numberofpoints = mesh.n_vertices();
    in.pointlist = new double[in.numberofpoints * 3];
    in.pointmarkerlist = new int[in.numberofpoints];
    for (auto v : mesh.vertices()) {
        Eigen::Vector3d p = mesh.point(v);
        for (int j = 0; j < 3; j++) {
            in.pointlist[3 * v.idx() + j] = p(j);
        }
        in.pointmarkerlist[v.idx()] = v.idx() + 1;
    }
    in.numberoffacets = mesh.n_faces();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    for (auto f : mesh.faces()) {
        tetgenio::polygon polygon;
        polygon.numberofvertices = 3;
        polygon.vertexlist = new int[polygon.numberofvertices];
        int j = 0;
        for (auto fv : mesh.fv_range(f)) {
            polygon.vertexlist[j] = fv.idx();
            j++;
        }
        tetgenio::facet facet;
        facet.numberofpolygons = 1;
        facet.polygonlist = new tetgenio::polygon[facet.numberofpolygons];
        facet.polygonlist[0] = polygon;
        facet.numberofholes = 0;
        facet.holelist = nullptr;
        in.facetlist[f.idx()] = facet;
    }
    tetgenio out;
    try {
        tetrahedralize(&behavior, &in, &out);
    } catch (int error_code) {
        std::cout << "TetGen failed, error code " << error_code << std::endl;
        return false;
    } catch (...) {
        std::cout << "TetGen failed" << std::endl;
        return false;
    }
    outMesh = VolumeMesh();
    auto Q = outMesh.property<Vertex, Vector3q>("Q");
    std::vector<Vertex> vertices;
    for (int i = 0; i < out.numberofpoints; i++) {
        Eigen::Vector3d p;
        for (int j = 0; j < 3; j++) {
            p(j) = out.pointlist[3 * i + j];
        }
        Vertex v = outMesh.add_vertex(p);
        if (out.pointmarkerlist[i] > 0) {
            int vIdx = out.pointmarkerlist[i] - 1;
            auto oldV = mesh.vertex_handle(vIdx);
            Q[v] = mesh.data(oldV).point_q();
        } else {
            Q[v] = p.cast<mpq_class>();
        }
        vertices.push_back(v);
    }
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        std::vector<Vertex> t;
        for (int j = 0; j < 4; j++) {
            t.push_back(vertices[out.tetrahedronlist[4 * i + j]]);
        }
        outMesh.add_cell(t, true);
    }
    for (auto c : outMesh.cells()) {
        if (outMesh.degenerate_or_inverted(c)) {
            std::cout << "TetGen generated degenerate or inverted tetrahedron" << std::endl;
            return false;
        }
    }
    return true;
}
