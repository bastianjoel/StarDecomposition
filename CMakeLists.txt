cmake_minimum_required(VERSION 3.1)
project(StarDecomposition)
set(CMAKE_CXX_STANDARD 17)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

option(GUI "GUI" ON)
option(SAVE_DEBUG_MESHES "SAVE_DEBUG_MESHES" ON)

# Eigen
find_package(Eigen3 REQUIRED CONFIG)

# CGAL
find_package(CGAL REQUIRED)

# GMP
find_package(GMP REQUIRED)
find_package(GMPXX REQUIRED)

# OpenVolumeMesh
add_subdirectory(ext/OpenVolumeMesh)

# OpenMesh
add_subdirectory(ext/OpenMesh)

# Exact predicates
add_library(predicates)
target_sources(predicates PRIVATE ext/predicates/predicates.c)
target_include_directories(predicates PUBLIC ext/predicates)

# TetGen
add_subdirectory(ext/tetgen)

if(GUI)
    # Dear ImGui
    add_library(glad)
    target_sources(glad PRIVATE ext/glad/src/glad.c)
    target_include_directories(glad PUBLIC ext/glad/include)
    add_subdirectory(ext/glfw)
    add_library(imgui)
    set(IMGUI_SOURCES
        ext/imgui/imgui.cpp
        ext/imgui/imgui_demo.cpp
        ext/imgui/imgui_draw.cpp
        ext/imgui/imgui_widgets.cpp
        ext/imgui/imgui_tables.cpp
        ext/imgui/backends/imgui_impl_glfw.cpp
        ext/imgui/backends/imgui_impl_opengl3.cpp
    )
    target_sources(imgui PRIVATE ${IMGUI_SOURCES})
    target_link_libraries(imgui PUBLIC glad glfw)
    target_include_directories(imgui PUBLIC ext/imgui ext/imgui/backends)

    # stb
    add_library(stb INTERFACE)
    target_include_directories(stb INTERFACE ext/stb)

    # tinyfiledialogs
    add_library(tinyfiledialogs)
    target_sources(tinyfiledialogs PRIVATE ext/tinyfiledialogs/tinyfiledialogs.c)
    target_include_directories(tinyfiledialogs PUBLIC ext/tinyfiledialogs)
endif()

add_library(StarDecomposition)
target_sources(StarDecomposition PRIVATE src/decompose.cc src/sd_boundary_lp.cc src/sd_boundary_chebyshev.cc src/sd_boundary.cc src/lp.cc src/main.cc src/mesh.cc src/volume_mesh.cc src/retet.cc src/sd.cc)
target_link_libraries(StarDecomposition Eigen3::Eigen CGAL::CGAL ${GMP_LIBRARIES} ${GMPXX_LIBRARIES} OpenVolumeMesh::OpenVolumeMesh OpenMeshCore predicates tet)
target_include_directories(StarDecomposition PUBLIC ${CGAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GMPXX_INCLUDE_DIR} ext/OpenVolumeMesh/src ${CMAKE_BINARY_DIR}/ext/OpenVolumeMesh/src ext/OpenMesh/src ext/tetgen src)

add_executable(sd)
target_sources(sd PRIVATE src/main.cc)
target_link_libraries(sd StarDecomposition)

if(GUI)
    target_sources(sd PRIVATE src/viewer.cc)
    target_link_libraries(sd imgui tinyfiledialogs stb)
    target_compile_definitions(sd PRIVATE GUI SAVE_DEBUG_MESHES)

    target_sources(StarDecomposition PRIVATE src/viewer.cc)
    target_link_libraries(StarDecomposition imgui tinyfiledialogs stb)
    target_compile_definitions(StarDecomposition PRIVATE GUI SAVE_DEBUG_MESHES)
endif()
