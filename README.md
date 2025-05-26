# Star Decomposition

## Prerequisites

The code can be used as either a **library** or an **executable**.

In either case, the following libraries need to be installed:

- [Eigen](https://eigen.tuxfamily.org)
- [CGAL](https://www.cgal.org)
- [GMP](https://gmplib.org)

The following libraries are included in the [ext](ext) directory:

- [OpenVolumeMesh](https://gitlab.vci.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git)
- [OpenMesh](https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git)
- [Exact predicates](https://www.cs.cmu.edu/~quake/robust.html)
- [TetGen](http://www.tetgen.org)

For the viewer:

- [Dear ImGui](https://github.com/ocornut/imgui.git)
  - [Glad](https://glad.dav1d.de)
  - [GLFW](https://github.com/glfw/glfw.git)
- [stb](https://github.com/nothings/stb)
- [tinyfiledialogs](https://sourceforge.net/projects/tinyfiledialogs/)

Remember to clone this repository recursively:\
```git clone --recursive https://github.com/bastianjoel/StarDecomposition.git```

## Library

### CMake

The library can be included via CMake like this:

```cmake
add_subdirectory(path/to/StarDecomposition)
target_link_libraries(yourTarget StarDecomposition)
```

### Usage

The function ```sd``` implements the method of the paper.
Input is an arbirary mesh $M$ represented using
[OpenVolumeMesh](https://www.graphics.rwth-aachen.de/software/openvolumemesh/) or
[OpenMesh](https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git).
In case of a volume mesh it is required to have no degenerate or inverted tetrahedra.
The function ```retetrahedrize``` using
[TetGen](http://www.tetgen.org)
can be used to achieve this.
Output are the star-shaped meshes $M_i$.
Their vertices are available in rational numbers in the property ```Q``` of type ```Vector3q```.

```cpp
#include <sd.h>

Mesh M = ...;

Mesh components = sd(M);
```

## Executable

### Building

- ```mkdir build```
- ```cd build```
- ```cmake [-SAVE_DEBUG_MESHES=OFF] [-DGUI=OFF] ..```
- ```make -j```

The option ```GUI``` controls if the program is built with the viewer.   
The option ```SAVE_DEBUG_MESHES``` controls if the program should output debug meshes to the directory ```debug```.


### Usage

```./sd <M> [-o <out_dir>] [-f] [-b] [-s <number>] [-a tet|boundary-lp|boundary]```

- ```<M>```:
Tetrahedral mesh $M$.
- ```-o <out_dir>```:
Optional argument to output the mesh pairs into the directory ```out_dir```.
Meshes are output in ```.ovm``` format
and have the property ```Q_string``` of type ```std::string```
containing the rational vertices as strings.
- ```-f```:
If set, the program will run a feasibility test on the resulting meshes.
The test checks if the component meshes are star-shaped.
- ```-b```:
If set, the program will run in benchmark mode. 
The program will not output any data except for a benchmark result.
The result is in the form of a CSV file with the following columns: 
  - filename: The name of the input file.
  - algorithm: The algorithm used for the star decomposition.
  - time: The time taken for the star decomposition.
  - result_feasible: Whether the star decomposition was feasible or not.
  - result_components: The number of components in the star decomposition.
  - seed: The seed used for the random number generator.
  - result_cell_counts: The number of tetrahedra in each component. Separated by ```|```.
  - result_boundary_face_counts: The number of boundary faces in each component. Separated by ```|```.
- ```-s```:
If set, the program will use the given seed for the random number generator.
- ```-a tet|boundary-lp|boundary```:
The algorithm used for the star decomposition.
  - ```tet```:
  The algorithm uses the decomposition algorithm from https://github.com/SteffenHinderink/StarDecompositionMaps
  - ```boundary-lp```:
  The algorithm variant that uses a linear program to determine the close vertex position.
  - ```boundary```:
  The algorithm variant that dynamically moves the close vertex position.

e.g. ```./sd ../meshes/open.vtk```

The viewer shows which faces are currently added to the component.
Faces that were added are colored in green, faces that could not be added at the 
current stage are colored in red and faces that are impossible to be added to 
the current star component are colored in blue.
When finished the components in different colors are displayed.

## Benchmark
The benchmark results are stored in the ```benchmarks``` directory.
To run benchmarks the `run-bench.sh` script can be used. 
It accepts two arguments:
- The first argument is the path to the directory containing the meshes.
- The second argument is a file name pattern without the extension (e.g. `.vtk`) to filter the meshes.

e.g. `./run-bench.sh ../meshes "*"`

Files of the format `*.vtk`, `*.off` and `*.stl` will be considered.
The benchmark results are stored in the file `results.csv`. 
Note that the `results.csv` file will be overwritten if it already exists.

To get results while the script is still running, you can use the outputs in `.output`:
```bash
cat .output/*.txt | less
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Note that any resources under `ext/` may have their own license, which may or may not be different from the MIT license.
