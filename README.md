# Diskify: a Combinatorial Approach to the Computation of Cut Networks on Surfaces
This repository contains the implementation of the methodologies described in the paper *Diskify: a Combinatorial Approach to the Computation of Cut Networks on Surfaces*
- DOI: TBA
- PDF: TBA

## Requirements
The repository depends on some external libraries, which are included as submodules. Please be sure to synchronize the modules with their most up to date commits.
```
git submodule update --init --recursive --remote
```

In order to compile the code in this repository, you need the following:
- A working C++ compiler compliant with the standard C++17.
- A working installation of CMake, at least version 3.23.0. Probably some older versions are also supported.
- (Optional) A working installation of OpenMP that can be found by CMake.

> :warning: In its current state, this repository has been extensively tested only on Windows with the MSVC compiler. No specific Windows/MSVC libraries have been included, so porting to other systems should be free/easy.

## Building instructions
The building process is carried out by CMake. You need to move into a fresh build directory
```sh
mkdir build
cd build
```
and the configure and build with CMake
```sh
cmake ..
cmake --build . --config release --parallel
```

The building process should produce the following tree structures in the build directory
```
build
|
+-- lib
|   |
|   +-- DFY.lib
|
+-- bin
|   |
|   +-- Diskify[.exe]
|   +-- UVPatches[.exe]
|   +-- GraphicalAbstract[.exe]
|   +-- VoronoiRefinementFigure[.exe]
|
+-- (CMake build data)
```

## Applications
The building process produces three executables.
- `Diskify`: the application that computes the cut to parametrize the surface into a topological disk.
- `UVPatches`: the application that computes the segmentation of the surface into topological disks.
- `GraphicalAbstract`: the application that produces the data used for the graphical abstract of the paper.
- `VoronoiRefinementFigure`: the application that produces the data used for the figure in the paper that discusses the refinement steps from the basic Voronoi decomposition to the merged regions.

Notice that, for compatibility with future extensions, the current implementation of the algorithm uses a `std::priority_queue` for computing spanning tree, resulting in the computation of the *maximum* spanning tree, rather than the *minimum* spanning tree, as discussed in the paper. This means that the metrics provided in the implementation are inverted with respect to those discussed in the paper. However, the results are unchanged.

### Diskify
The application can be executed with the following syntax
```sh
./bin/Diskify input_mesh [OPTIONS]
```
where `input_mesh` is the path to a triangular mesh in a format supported by [libigl](https://github.com/libigl/libigl) and the options are as follows.

|Short|Long|Description|Default|
|-|-|-|-|
|`-n`|`--num-samples`|Number of initial samples for the Voronoi decomposition. If <= 0, 2 * (genus + 1) is used.|-1|
|`-d`|`--disk-samples`|Number of samples to use in the subsampling process.|5|
|`-a`|`--algorithm`|The UV parametrization algorithm to use once the cut has been computed. See acceptable values below.|`tutte`|
|`-m`|`--metric`|The metric used to weight the edges of the dual graph and the boundary components. See acceptable values below.|`angular`|
|`-o`|`--output`|The path to the output file.|`input_mesh-uv.obj`|
|`-s`|`--smooth`|Flat or smooth triangle normals in the output file.|**OFF**|
|`-v`|`--verbosity`|Regulated the amount of output text. Values in range `[0, 2]`.|1|
|`-h`|`--help`|Suppress the execution and outputs the usage text.|**OFF**|

Currently, the application implements the following parametrization algorithms (`-a` option):
- `tutte`: Tutte's embedding (no weights, boundary to the `[0, 1]` square)
- `harmonic`: Harmonic mapping (edge lengths as weights, boundary to the `[0, 1]` square)
- `arap`: ARAP deformation of the Tutte's embedding (no boundary constraints)
- `conformal`: least squares conformal embedding (no boundary constraints)

The supported metrics (`-m` option) are the following:
- `euclidean`: edges in the dual graph are weighted by the *Euclidean distance* between the barycenters of the adjacent triangles. Boundary components are weighted by their *total edge length*.
- `geodesic`: edges in the dual graph are weighted by the *geodesic distance* between the barycenters of the adjacent triangles. Boundary components are weighted by their *total edge length*.
- `angular`: edges in the dual graph are weighted by the *dihedral angle* between the adjacent triangles. Boundary components are weighted by the *absolute value of their mean and Gaussian curvatures*.


### UVPatches
The application can be executed with the following syntax
```sh
./bin/UVPatches input_mesh [OPTIONS]
```
where `input_mesh` is the path to a triangular mesh in a format supported by [libigl](https://github.com/libigl/libigl) and the options are as follows.

|Short|Long|Description|Default|
|-|-|-|-|
|`-n`|`--num-samples`|Number of initial samples for the Voronoi decomposition. If <= 0, 2 * (genus + 1) is used.|-1|
|`-d`|`--disk-samples`|Number of samples to use in the subsampling process.|5|
|`-a`|`--algorithm`|The UV parametrization algorithm to use once the cut has been computed. See acceptable values below.|`tutte`|
|`-m`|`--metric`|The metric used to weight the edges of the dual graph and the boundary components. See acceptable values below.|`angular`|
|`-t`|`--threshold`|The threshold value to ignore merging two regions.|0.76|
|`-o`|`--output`|The path to the output file.|`input_mesh-uv.obj`|
|`-s`|`--smooth`|Flat or smooth triangle normals in the output file.|**OFF**|
|`-p`|`--packing`|UV islands are placed so that they do not overlap.|**OFF**|
|`-v`|`--verbosity`|Regulated the amount of output text. Values in range `[0, 2]`.|1|
|`-h`|`--help`|Suppress the execution and outputs the usage text.|**OFF**|

Currently, the application implements the following parametrization algorithms (`-a` option):
- `tutte`: Tutte's embedding (no weights, boundary to the `[0, 1]` square)
- `harmonic`: Harmonic mapping (edge lengths as weights, boundary to the `[0, 1]` square)
- `arap`: ARAP deformation of the Tutte's embedding (no boundary constraints)
- `conformal`: least squares conformal embedding (no boundary constraints)

The supported metrics (`-m` option) are the following:
- `euclidean`: edges in the dual graph are weighted by the *Euclidean distance* between the barycenters of the adjacent triangles. Merging score tries to minimize the *perimeter-to-area ratio*.
- `geodesic`: edges in the dual graph are weighted by the *geodesic distance* between the barycenters of the adjacent triangles. Merging score tries to minimize the *perimeter-to-area ratio*.
- `angular`: edges in the dual graph are weighted by the *dihedral angle* between the adjacent triangles. Merging score tries to minimize the *absolute value of mean and Gaussian curvatures*.


### GraphicalAbstract and VoronoiRefinementFigure
The applications can be executed with the following syntax
```sh
./bin/GraphicalAbstract
./bin/VoronoiRefinementFigure
```
to reproduce the data used for the graphical abstract and the paper's figure discussing the step of the refinement from the basic Voronoi decomposition to the merged regions.