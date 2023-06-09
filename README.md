# Bézier Guarding: Precise Higher-Order Meshing of Curved 2D Domains

Original code + cmake changes to build dependencies in place

## original README.md

------------------------------------------------

### OVERVIEW

This is an example implementation of the Bézier Guarding algorithm described in the article:

Bézier Guarding: Precise Higher-Order Meshing of Curved 2D Domains.
Manish Mandad and Marcel Campen, Osnabrück University.
ACM Transactions on Graphics 39(4), 2020. (Proc. SIGGRAPH 2020)

It furthermore includes various algorithmic preprocessing steps that attempt to verify the validity of the provided input domain curves before the actual Bézier Guarding algorithm is started. This includes a step that, in the case of intersecting input (which violates the input assumptions) attempts to resolve these intersections in the input by splitting the intersecting curves at the approximate points of intersection. Note that these prechecking and preprocessing steps are not guaranteed to succeed in all cases, but they give valuable hints in case some input violates the assumptions.

The code includes the EXACT as well as the FLOAT variant of the algorithm, as explained in the above article. It can be configured as either by setting the ALL_EXACT flag to True or False in CMake.

The input is expected in a simple .json format (see the included example data). In the end the algorithm reports success or failure and optionally (see the command line parameter description below) outputs the generated higher-order mesh in .msh format and/or a visualization in .eps format. Note that these output formats require rounding to standard limited precision floating point numbers, which may cause numerical degeneracies even if the code was run in EXACT mode. While the code does not include a full-fledged mesh optimization stage for the initially generated valid mesh, with the o-flag (see below) a very basic mesh relaxation postprocess can be enabled, increasing the chance that afterwards also complicated results can be written to file without numerical issues.

See the LICENSE file for license information.

Version: 1.1



### HOW TO USE

## Requirements

* CMake (v >= 3.1)
* GMP library
* GMPXX library
* CGAL (only needed for EXACT mode)

## Building

In project root directory:
* mkdir build
* cd build
* cmake .. -DALL_EXACT=True (to build the EXACT version)
* OR: cmake .. -DALL_EXACT=False (to build the FLOAT version)
* make

## Running

Need to specify 1 (optionally 2) command line parameters
  * path to the input file in .json format
  * string with char-flags:
     * 'b' to add bounding box curves
     * 'e' to enable eps output
     * 'g' to enable gmsh output
     * 'o' to enable basic mesh relaxation postprocessing


In directory containing binary "bezmeshCLI", e.g.:
```
./bezmeshCLI "path/to/your/curve_file.json" be
```
