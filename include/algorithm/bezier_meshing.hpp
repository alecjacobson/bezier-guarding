#ifndef BZMSH_BEZIERMESHING_HPP
#define BZMSH_BEZIERMESHING_HPP

#include <string>

namespace bzmsh
{
/**
 * @brief Executes the bezier meshing algorithm on the curves contained
 *        in the specified .json file.
 *
 * Algorithm stages:
 * - Read curves from file
 * - Convert curves to same maximum degree and reparametrize to t=0 -> t=1
 * - Detect irregular curves (quit on encounter)
 * - Detect intersections and try to resolve by splitting at intersections (quit if non-resolvable)
 * - Detect curves with zero-angles (quit on encounter)
 *
 * - Step (1): Split non-guardable curves
 * - Step (2a): Guard each Bezier curve with a quad envelope
 * - Step (2b): Keep splitting at envelope intersections and re-guard until no more intersections
 * - Step (3): Triangulate the space between the quads
 * - Steps (4+5): Compute regular control point positions
 *
 * - Test for Irregularities
 * - Output resulting (higher-order) triangulation
 *
 * @param filename .json file containing curve information
 * @param addBbox whether a box surrounding all curves should be added (so the space between closed
 *        curves gets triangulated as well)
 * @param epsOutput whether a vector graphic .eps visualization of the generated mesh should be
 *        printed
 * @param gmshOutput whether a gmsh format .msh mesh file should be created from the generated mesh
 * @return true if successful
 * @return false if unresolvable error encountered
 */
bool bezierMeshing(std::string filename,
                   bool addBbox = true,
                   bool epsOutput = true,
                   bool gmshOutput = true,
                   bool optimizeMesh = false, int exp_type=0);

} // namespace bzmsh

#endif
