#ifndef BZMSH_BUILD_TRIANGULATION_HPP
#define BZMSH_BUILD_TRIANGULATION_HPP

#include "datastructure/guarded_bezier_curve.hpp"
#include "datastructure/triangulation_constraints.hpp"

namespace bzmsh
{

/**
 * @brief Perform constrained delaunay triangulation given vertices, segment constraints
 *        and segment markers (signalling curved segments). Optionally a bounding
 *        box surrounding the curves is generated before triangulation. Any triangles
 *        not enclosed in either the bbox or closed curves will be removed (curves remain).
 *
 * @param in constraints (vertices, segments and segment markers)
 * @param out resulting (Bezier)-triangle mesh (triangle and vertex lists, no control points yet)
 * @param addBbox whether a bbox should be generated
 * @return true if during CDT no segments had to be merged or split (i.e. num of segments and
 *         vertices unchanged)
 * @return false else (meaning there must have been intersections or ill-defined segments)
 */

bool performCDT(TriangulationConstraints<double>& in, BezierMesh<double>& out, bool addBbox = false);
bool performCDT(TriangulationConstraints<mpq_class>& in, BezierMesh<mpq_class>& out, bool addBbox = false);

/**
 * @brief Construct control points for edges and triangles of a bezier mesh given
 *        the original bezier curves.
 *
 * @param guardedCurves Curves that were used in mesh generation
 * @param mesh mesh generated by performCDT()
 * @return true if successful
 * @return false else
 */
template <typename T>
bool insertMeshCtrlpts(const vector<GuardedBezierCurve<T>>& guardedCurves, BezierMesh<T>& mesh);

/**
 * @brief Uniformly distribute control points inside of triangle defined by 3 vertices
 *
 * @param vertices triangle vertices
 * @param degree triangle degree / order
 * @param ctrlpts resulting control points
 */
template <typename T>
void interpolateInnerCtrlpts(const vector<Vec2<T>>& vertices,
                             int degree,
                             vector<Vec2<T>>& ctrlpts);

/**
 * @brief Construct control points of a given triangle that contains a curve.
 *
 * @param gc guarded bezier curve that is part of the triangle
 * @param tri triangle
 * @param mesh entire mesh
 */
template <typename T>
void calculateInnerCtrlpts(const GuardedBezierCurve<T>& gc, Triangle& tri, BezierMesh<T>& mesh);

/**
 * @brief Split all the triangles of a given mesh (with valid control points!) into 4 subtriangles
 * each.
 *
 * @param mesh triangle mesh containing triangles with properly distributed control points
 * @param subTriangulation resulting triangulation
 */
template <typename T>
void subdivide(const BezierMesh<T>& mesh, BezierMesh<T>& subTriangulation);

} // namespace bzmsh


#endif
